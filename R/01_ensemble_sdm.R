# Ensemble species-distribution (habitat-suitability) workflow with the sdm package.
#
# Uses ONLY the example data that ships with the sdm package (public), so the script is
# fully reproducible. The structure mirrors a production workflow:
#   presence points + background points -> collinearity check -> several algorithms
#   -> bootstrap / cross-validation replicates -> evaluation -> variable importance
#   -> AUC-weighted ensemble -> tile-wise prediction of a raster.
#
# Run from the repository root:   Rscript R/01_ensemble_sdm.R

suppressPackageStartupMessages({
  library(sdm)
  library(terra)
  library(raster)
  library(usdm)
})
source("R/helpers.R")

set.seed(42)
dir.create("outputs", showWarnings = FALSE)

## 1. Data (shipped with the sdm package) ------------------------------------
sp_file   <- system.file("external/species.shp", package = "sdm")
pts       <- terra::vect(sp_file)
pres      <- pts[pts$Occurrence == 1, ]                       # keep presences only
env_files <- list.files(system.file("external", package = "sdm"),
                        pattern = "asc$", full.names = TRUE)
env       <- terra::rast(env_files)
message("Presence points: ", nrow(pres), " | predictors: ", paste(names(env), collapse = ", "))

## 2. Collinearity screening (variance inflation factor) ----------------------
env_df  <- na.omit(as.data.frame(env))
vif_tab <- usdm::vif(env_df)                       # VIF of all predictors
print(vif_tab)
write.csv(vif_tab, "outputs/vif.csv", row.names = FALSE)

vs   <- usdm::vifstep(env_df, th = 10)             # drop predictors with VIF > 10, one at a time
keep <- vs@results$Variables
message("Predictors kept after VIF screening: ", paste(keep, collapse = ", "))
env  <- env[[keep]]

## 3. Data object: presence-only points + 1000 background points --------------
pres$Occurrence <- 1
d <- sdmData(Occurrence ~ ., train = as(pres, "Spatial"),
             predictors = raster::stack(env),
             bg = list(n = 1000, method = "gRandom"))

## 4. Fit several algorithms with repeated resampling -------------------------
m <- sdm(Occurrence ~ ., data = d,
         methods = c("glm", "rf", "svm"),
         replication = c("boot", "cv"), cv.folds = 5, test.percent = 30, n = 1)
print(m)

## 5. Evaluation and variable importance --------------------------------------
ev <- getEvaluation(m, stat = c("AUC", "TSS", "COR"), opt = 2)
write.csv(ev, "outputs/evaluation.csv", row.names = FALSE)
ev_m <- merge(ev, getModelInfo(m)[, c("modelID", "method")], by = "modelID")
ev_summary <- aggregate(cbind(AUC, TSS, COR) ~ method, data = ev_m, FUN = mean)
print(ev_summary)                                    # mean performance per algorithm
write.csv(ev_summary, "outputs/evaluation_by_method.csv", row.names = FALSE)

vi <- getVarImp(m)
png("outputs/variable_importance.png", width = 900, height = 600, res = 120)
print(plot(vi, "auc"))
dev.off()

## 6. AUC-weighted ensemble over the whole raster -----------------------------
en <- ensemble(m, raster::stack(env),
               setting = list(method = "weighted", stat = "AUC"))

## 7. Same ensemble, predicted tile by tile and mosaicked ---------------------
## (identical result; keeps memory flat when the raster has tens of millions of cells)
en_tiled <- predict_ensemble_tiled(m, env, n_row = 2, n_col = 2, weight = "AUC")
en_tiled <- terra::resample(en_tiled, terra::rast(en))
agree <- stats::cor(terra::values(terra::rast(en)), terra::values(en_tiled),
                    use = "complete.obs")
message("Correlation between direct and tile-wise ensemble: ", round(agree, 4))

png("outputs/ensemble_suitability.png", width = 900, height = 700, res = 120)
terra::plot(terra::rast(en), main = "AUC-weighted ensemble: habitat suitability")
terra::points(pres, pch = 16, cex = 0.5)
dev.off()

writeLines(capture.output(sessionInfo()), "outputs/sessionInfo.txt")
message("Done. Results are in outputs/")
