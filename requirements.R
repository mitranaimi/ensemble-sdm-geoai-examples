# Install the R packages used by R/01_ensemble_sdm.R
install.packages(c("sdm", "terra", "raster", "sp", "usdm",
                   "randomForest", "kernlab", "gbm", "mgcv"))
sdm::installAll()   # installs the modelling backends used by sdm
