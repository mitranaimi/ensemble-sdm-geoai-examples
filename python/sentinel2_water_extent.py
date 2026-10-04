"""Water-extent mapping from Sentinel-2 L2A with the MNDWI index (public data, no account needed).

Data come from the Microsoft Planetary Computer STAC API (open Sentinel-2 L2A archive).
For two date windows the script builds cloud-masked median composites, computes

    MNDWI = (Green - SWIR1) / (Green + SWIR1)      (Sentinel-2 bands B03 and B11)

thresholds it into a water mask, and derives the two layers that are typically used to
summarise repeated flood or inundation events:
    * maximum extent  = union of the water masks
    * common extent   = intersection of the water masks

Usage (defaults show a seasonal water-body example):
    python python/sentinel2_water_extent.py
    python python/sentinel2_water_extent.py --bbox 45.25 37.35 45.40 37.45 \
        --period-a 2023-04-01/2023-05-31 --period-b 2023-09-15/2023-10-31
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import planetary_computer
import pystac_client
import stackstac

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
# Scene Classification Layer values to mask out: 1 saturated, 3 cloud shadow, 8/9 cloud, 10 cirrus
BAD_SCL = [1, 3, 8, 9, 10]


def item_epsg(item):
    """EPSG code of a STAC item (newer STAC uses 'proj:code', older 'proj:epsg')."""
    props = item.properties
    if "proj:epsg" in props:
        return int(props["proj:epsg"])
    return int(str(props["proj:code"]).split(":")[-1])


def composite(bbox, period, max_cloud=30, resolution=20, max_scenes=8):
    """Cloud-masked median composite of B03, B11 for one date window."""
    catalog = pystac_client.Client.open(STAC_URL, modifier=planetary_computer.sign_inplace)
    items = catalog.search(collections=["sentinel-2-l2a"], bbox=bbox, datetime=period,
                           query={"eo:cloud_cover": {"lt": max_cloud}}).item_collection()
    if len(items) == 0:
        raise SystemExit(f"No Sentinel-2 scenes found for {period} in {bbox}")
    # keep the clearest scenes only: keeps the example fast and memory-light
    items = sorted(items, key=lambda i: i.properties["eo:cloud_cover"])[:max_scenes]
    print(f"{period}: using {len(items)} clearest scenes")
    da = stackstac.stack(items, assets=["B03", "B11", "SCL"], bounds_latlon=bbox,
                         resolution=resolution, epsg=item_epsg(items[0]),
                         dtype="float64", fill_value=np.nan, rescale=False)
    good = ~da.sel(band="SCL").isin(BAD_SCL)
    green = da.sel(band="B03").where(good)
    swir = da.sel(band="B11").where(good)
    return green.median("time").compute(), swir.median("time").compute()


def water_mask(green, swir, threshold=0.0):
    mndwi = (green - swir) / (green + swir)
    return mndwi, (mndwi > threshold)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bbox", nargs=4, type=float, default=[45.25, 37.35, 45.40, 37.45],
                    metavar=("W", "S", "E", "N"))
    ap.add_argument("--period-a", default="2023-04-01/2023-05-31")
    ap.add_argument("--period-b", default="2023-09-15/2023-10-31")
    ap.add_argument("--threshold", type=float, default=0.0, help="MNDWI water threshold")
    ap.add_argument("--out", default="outputs")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ga, sa = composite(args.bbox, args.period_a)
    gb, sb = composite(args.bbox, args.period_b)
    mndwi_a, wa = water_mask(ga, sa, args.threshold)
    mndwi_b, wb = water_mask(gb, sb, args.threshold)

    union = wa | wb          # maximum extent
    common = wa & wb         # common extent
    px_km2 = (20 * 20) / 1e6
    stats = {"water_a_km2": float(wa.sum()) * px_km2, "water_b_km2": float(wb.sum()) * px_km2,
             "maximum_extent_km2": float(union.sum()) * px_km2, "common_extent_km2": float(common.sum()) * px_km2}
    with open(os.path.join(args.out, "water_extent_stats.csv"), "w") as f:
        f.write("metric,value\n" + "\n".join(f"{k},{v:.2f}" for k, v in stats.items()) + "\n")
    print(stats)

    fig, ax = plt.subplots(1, 4, figsize=(16, 4.5))
    for a, arr, title in zip(ax, [wa, wb, union, common],
                             [f"Water {args.period_a}", f"Water {args.period_b}", "Maximum extent (union)", "Common extent (intersection)"]):
        a.imshow(np.asarray(arr), cmap="Blues", interpolation="nearest")
        a.set_title(title, fontsize=9)
        a.axis("off")
    fig.suptitle("Sentinel-2 L2A, MNDWI water masks", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "water_extent.png"), dpi=130)
    print("Saved outputs to", args.out)


if __name__ == "__main__":
    main()
