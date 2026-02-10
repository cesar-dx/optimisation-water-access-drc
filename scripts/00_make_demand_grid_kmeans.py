# scripts/00_make_demand_grid_kmeans.py
# K-MEANS CLUSTERING FOR POPULATION-DENSITY-ADAPTIVE DEMAND GRID (FAST AGGREGATION)
# Outputs: CSV with columns id, lat, lon, w (same schema as before)

from sklearn.cluster import MiniBatchKMeans
import numpy as np
import pandas as pd
from osgeo import gdal, ogr


def create_kmeans_demand_grid(tif_path, poly_path, out_csv, n_clusters=100000):
    """
    Create demand grid using K-means clustering on population-weighted pixels.

    - Uses sample_weight (no point replication).
    - Aggregates populations per cluster WITHOUT Python looping over all k.
    - Output schema: id, lat, lon, w  (w = cluster population)

    Args:
        tif_path: Path to WorldPop raster (e.g., .tif)
        poly_path: Path to polygon boundary (GeoJSON/Shapefile)
        out_csv: Output CSV path
        n_clusters: Number of clusters / demand cells

    Returns:
        DataFrame with columns: id, lat, lon, w
    """
    print(f"Creating {n_clusters:,} demand cells using K-means clustering...")

    # ------------------ 1) Load raster ------------------
    ds = gdal.Open(tif_path)
    if ds is None:
        raise RuntimeError(f"Cannot open raster: {tif_path}")

    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(float)
    nodata = band.GetNoDataValue()
    gt = ds.GetGeoTransform()
    height, width = arr.shape

    if nodata is not None:
        arr[arr == nodata] = np.nan

    # ------------------ 2) Rasterize polygon mask ------------------
    v_ds = ogr.Open(poly_path)
    if v_ds is None:
        raise RuntimeError(f"Cannot open polygon: {poly_path}")

    layer = v_ds.GetLayer()
    mem_drv = gdal.GetDriverByName("MEM")
    mask_ds = mem_drv.Create("", width, height, 1, gdal.GDT_Byte)
    mask_ds.SetGeoTransform(gt)
    mask_ds.SetProjection(ds.GetProjection())
    gdal.RasterizeLayer(mask_ds, [1], layer, burn_values=[1])
    poly_mask = mask_ds.GetRasterBand(1).ReadAsArray().astype(bool)

    # ------------------ 3) Extract populated pixels ------------------
    valid_mask = poly_mask & np.isfinite(arr) & (arr > 0)
    row_idx, col_idx = np.where(valid_mask)
    populations = arr[valid_mask].astype(float)

    # Convert pixel indices -> lon/lat (centroid of pixel)
    x0, px, rx, y0, ry, py = gt
    lons = x0 + (col_idx + 0.5) * px + (row_idx + 0.5) * rx
    lats = y0 + (col_idx + 0.5) * ry + (row_idx + 0.5) * py

    total_pop = populations.sum()
    print(f"Extracted {len(populations):,} populated pixels (total pop: {total_pop:,.0f})")

    coords = np.column_stack([lons, lats])
    weights = populations  # sample_weight

    # ------------------ 4) Fit MiniBatchKMeans ------------------
    print(f"Running MiniBatchKMeans with {n_clusters:,} clusters (sample_weight enabled)...")
    kmeans = MiniBatchKMeans(
        n_clusters=n_clusters,
        batch_size=50000,
        max_iter=50,
        random_state=42,
        verbose=1,
        n_init=3
    )
    kmeans.fit(coords, sample_weight=weights)

    # Assign points to clusters
    labels = kmeans.predict(coords)

    # ------------------ 5) FAST aggregation (NO Python loop over k) ------------------
    # Sum population per cluster id (0..n_clusters-1)
    pop_by_cluster = np.bincount(labels, weights=populations, minlength=n_clusters)

    # Keep only non-empty clusters (some can be empty)
    nonempty = pop_by_cluster > 0
    nonempty_idx = np.where(nonempty)[0]

    # Centroids from kmeans
    centers = kmeans.cluster_centers_[nonempty_idx]  # shape: (#nonempty, 2)

    df_out = pd.DataFrame({
        "id": nonempty_idx + 1,           # 1-indexed ids
        "lat": centers[:, 1],
        "lon": centers[:, 0],
        "w": pop_by_cluster[nonempty_idx].astype(float),
    })

    # ------------------ 6) Save ------------------
    import os
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df_out.to_csv(out_csv, index=False)

    print("\n✅ K-means clustering complete!")
    print(f"   Wrote {len(df_out):,} demand cells -> {out_csv}")
    print(f"   Total population captured: {df_out['w'].sum():,.0f}")
    print(f"   Non-empty clusters: {len(df_out)} / {n_clusters}")

    return df_out


if __name__ == "__main__":
    # Example usage
    kmeans_df = create_kmeans_demand_grid(
        tif_path="cod_ppp_2020_1km_Aggregated_UNadj.tif",
        poly_path="geoboundaries.geojson",
        out_csv="data/demand_grid_kmeans.csv",
        n_clusters=100000
    )

    print("\n📊 Summary:")
    print(f"   Mean population per cell: {kmeans_df['w'].mean():,.0f}")
    print(f"   Median population per cell: {kmeans_df['w'].median():,.0f}")
    print(f"   Min/Max population: {kmeans_df['w'].min():,.0f} / {kmeans_df['w'].max():,.0f}")