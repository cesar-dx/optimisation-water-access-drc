import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------
# Repo paths
# ---------------------------
ROOT = Path(__file__).resolve().parents[1]

CENTROIDS_CSV = ROOT / "data/demand_grid_kmeans.csv"
JF_CSV = ROOT / "data/JF_functional.csv"
JR_CSV = ROOT / "data/JR_nonfunctional.csv"
JN_CSV = ROOT / "data/JN_candidates.csv"

PASS2_DIR = ROOT / "out"
OUT_DIR = ROOT / "figs/maps"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------
# Budgets to plot
# ---------------------------
BUDGETS_M = [0, 2.5, 5, 7.5, 10, 15, 20, 30]

# ---------------------------
# Helpers
# ---------------------------
def safe_read(path):
    if not path.exists():
        print(f"⚠️ Missing: {path}")
        return pd.DataFrame()
    df = pd.read_csv(path)
    return df if not df.empty else pd.DataFrame()

def load_pass2_ids(kind, budget_m):
    """
    Looks for:
      out/pass2_repairs_b15.csv
      out/pass2_new_b15.csv
    """
    b = int(budget_m) if float(budget_m).is_integer() else str(budget_m).replace(".", "p")
    fname = PASS2_DIR / f"pass2_{kind}_b{b}.csv"
    return safe_read(fname)

# ---------------------------
# Load static data
# ---------------------------
centroids = pd.read_csv(CENTROIDS_CSV)
JF = pd.read_csv(JF_CSV)
JR = pd.read_csv(JR_CSV)
JN = pd.read_csv(JN_CSV)

# Normalize ID columns
for df in (JF, JR, JN):
    df["id"] = df["id"].astype(str)

# Precompute centroid sizes
centroid_sizes = (centroids["w"] / centroids["w"].max()) * 40

# ---------------------------
# Plot loop
# ---------------------------
for b in BUDGETS_M:
    df_rep_ids = load_pass2_ids("repairs", b)
    df_new_ids = load_pass2_ids("new", b)

    if df_rep_ids.empty and df_new_ids.empty:
        print(f"⚠️ Skipping {b}M (no pass2 outputs)")
        continue

    # Merge IDs → coordinates
    repairs = df_rep_ids.merge(JR, on="id", how="left") if not df_rep_ids.empty else pd.DataFrame()
    new = df_new_ids.merge(JN, on="id", how="left") if not df_new_ids.empty else pd.DataFrame()

    plt.figure(figsize=(12, 10))

    # Population centroids
    sc = plt.scatter(
        centroids["lon"], centroids["lat"],
        c=centroids["w"],
        cmap="viridis",
        s=centroid_sizes,
        alpha=0.6,
        label="Population centroids"
    )
    plt.colorbar(sc, label="Population")

    # Functional wells
    plt.scatter(
        JF["lon"], JF["lat"],
        c="blue", s=40, label="Functional wells"
    )

    # Repairs
    if not repairs.empty:
        plt.scatter(
            repairs["lon"], repairs["lat"],
            c="orange", marker="x", s=60, label="Repair wells"
        )

    # New wells
    if not new.empty:
        plt.scatter(
            new["lon"], new["lat"],
            c="green", marker="^", s=60, label="New wells"
        )

    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title(f"Selected Well Locations (Budget = ${b}M)")
    plt.legend()
    plt.tight_layout()

    out = OUT_DIR / f"well_locations_budget_{str(b).replace('.', 'p')}M.png"
    plt.savefig(out, dpi=300)
    plt.close()

    print(f"✅ Saved {out}")