using CSV, DataFrames, Statistics
# scripts/01_neighbor_stats.jl
# Purpose:
# - Compute baseline coverage by functional wells
# - Measure neighbor density and arc counts for a given R_m
# - Sanity-check model size before optimization

# ---------------- Paths & params ----------------
DEMAND_CSV = "data/demand_grid_kmeans.csv"   # id,lat,lon,w
JF_CSV     = "data/JF_functional.csv"     # id,lat,lon
JR_CSV     = "data/JR_nonfunctional.csv"  # id,lat,lon
JN_CSV     = "data/JN_candidates.csv"     # id or nid, lat, lon
R_m        = 1500.0                        # meters

# ---------------- Load data ----------------
demand = CSV.read(DEMAND_CSV, DataFrame)
JF     = CSV.read(JF_CSV, DataFrame)
JR     = CSV.read(JR_CSV, DataFrame)
JN     = CSV.read(JN_CSV, DataFrame)
if :id ∉ names(JN) && :nid ∈ names(JN)
    rename!(JN, :nid => :id)
end

# Ensure numeric
for df in (demand, JF, JR, JN)
    df.lat = Float64.(df.lat)
    df.lon = Float64.(df.lon)
end
demand.w = Float64.(demand.w)

# ---------------- Geo helper ----------------
const EARTH_R = 6_371_000.0
haversine_m(lat1, lon1, lat2, lon2) = begin
    φ1 = deg2rad(lat1); λ1 = deg2rad(lon1)
    φ2 = deg2rad(lat2); λ2 = deg2rad(lon2)
    dφ = φ2 - φ1; dλ = λ2 - λ1
    a = sin(dφ/2)^2 + cos(φ1)*cos(φ2)*sin(dλ/2)^2
    2*EARTH_R*asin(sqrt(a))
end

# ---------------- Neighbor counting (with quick bbox prefilter) ----------------
"""
Return vector n where n[i] = # of wells within R_m of demand i.
Uses a fast lat/lon bounding-box prefilter before haversine.
"""
function count_neighbors(demand::DataFrame, wells::DataFrame, R_m::Float64)
    nI = nrow(demand)
    n  = zeros(Int, nI)
    latW = wells.lat; lonW = wells.lon
    for i in 1:nI
        li = demand.lat[i]; loi = demand.lon[i]
        # ~meters-to-degrees
        dlat = R_m / 111_320.0
        dlon = R_m / (111_320.0 * max(1e-6, cosd(li)))
        # prefilter candidates in bbox
        idx = findall(@. (latW ≥ li - dlat) & (latW ≤ li + dlat) & (lonW ≥ loi - dlon) & (lonW ≤ loi + dlon))
        # exact distance check
        c = 0
        @inbounds for j in idx
            if haversine_m(li, loi, latW[j], lonW[j]) ≤ R_m
                c += 1
            end
        end
        n[i] = c
    end
    return n
end

# Count neighbors per set
nF = count_neighbors(demand, JF, R_m)
nR = count_neighbors(demand, JR, R_m)
nN = count_neighbors(demand, JN, R_m)
nAll = nF .+ nR .+ nN

# ---------------- Baseline coverage by functional wells ----------------
served = nF .> 0
served_cells = count(==(true), served)
served_pop   = sum(demand.w[served])
total_pop    = sum(demand.w)

println("Baseline coverage by functional wells (R_m = $(R_m) m):")
println("  Served demand cells: $(served_cells) / $(nrow(demand)) ($(round(100*served_cells/nrow(demand), digits=2))%)")
println("  Served population: $(round(served_pop; digits=0)) / $(round(total_pop; digits=0)) ($(round(100*served_pop/total_pop, digits=2))%)")

# ---------------- Neighbor stats & arc counts ----------------
avgNF   = mean(nF)
avgNR   = mean(nR)
avgNN   = mean(nN)
avgNAll = mean(nAll)

arcs_F   = sum(nF)      # total arcs from demand to JF within R_m
arcs_R   = sum(nR)      # ... to JR
arcs_N   = sum(nN)      # ... to JN
arcs_all = sum(nAll)    # total arcs used by the model (given R_m)

println("\nNeighbor stats (R_m = $(R_m/1000) km):")
println("  avg |N_F(i)| = $(round(avgNF, digits=3))   arcs_F = $(arcs_F)")
println("  avg |N_R(i)| = $(round(avgNR, digits=3))   arcs_R = $(arcs_R)")
println("  avg |N_N(i)| = $(round(avgNN, digits=3))   arcs_N = $(arcs_N)")
println("  avg |N(i)|   = $(round(avgNAll, digits=3)) arcs_all = $(arcs_all)")

# ---------------- Save served flag (as before) ----------------
demand.served_by_JF = served
CSV.write("data/demand_coverage_baseline.csv", demand)
println("\nWrote coverage file: data/demand_coverage_baseline.csv")