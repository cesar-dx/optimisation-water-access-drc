using CSV, DataFrames, Statistics

# ------------------ FILES ------------------
DEMAND_CSV = "data/demand_grid_kmeans.csv"     # must have: id, lat, lon, w
JF_CSV     = "data/JF_functional.csv"          # must have: id, lat, lon
JR_CSV     = "data/JR_nonfunctional.csv"       # must have: id, lat, lon
JN_CSV     = "data/JN_candidates.csv"          # must have: id (or nid), lat, lon

R_m = 1500.0  # meters

# ------------------ LOAD (STRING COLS) ------------------
demand = CSV.read(DEMAND_CSV, DataFrame)
JF     = CSV.read(JF_CSV, DataFrame)
JR     = CSV.read(JR_CSV, DataFrame)
JN     = CSV.read(JN_CSV, DataFrame)

# If JN uses "nid" instead of "id", rename it (strings!)
if !("id" in names(JN)) && ("nid" in names(JN))
    rename!(JN, "nid" => "id")
end

println("Demand columns: ", names(demand))
println("JF columns: ", names(JF))
println("JR columns: ", names(JR))
println("JN columns: ", names(JN))

# ------------------ ASSERTIONS (STRING) ------------------
@assert all(in.(["id","lat","lon","w"], Ref(names(demand)))) "demand needs id, lat, lon, w"
@assert all(in.(["id","lat","lon"],     Ref(names(JF))))     "JF needs id, lat, lon"
@assert all(in.(["id","lat","lon"],     Ref(names(JR))))     "JR needs id, lat, lon"
@assert all(in.(["id","lat","lon"],     Ref(names(JN))))     "JN needs id, lat, lon"

# ------------------ ENSURE NUMERIC TYPES ------------------
demand[!, "lat"] = Float64.(demand[!, "lat"])
demand[!, "lon"] = Float64.(demand[!, "lon"])
demand[!, "w"]   = Float64.(demand[!, "w"])

for df in (JF, JR, JN)
    df[!, "lat"] = Float64.(df[!, "lat"])
    df[!, "lon"] = Float64.(df[!, "lon"])
end

# ------------------ SETS (keep IDs as strings) ------------------
I      = demand[!, "id"]                         # could be Ints; that’s fine
JF_ids = string.(JF[!, "id"])
JR_ids = string.(JR[!, "id"])
JN_ids = string.(JN[!, "id"])

# index maps for fast coordinate lookup
lat_i = Dict(I .=> demand[!, "lat"])
lon_i = Dict(I .=> demand[!, "lon"])
w_i   = Dict(I .=> demand[!, "w"])

lat_jf = Dict(JF_ids .=> JF[!, "lat"]);  lon_jf = Dict(JF_ids .=> JF[!, "lon"])
lat_jr = Dict(JR_ids .=> JR[!, "lat"]);  lon_jr = Dict(JR_ids .=> JR[!, "lon"])
lat_jn = Dict(JN_ids .=> JN[!, "lat"]);  lon_jn = Dict(JN_ids .=> JN[!, "lon"])

# ------------------ HAVERSINE (meters) ------------------
const EARTH_R = 6_371_000.0
haversine_m(lat1, lon1, lat2, lon2) = begin
    φ1 = deg2rad(lat1); λ1 = deg2rad(lon1)
    φ2 = deg2rad(lat2); λ2 = deg2rad(lon2)
    dφ = φ2 - φ1; dλ = λ2 - λ1
    a = sin(dφ/2)^2 + cos(φ1)*cos(φ2)*sin(dλ/2)^2
    2 * EARTH_R * asin(sqrt(a))
end

# ------------------ NEIGHBOR SETS ------------------
N_F = Dict{eltype(I), Vector{String}}()
N_R = Dict{eltype(I), Vector{String}}()
N_N = Dict{eltype(I), Vector{String}}()

for i in I
    li, lo = lat_i[i], lon_i[i]

    N_F[i] = [j for j in JF_ids if haversine_m(li, lo, lat_jf[j], lon_jf[j]) <= R_m]
    N_R[i] = [j for j in JR_ids if haversine_m(li, lo, lat_jr[j], lon_jr[j]) <= R_m]
    N_N[i] = [k for k in JN_ids if haversine_m(li, lo, lat_jn[k], lon_jn[k]) <= R_m]
end

N_all = Dict{eltype(I), Vector{String}}()
for i in I
    N_all[i] = unique(vcat(N_F[i], N_R[i], N_N[i]))
end

# ------------------ ARC COUNTS + AVG NEIGHBORS ------------------
avgNF = mean(length.(values(N_F)))
avgNR = mean(length.(values(N_R)))
avgNN = mean(length.(values(N_N)))
avgN  = mean(length.(values(N_all)))

arcsF = sum(length.(values(N_F)))
arcsR = sum(length.(values(N_R)))
arcsN = sum(length.(values(N_N)))
arcsA = sum(length.(values(N_all)))

println("\nNeighbor stats (R_m = $(R_m/1000) km):")
println("  avg |N_F(i)| = ", round(avgNF, digits=3), "   arcs_F = ", arcsF)
println("  avg |N_R(i)| = ", round(avgNR, digits=3), "   arcs_R = ", arcsR)
println("  avg |N_N(i)| = ", round(avgNN, digits=3), "   arcs_N = ", arcsN)
println("  avg |N(i)|   = ", round(avgN,  digits=3), "   arcs_all = ", arcsA)


CSV.write("out/debug_neighbor_stats.csv", DataFrame(
     R_m = R_m, avgNF = avgNF, avgNR = avgNR, avgNN = avgNN, avgN = avgN,
     arcsF = arcsF, arcsR = arcsR, arcsN = arcsN, arcsA = arcsA
 ))