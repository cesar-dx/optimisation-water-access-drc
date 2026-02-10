SSTAR = 3.0675720882843688e7       # <-- plug in Pass-1 optimal covered population
B     = 3e7          # <-- plug in the SAME budget used in Pass-1

using JuMP, Gurobi, CSV, DataFrames, Statistics

# ------------------ Load data ------------------
demand = CSV.read("data/demand_grid_kmeans.csv", DataFrame)   # id, lat, lon, w
JF     = CSV.read("data/JF_functional.csv",     DataFrame)    # id, lat, lon
JR     = CSV.read("data/JR_nonfunctional.csv",  DataFrame)    # id, lat, lon
JN     = CSV.read("data/JN_candidates.csv",     DataFrame)    # id (or nid), lat, lon

# If your JN file uses :nid instead of :id, normalize it:
if :id ∉ names(JN) && :nid ∈ names(JN)
    rename!(JN, :nid => :id)
end

# Ensure numeric columns
for df in (demand, JF, JR, JN)
    df.lat = Float64.(df.lat)
    df.lon = Float64.(df.lon)
end
demand.w = Float64.(demand.w)

# ------------------ Sets & maps ------------------
I      = collect(demand.id)
JF_ids = string.(JF.id)
JR_ids = string.(JR.id)
JN_ids = string.(JN.id)

lat_i = Dict(demand.id .=> demand.lat);  lon_i = Dict(demand.id .=> demand.lon);  w = Dict(demand.id .=> demand.w)
lat_jf = Dict(string.(JF.id) .=> JF.lat); lon_jf = Dict(string.(JF.id) .=> JF.lon)
lat_jr = Dict(string.(JR.id) .=> JR.lat); lon_jr = Dict(string.(JR.id) .=> JR.lon)
lat_jn = Dict(string.(JN.id) .=> JN.lat); lon_jn = Dict(string.(JN.id) .=> JN.lon)

# ------------------ Geo helpers ------------------
const EARTH_R = 6_371_000.0
haversine_m(lat1, lon1, lat2, lon2) = begin
    φ1 = deg2rad(lat1); λ1 = deg2rad(lon1)
    φ2 = deg2rad(lat2); λ2 = deg2rad(lon2)
    dφ = φ2 - φ1; dλ = λ2 - λ1
    a = sin(dφ/2)^2 + cos(φ1)*cos(φ2)*sin(dλ/2)^2
    2*EARTH_R*asin(sqrt(a))
end

# ------------------ Policy / cost params ------------------
R_m    = 1500.0                         # access radius (meters) -- keep same as Pass-1
cR_val = 2200.0                         # repair cost per site (edit if you have per-site)
cN_val = 8000.0                         # new-site cost per site

# ------------------ Neighbor sets (from R_m) ------------------
N_F = Dict{Int, Vector{String}}()
N_R = Dict{Int, Vector{String}}()
N_N = Dict{Int, Vector{String}}()
for i in I
    li, lo = lat_i[i], lon_i[i]
    N_F[i] = [j for j in JF_ids if haversine_m(li, lo, lat_jf[j], lon_jf[j]) <= R_m]
    N_R[i] = [j for j in JR_ids if haversine_m(li, lo, lat_jr[j], lon_jr[j]) <= R_m]
    N_N[i] = [k for k in JN_ids if haversine_m(li, lo, lat_jn[k], lon_jn[k]) <= R_m]
end

# Combine into N_all (unique vector of site ids per demand i)
N_all = Dict{Int, Vector{String}}()
for i in I
    N_all[i] = unique(vcat(N_F[i], N_R[i], N_N[i]))
end

# Precompute distances only for arcs within R_m
d = Dict{Tuple{Int,String},Float64}()
for i in I, j in N_all[i]
    li, lo = lat_i[i], lon_i[i]
    if j in JF_ids
        d[(i,j)] = haversine_m(li, lo, lat_jf[j], lon_jf[j])
    elseif j in JR_ids
        d[(i,j)] = haversine_m(li, lo, lat_jr[j], lon_jr[j])
    else # j in JN_ids
        d[(i,j)] = haversine_m(li, lo, lat_jn[j], lon_jn[j])
    end
end

println("PASS-2 setup: R = $(R_m/1000) km, avg |N_all(i)| = ",
        round(mean([length(N_all[i]) for i in I]), digits=2))

# ------------------ PASS 2 model: Min distance, keep ≥ S* and same budget ------------------
# Build arc list
A = [(i,j) for i in I for j in N_all[i] if haskey(d, (i,j))]

m2 = Model(Gurobi.Optimizer)
set_optimizer_attributes(m2, "OutputFlag"=>1, "MIPGap"=>1e-3, "Threads"=>8)

@variable(m2, x[A] >= 0)               # assignment share of cell i to site j
@variable(m2, z[I], Bin)               # 1 if demand cell i counted as covered
@variable(m2, yR[JR_ids], Bin)         # repair decisions
@variable(m2, yN[JN_ids], Bin)         # new-site decisions

# Each covered cell must be fully assigned across its neighbors
@constraint(m2, [i in I], sum(x[(i,j)] for j in N_all[i]) == z[i])

# Linking: assignments to JR/JN require opening; JF are implicitly open
@constraint(m2, [i in I, j in N_R[i]], x[(i,j)] <= yR[j])
@constraint(m2, [i in I, k in N_N[i]], x[(i,k)] <= yN[k])

# Do not sacrifice Pass-1 coverage
@constraint(m2, sum(w[i]*z[i] for i in I) >= 0.999999 * SSTAR)

# SAME budget as Pass-1
@constraint(m2, sum(cR_val * yR[j] for j in JR_ids) +
               sum(cN_val * yN[k] for k in JN_ids) <= B)

# Objective: minimize population-weighted distance
@objective(m2, Min, sum(w[i] * d[(i,j)] * x[(i,j)] for (i,j) in A))

optimize!(m2)

# ------------------ Report ------------------
obj_pm    = objective_value(m2)                              # person-meters
coveredW  = sum(w[i] * (value(z[i]) >= 0.5 ? 1.0 : 0.0) for i in I)
avg_distm = obj_pm / max(coveredW, 1e-9)

println("\nPASS 2 complete")
println("  Covered population: ", round(coveredW), "  (target S* = ", round(SSTAR), ")")
println("  Avg distance among covered: ", round(avg_distm, digits=2), " m")


CSV.write("out/pass2_repairs.csv", DataFrame(id=[j for j in JR_ids if value(yR[j]) >= 0.5]))
CSV.write("out/pass2_new.csv",     DataFrame(id=[k for k in JN_ids if value(yN[k]) >= 0.5]))