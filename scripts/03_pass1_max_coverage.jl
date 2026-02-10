using JuMP, Gurobi, CSV, DataFrames

# ------------------ Inputs ------------------
DEMAND_CSV = "data/demand_grid_kmeans.csv"
JF_CSV     = "data/JF_functional.csv"
JR_CSV     = "data/JR_nonfunctional.csv"
JN_CSV     = "data/JN_candidates.csv"

R_m = 1500.0
B   = 30_000_000.0

cR = 2200.0
cN = 8000.0

# ------------------ Load data ------------------
demand = CSV.read(DEMAND_CSV, DataFrame)
JF     = CSV.read(JF_CSV, DataFrame)
JR     = CSV.read(JR_CSV, DataFrame)
JN     = CSV.read(JN_CSV, DataFrame)

# normalize id column for candidates
if :id ∉ names(JN) && :nid ∈ names(JN)
    rename!(JN, :nid => :id)
end

# ensure numeric columns
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

lat_i = Dict(demand.id .=> demand.lat)
lon_i = Dict(demand.id .=> demand.lon)
w     = Dict(demand.id .=> demand.w)

lat_jf = Dict(string.(JF.id) .=> JF.lat); lon_jf = Dict(string.(JF.id) .=> JF.lon)
lat_jr = Dict(string.(JR.id) .=> JR.lat); lon_jr = Dict(string.(JR.id) .=> JR.lon)
lat_jn = Dict(string.(JN.id) .=> JN.lat); lon_jn = Dict(string.(JN.id) .=> JN.lon)

# ------------------ Geo helper ------------------
const EARTH_R = 6_371_000.0
haversine_m(lat1, lon1, lat2, lon2) = begin
    φ1 = deg2rad(lat1); λ1 = deg2rad(lon1)
    φ2 = deg2rad(lat2); λ2 = deg2rad(lon2)
    dφ = φ2 - φ1; dλ = λ2 - λ1
    a = sin(dφ/2)^2 + cos(φ1)*cos(φ2)*sin(dλ/2)^2
    2*EARTH_R*asin(sqrt(a))
end

# ------------------ Neighbor sets within R_m ------------------
N_R = Dict{Int, Vector{String}}()
N_N = Dict{Int, Vector{String}}()
f_cov = Dict{Int, Int}()

for i in I
    li, lo = lat_i[i], lon_i[i]

    # functional neighbors count
    f_cov[i] = count(j -> haversine_m(li, lo, lat_jf[j], lon_jf[j]) <= R_m, JF_ids)

    # candidate neighbors
    N_R[i] = [j for j in JR_ids if haversine_m(li, lo, lat_jr[j], lon_jr[j]) <= R_m]
    N_N[i] = [k for k in JN_ids if haversine_m(li, lo, lat_jn[k], lon_jn[k]) <= R_m]
end

# ------------------ Pass 1 model ------------------
m1 = Model(Gurobi.Optimizer)
set_optimizer_attributes(m1, "OutputFlag"=>1, "MIPGap"=>1e-3, "Threads"=>8)

@variable(m1, z[I], Bin)          # covered demand cell?
@variable(m1, yR[JR_ids], Bin)    # repair this nonfunctional well?
@variable(m1, yN[JN_ids], Bin)    # build new well at candidate site?

# Coverage logic (correct)
@constraint(m1, [i in I],
    z[i] <= (f_cov[i] > 0 ? 1 : 0) + sum(yR[j] for j in N_R[i]) + sum(yN[k] for k in N_N[i])
)

# Budget constraint
@constraint(m1, sum(cR * yR[j] for j in JR_ids) + sum(cN * yN[k] for k in JN_ids) <= B)

# Objective: maximize population covered
@objective(m1, Max, sum(w[i] * z[i] for i in I))

optimize!(m1)

Sstar = objective_value(m1)
println("\nPASS 1 complete")
println("  Budget B = $B")
println("  Covered population S* = ", round(Sstar))
println("  Repairs chosen = ", count(j -> value(yR[j]) > 0.5, JR_ids))
println("  New chosen = ", count(k -> value(yN[k]) > 0.5, JN_ids))