using CSV, DataFrames

# ---------- Inputs ----------
demand = CSV.read("data/demand_grid_kmeans.csv", DataFrame)    # id, lat, lon, w
JF     = CSV.read("data/JF_functional.csv", DataFrame)      # id, lat, lon

# ---------- Geo helpers ----------
const EARTH_R = 6_371_000.0
haversine_m(lat1, lon1, lat2, lon2) = begin
    φ1 = deg2rad(lat1); λ1 = deg2rad(lon1)
    φ2 = deg2rad(lat2); λ2 = deg2rad(lon2)
    dφ = φ2 - φ1; dλ = λ2 - λ1
    a = sin(dφ/2)^2 + cos(φ1)*cos(φ2)*sin(dλ/2)^2
    2*EARTH_R*asin(sqrt(a))
end

# ---------- Params ----------
R_m   = 1500     # coverage radius
K     = 3000       # how many new-site candidates to keep
Sep_m = 1.5*R_m     # min spacing between candidates (greedy thinning)

# ---------- Uncovered test: NOT within R of any JF or JR ----------
function is_uncovered(lat, lon, JF::DataFrame, R_m)
    for r in eachrow(JF)
        if haversine_m(lat, lon, r.lat, r.lon) <= R_m
            return false
        end
    end
    return true
end

demand.:uncovered = map(r -> is_uncovered(r.lat, r.lon, JF, R_m), eachrow(demand))

# ---------- Sort uncovered by population and enforce spacing ----------
cand = demand[demand.uncovered .== true, [:id,:lat,:lon,:w]]
sort!(cand, :w, rev=true)

selected = DataFrame(id=Int[], lat=Float64[], lon=Float64[], w=Float64[])
for r in eachrow(cand)
    good = true
    for s in eachrow(selected)
        if haversine_m(r.lat, r.lon, s.lat, s.lon) < Sep_m
            good = false; break
        end
    end
    if good
        push!(selected, (r.id, r.lat, r.lon, r.w))
        if nrow(selected) >= K; break; end
    end
end

rename!(selected, :id => :nid)   # optional but nice: distinguish new-site ids
CSV.write("data/JN_candidates.csv", selected[:, [:nid, :lat, :lon]])
println("Wrote ", nrow(selected), " candidates to data/JN_candidates.csv")