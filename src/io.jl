using DataFrames

function normalize_id!(df::DataFrame)
    if :id ∉ names(df) && :nid ∈ names(df)
        rename!(df, :nid => :id)
    end
    return df
end