"""
$(TYPEDEF)

Vector with constant value
"""
struct VectorOfConstants{T}
    val::T
    len::Int64
end

"""
$(TYPEDSIGNATURES)

Length
"""
Base.length(v::VectorOfConstants)=v.len

"""
$(TYPEDSIGNATURES)

Size
"""
Base.size(v::VectorOfConstants)=(v.len,)

"""
$(TYPEDSIGNATURES)

Access
"""
Base.getindex(v::VectorOfConstants,i) = v.val

"""
$(TYPEDSIGNATURES)

Return vector of unique values
"""
Base.unique(v::VectorOfConstants)  = [v.val]

