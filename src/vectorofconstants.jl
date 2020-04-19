struct VectorOfConstants{T}
    val::T
    len::Int64
end
Base.length(v::VectorOfConstants)=v.len
Base.size(v::VectorOfConstants)=(v.len,)
Base.getindex(v::VectorOfConstants{T},i) where T =v.val
