mutable struct Sub2Idx{N} <: AbstractArray{UInt64, N}
    dims::NTuple{N,Int64}
    ndims::Int64
    function Sub2Idx(dims::Int64...)
        this = new{dims|>length}()
        this.dims = dims
        this.ndims = dims |>length
    return this end
end

Base.size(this::Sub2Idx) = this.dims
Base.ndims(this::Sub2Idx) = this.ndims
Base.IndexStyle(::Sub2Idx) = IndexLinear()
Base.getindex(::Sub2Idx, i::Int) = i


mutable struct 
