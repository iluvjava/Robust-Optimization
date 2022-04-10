include("utilities.jl")

### ============================================================================
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
Base.length(this::Sub2Idx) = prod(this.dims)
Base.IndexStyle(::Sub2Idx) = IndexLinear()
Base.getindex(::Sub2Idx, i::Int) = i



### ============================================================================
mutable struct VariableCoefManager{N} <: AbstractArray{Number, N}
    dims::NTuple{N,Int64}
    ndims::Int64
    v::Symbol           # the name of the variable. 
    coefs::Dict{NTuple{N, Int64}, Number}
    function VariableCoefManager(v::Symbol, dims::Int64...)
        N = dims |>length
        this = new{N}()
        this.dims = dims
        this.ndims = dims |>length
        this.v = v
        this.coefs = Dict{Tuple, Number}()
    return this end
end

Base.size(this::VariableCoefManager) = this.dims
Base.ndims(this::VariableCoefManager) = this.ndims
Base.length(this::VariableCoefManager) = prod(this.dims)
Base.IndexStyle(::VariableCoefManager) = IndexCartesian()

function Base.getindex(this::VariableCoefManager, idx::Int...)
    if idx in keys(this.coefs)
        return this.coefs[idx]
    end
    
return 0 end


function Base.setindex!(this::VariableCoefManager, val::Number, idx::Int...)
    if idx > this.dims
        error("Index out of range. index given is: $(idx), but ndims for variable coef manager is $(this.ndims)")
    end
    this.coefs[idx] = val
end

function Base.empty!(this::VariableCoefManager)
    this.empty!(this.coefs)
end

function Base.show(io::IO, this::VariableCoefManager)
    println(io, "indexer |-> Coefficients")
    for kv in this.coefs
        println(io, "$(kv[1]) |-> $(kv[2])")
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, this::VariableCoefManager)
    show(io, this)
end



### ============================================================================
### Coefficient Matrix Manager
### ============================================================================
"""
    Register variables names, and indices range. Matrix manage 
    allows efficient assignment of coefficients for the variable 
    you defined. 
"""
mutable struct CoefficientMatrixManager
    n::Int64
    m::Int64
    var_col::Dict{Symbol, Int64}    # variable symbol maps to starting column index.
    var_sub::Dict{Symbol, Sub2Idx}   # converter for the variable subscript to linear index. 
    
    # For sparse array conversion. 
    col::Vector{Int64}
    row::Vector{Int64}
    val::Vector{Float64}
    
    function CoefficientMatrixManager()
        this = new()
        this.n = 1
        this.m = 1
        this.var_col = Dict{Symbol,Int64}()
        this.var_sub = Dict{Symbol, Sub2Idx}()

        this.col = Vector{Int64}()
        this.row = Vector{Int64}()
        this.val = Vector{Float64}()
    return this end
end

"""
    Register a variable given it's indexer dimension. 
    For example: x[1:3, 1:5], then register like: 
    RegisterVar(this, :x, 3, 5)

"""
function RegisterVar(
    this::CoefficientMatrixManager, v::Symbol, i::Int64...
)
    if v in keys(this.var_col)
        error("symbol: $(v) already registered starting at column: $(this.var_col[v]), can't register it again")
    end
    idx = Sub2Idx(i...)
    this.var_col[v] = this.n
    this.var_sub[v] = idx
    this.n += length(idx)
return end

"""
    Register the variable, using the an instance of the variable coefficients manager. 
"""
function RegisterVar(this::CoefficientMatrixManager, v::VariableCoefManager)
    RegisterVar(this, v.v, v.dims...)
return end

"""
    Get the indexer, the subscript converter for a given symbol. 
"""
function GetIndexer(this::CoefficientMatrixManager, v::Symbol) return this.var_sub[v] end

function RegisterRow(
        this::CoefficientMatrixManager,
        v::Symbol, 
        d::Dict
    )
    if v in keys(this.var_col)
        colStart = this.var_col[v] - 1
        for (t, n) in d
            vIdx = this.var_sub[v][t...]
            push!(this.col, vIdx + colStart)
            push!(this.row, this.m)
            push!(this.val, n)
        end
        this.n += 1
    return end
    error("variable with symbol: \"$(v)\" never registered, can't add row for this variable yet.")
return end

function RegisterVarableCoefficients(
    this::CoefficientMatrixManager, 
    v::VariableCoefManager
)
    RegisterRow(this, v.v, v.coefs)
    
return end

function GetMatrix(this::CoefficientMatrixManager)
    
return sparse(this.row, this.col, this.val) end


