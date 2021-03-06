using SparseArrays, LinearAlgebra

### VariableCoefficientHolder ------------------------------------------------------------------------------------------
###     A class that holds the coefficients using the subscripts of a variable
###     * Subscripts can be tensor indexer. 
mutable struct VariableCoefficientHolder{N} <: AbstractArray{Number, N}
    dims::NTuple{N,Int64}
    ndims::Int64
    v::Symbol           # the name of the variable. 
    coefs::Dict{Int, Number}
    function VariableCoefficientHolder(v::Symbol, dims::Int64...)
        N = dims |>length
        this = new{N}()
        this.dims = dims
        this.ndims = dims |>length
        this.v = v
        this.coefs = Dict{Int, Number}()
    return this end

end

Base.size(this::VariableCoefficientHolder) = this.dims
Base.ndims(this::VariableCoefficientHolder) = this.ndims
Base.length(this::VariableCoefficientHolder) = prod(this.dims)
Base.IndexStyle(::VariableCoefficientHolder) = IndexLinear()

function Base.getindex(this::VariableCoefficientHolder, idx::Int)
    if idx in keys(this.coefs)
        return this.coefs[idx]
    end
return 0 end

function Base.setindex!(
    this::VariableCoefficientHolder, 
    val::Number, 
    idx::Int
)
    if idx > this|>length
        error("Index out of range. index given is: $(idx)"*
        ", but ndims for variable coef manager is $(this.ndims)")
    end
    this.coefs[idx] = val
return end

function Base.empty!(this::VariableCoefficientHolder)
    empty!(this.coefs)
return end

function (this::VariableCoefficientHolder)()
    empty!(this)
return this end


### COEFFICIENT MATRIX -------------------------------------------------------------------------------------------------
###     Keep track of the coefficient of the variables and the symbol 
###     representing the 
###     variable for each individual constraints

mutable struct CoefficientMatrix
    n::Int64        # number of columns in the matrix. 
    m::Int64        # bumber of rows in the matrix. 
    var_posi::Dict{Symbol, Int64}   # which column is the variable assigned to.
    var_dims::Vector{Pair}
    
    coefs::Vector{NTuple{3, Float64}}
    function CoefficientMatrix()
        this = new()
        this.n = 0
        this.m = 1
        this.var_posi = Dict{Symbol, Int64}()
        this.coefs = Vector{Tuple}()
        this.var_dims = Vector{Pair}()
    return this end
end

function Register(this::CoefficientMatrix, v::VariableCoefficientHolder)
    if v.v in this.var_posi |> keys
        error("Symbol: $(v.v) already register to the matrix, can't do it again")
    end
    this.var_posi[v.v] = this.n
    this.n += v|>length
    push!(this.var_dims, v.v => v.dims)
return end

function (this::CoefficientMatrix)(vh::VariableCoefficientHolder)
    if !(vh in this)
        Register(this, vh)
    end
    AddCoefficientsFor(this, vh)
    empty!(vh)
return this end

function (this::CoefficientMatrix)(vh::VariableCoefficientHolder...)
    for II in vh
        this(II)
    end
return this end

function AddCoefficientsFor(this::CoefficientMatrix, var::VariableCoefficientHolder)
    if var.v in this.var_posi|>keys
        for (k, v) in var.coefs
            push!(this.coefs, (this.m, k + this.var_posi[var.v], v))
        end
        return
    end
    error("Variable with symbol: \"$(var.v)\" never registered. ")
end

function Base.in(var::VariableCoefficientHolder, this::CoefficientMatrix)
return var.v in this.var_posi|>keys end

function (this::CoefficientMatrix)()
    this|>NextRow!
return this end

function NextRow!(this::CoefficientMatrix)
    this.m += 1
return end

function JumpToRow(this::CoefficientMatrix, row::Int)
    @assert row > 0 "Expect jumping to a row that is strictly positive"*
    ", but got: $(row) <= 0"
    this.m = row
end

function GetMatrix(this::CoefficientMatrix)
    if this.m == 1
        error("Can't get matrix yet, cursor is on the first row, you need to use 'C()'"*
        " to complete the first row of the coefficient matrix. ")
    end
return sparse(
    [convert(Float64, e[1]) for e in this.coefs], 
    [convert(Float64, e[2]) for e in this.coefs], 
    [convert(Float64, e[3]) for e in this.coefs], 
    this.m - 1, 
    this.n
) end

function VariableList(this::CoefficientMatrix)
return [s[1] for s in this.var_dims] end


"""
    Given a list of matrices, sync the cursor to the maximum row number on all 
    of these matrices. 

    Returns nothing. 
"""
function SyncRow(matrices::CoefficientMatrix...)
    MaxRow = maximum([m.m for m in matrices])
    for m in matrices
        JumpToRow(m, MaxRow)
    end
return MaxRow end
