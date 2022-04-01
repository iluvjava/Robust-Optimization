include("./problem_parameters.jl")
include("./utilities.jl")


mutable struct MP
    M::Model
    w::Array{VariableRef, 3}
    G::Int64; T::Int64

    function MP(M::Model)
        this = new(M)
        this.w = 
        @variable(
            this.M, 
            w[
                [:x, :y, :z],  
                1:CONST_PROBLEM_PARAMETERS.GENERATOR_PRIMARY, 
                1:CONST_PROBLEM_PARAMETERS.HORIZON
            ], Bin
        )
        # -- instnace copy of global variable: 
        this.G = CONST_PROBLEM_PARAMETERS.GENERATOR_PRIMARY
        this.T = CONST_PROBLEM_PARAMETERS.HORIZON
        AddConstraint2!(this)
        AddConstraint3!(this)
    return this end
    function MP() return MP(Model(GLPK.Optimizer)) end 
    
end

function AddConstraint2!(this::MP)
    w = (this|>GetModel)[:w]
    for n in 1: this.G, t in 1: this.T
        @constraint(
            this|>GetModel, w[:x, n, t] + w[:z, n, t] <= 1
        )
    end
return end

function AddConstraint3!(this::MP)
    w = (this|>GetModel)[:w]
    for n in 1: this.G, t in 2: this.T
        @constraint(
            this|>GetModel,
            w[:y, n, t] - w[:y, n, t - 1] == w[:x, n, t] - w[:z, n, t]
        )
    end
return end

function AddConstraint4!(this::MP)

return end

function GetModel(this::MP) return this.M end


### Simple Tests ===============================================================

mp = MP()
display(mp|>GetModel)