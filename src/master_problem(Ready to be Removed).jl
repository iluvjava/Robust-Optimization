include("./problem_parameters.jl")
include("./utilities.jl")

abstract type Problem
    # has a JuMP model in it. 
end

mutable struct MP<:Problem
    M::Model
    w::Array{VariableRef, 3}
    G::Int64; T::Int64
    Tmind::Array{Int}  # Min up down
    Tminu::Array{Int}  # Min up time

    # settings and options stuff: 
    Verbose::Bool

    function MP(M::Model)
        this = new(M)
        this.w = 
        @variable(
            this.M, 
            w[
                [:x, :y, :z],  
                1:PRIMARY_GENERATORS.generator_count, 
                1:CONST_PROBLEM_PARAMETERS.HORIZON
            ], Bin
        )
        # Scaler Continuous variables for demand interval. 
        @variable(
            this.M, 
            γ
        )
        # -- copy of global variables: 
        this.G = PRIMARY_GENERATORS.generator_count
        this.T = CONST_PROBLEM_PARAMETERS.HORIZON
        this.Tmind = PRIMARY_GENERATORS.Tmind
        @assert any(this.T .> this.Tmind) "Tmind: Minimum down time has to be less than "*
        "time horizon"
        this.Tminu = PRIMARY_GENERATORS.Tminu
        @assert any(this.T .> this.Tminu) "Tmind: Minimum up time has to be less than "*
        "time horizon"
        AddConstraint2!(this)
        AddConstraint3!(this)
        AddConstraint4!(this)
        AddConstraint5!(this)
    return this end

    function MP() return MP(Model(HiGHS.Optimizer)) end
    
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
    w = (this|>GetModel)[:w]
    for n in 1:this.G, t in this.Tminu[n]: this.T
        @constraint(
            this|>GetModel, 
            sum(
                w[:x, n, tau] for tau in t - this.Tminu[n] + 1: t
            ) 
            <= 
            w[:y, n, t]
        )
    end
return end

function AddConstraint5!(this::MP)
    w = (this|>GetModel)[:w]
    for n in 1:this.G, t in this.Tmind[n]: this.T
        @constraint(
            this|>GetModel, 
            sum(
                w[:z, n, tau] for tau in t - this.Tmind[n] + 1: t
            ) 
            <= 
            1 - w[:y, n, t]
        )
    end
return end

"""
    The initial MP problem requires a bounded feasible solution, we must bound
    the demain interval γ by some large number. 
"""
function AddGammaConstraint!(this::MP)
    γ = (this|>GetModel)[:γ]
    @constraint(
        this|>GetModel, 
        γ <= 1e4  # a large number.  
    )
return end

"""
    creates the objective value for it. 
"""
function AddObjective!(this::MP)
    m = this|>GetModel
    γ = m[:γ]
    @objective(m, Max, γ)
return end

"""
    Get the primary discrete decision variable as a vector of numbers for 
        the CCGA algorithm. 
    It returns vector as decision variable
"""
function Getw(this::MP)
    Vec = Vector{JuMP.VariableRef}()
    push!(Vec, this.w[:x, :, :])
    push!(Vec, this.w[:y, :, :])
    push!(Vec, this.w[:z, :, :])
return Vec end




### Simple Tests ===============================================================

# mp = MP()
# display(mp|>GetModel)