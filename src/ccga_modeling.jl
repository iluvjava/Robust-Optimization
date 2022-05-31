include("utilities.jl")
include("matrix_construction_export.jl")
using Logging
using .RobustOptim

### Important functions for variables constructions for JuMP models. 


"""
    Does a cartesian outter product on a list of vectors 
    and then return all the products in the form of a tuple packed into 
    a long long list. 
"""
function CartesianOutterProductList(v_list::AbstractVector...)
    result = Iterators.product(v_list...)|>collect
return result[:] end


"""
    Given a variable coefficient holder, using its dimension to convert 
    it to a list of tuples representing all the possible indexing of this 
    multi-dimensional variable. 
        * k: The extra forloop indexing from the CCGA algorithm. 
    
    Returns: 
        All possible indices for the decision variables represented by the 
        instance of the coefficient holder. 
"""
function IndicesList(
    holder::RobustOptim.VariableCoefficientHolder, 
    k::Union{Nothing, Int}=nothing
)
    NdimList = [1:III for III in size(holder)]
    if k === nothing 
        return CartesianOutterProductList(NdimList...)
    end
    push!(NdimList, k:k)
return CartesianOutterProductList(NdimList...) end


abstract type Problem
    # has a JuMP model in it. 
    # MUST IMPLEMENT: GetModel Method. 
end


### ====================================================================================================================
### Master Problem: 
###   * Search for primary feasible solution given branch cut and bounds. 
### ====================================================================================================================
mutable struct MP <: Problem
    M::Model
    w::Array{VariableRef, 3}
    G::Int64; T::Int64
    Tmind::Array{Int}  # Min up down
    Tminu::Array{Int}  # Min up time

    # settings and options stuff: 
    Verbose::Bool

    function MP(M::Model, gamma_upper=1e4)
        this = new(M)
        this.G = RobustOptim.PRIMARY_GENERATORS.generator_count
        this.T = RobustOptim.CONST_PROBLEM_PARAMETERS.HORIZON
        this.Tmind = RobustOptim.PRIMARY_GENERATORS.Tmind
        @assert any(this.T .> this.Tmind) "Tmind: Minimum down time has to be less than "*
        "time horizon"
        this.Tminu = RobustOptim.PRIMARY_GENERATORS.Tminu
        @assert any(this.T .> this.Tminu) "Tmind: Minimum up time has to be less than "*
        "time horizon"
        this.w = 
        @variable(
            this.M, 
            w[
                [:x, :y, :z],  
                1:this.G, 
                1:this.T
            ], Bin
        )
        # Scaler Continuous variables for demand interval. 
        @variable(
            this.M, 
            γ, 
            lower_bound=0, upper_bound=gamma_upper
        )
        # -- copy of global variables: 
        
        AddConstraint2!(this)
        AddConstraint3!(this)
        AddConstraint4!(this)
        AddConstraint5!(this)
        AddObjective!(this)
    return this end

    function MP(gamma_upper=1e4) return MP(Model(HiGHS.Optimizer), gamma_upper) end
    
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
    Introduce feasibility cut for the master problem which comes from the CCGA results.  
        * Continuous decision variable from FSP: u
        * Discrete decision variable from FMP: q
        * Adversarial demands from FMP: d
"""
function IntroduceCut(
    this::MP, 
    u::Vector{JuMP.VariableRef},
    q::Vector{JuMP.VariableRef}, 
    d::Vector{JuMP.VariableRef}
)
    model = this|>GetModel
    u = u.|>value
    q = q.|>value
    d = d.|>value
    w = Getw(this)
    B = RobustOptim.B
    h = RobustOptim.h
    G = RobustOptim.G
    C = RobustOptim.C
    @constraint(model, C*u + G*q .<= H*d + h - B*w)
return this end


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
    * It returns vector as decision variable
"""
function Getw(this::MP)
    Vec = Vector{JuMP.VariableRef}()
    push!(Vec, this.w[1, :, :][:]...)
    push!(Vec, this.w[2, :, :][:]...)
    push!(Vec, this.w[3, :, :][:]...)
return Vec end

"""
    Get the gamma, the objective of the LP function. 
"""
function GetGamma(this::MP)
    model = GetModel(this)
return model[:γ] end


### ====================================================================================================================
# FSP: Lower bound searcher! 
#   * Takes in a demands and it tries to satisfies it by choosing the secondary
#   * discrete and contiuous variables. 
### ====================================================================================================================
"""
    Given demand from FMP, test how feasible it is to determine an Lower Bound 
    for the feasibility slack variables. 
"""
mutable struct FSP <: Problem
    model::JuMP.Model
    u::Vector{JuMP.VariableRef}
    v::Vector{JuMP.VariableRef}
    q::Vector{JuMP.VariableRef}
    
    # Given parameters 
    w::Vector{Float64}
    d_star::Vector{Float64}
    gamma::Number
    
    function FSP()
        @warn("This construction only exists for testing purposes!")
        return FSP(
            zeros(size(RobustOptim.B, 2)),
            10, 
            ones(RobustOptim.d|>length)
        )
    end

    """
        Constrct the FSP, given the parameters setting the primary 
        variables, and then the adversarial demands vector. 
    """
    function FSP(w, gamma, d_star)
        this = new()
        this.model = Model(HiGHS.Optimizer)
        this.d_star = d_star
        this.gamma = gamma
        this.w = w
        this.d_star = d_star
        this|>PrepareVariables!
        this|>AddConstraints!
        @objective(this.model, Min, sum(this.v))
    return this end
    
    """
        Prepare all the variables for the FSP model, decision variables: 
            * u
            * v
            * q
        As a list of flattened JuMP variable refs, it will store the flatten 
        variables to the field of this struct. 
    """
    function PrepareVariables!(this::FSP)
        # Preparing variable: u into JuMP model. 
        u = RobustOptim.u
        @info "Prepareing variables for current FSP model: "
        @variable(this.model, c[u[1]|>IndicesList], lower_bound=0)
        @variable(this.model, c′[u[2]|>IndicesList], lower_bound=0)
        @variable(this.model, p[u[3]|>IndicesList], lower_bound=0)
        @variable(this.model, p′[u[4]|>IndicesList], lower_bound=0)
        @variable(this.model, regu[u[5]|>IndicesList], lower_bound=0)
        @variable(this.model, regu′[u[6]|>IndicesList], lower_bound=0)
        @variable(this.model, regd[u[7]|>IndicesList], lower_bound=0)
        @variable(this.model, regd′[u[8]|>IndicesList], lower_bound=0)
        @variable(this.model, sr[u[9]|>IndicesList], lower_bound=0)
        @variable(this.model, sr′[u[10]|>IndicesList], lower_bound=0)
        @variable(this.model, h[u[11]|>IndicesList], lower_bound=0)
        @variable(this.model, g_plus[u[12]|>IndicesList], lower_bound=0)
        @variable(this.model, g_minus[u[13]|>IndicesList], lower_bound=0)
        @variable(this.model, nsp[u[14]|>IndicesList], lower_bound=0)
        @variable(this.model, nsp′[u[15]|>IndicesList], lower_bound=0)
        this.u = [
            this.model[:c]...,
            this.model[:c′]...,
            this.model[:p]...,
            this.model[:p′]...,
            this.model[:regu]...,
            this.model[:regu′]...,
            this.model[:regd]...,
            this.model[:regd′]...,
            this.model[:sr]...,
            this.model[:sr′]...,
            this.model[:h]...,
            this.model[:g_plus]...,
            this.model[:g_minus]...,
            this.model[:nsp]...,
            this.model[:nsp′]...
        ]
        # Prepare for variable q
        q = RobustOptim.q
        @variable(this.model, x′[q[1]|>IndicesList], binary=true)
        @variable(this.model, y′[q[2]|>IndicesList], binary=true)
        @variable(this.model, z′[q[3]|>IndicesList], binary=true)
        this.q = [this.model[:x′]..., this.model[:y′]..., this.model[:z′]...]

        # Prepare for variable v, the slack. 
        @variable(this.model, v[1:length(RobustOptim.h)], lower_bound=0)
        this.v = this.model[:v]
        
    return end

    function AddConstraints!(this::FSP)
        u = this.u 
        v = this.v
        q = this.q
        h = RobustOptim.h
        d = this.d_star 
        w = this.w
        C = RobustOptim.C
        H = RobustOptim.H
        B = RobustOptim.B
        G = RobustOptim.G
        @info "Preparing constraints for the FSP model. "
        @constraint(this.model, c1, C*u - v .<= H*d + h - B*w - G*q)
    return end


    function Base.getindex(this::FSP, index::Symbol)
    return this.model[index] end
    
end

### FSP methods clusters -----------------------------------------------------------------------------------------------

function Getq(this::FSP)
return this.q end

function Getu(this::FSP)
return this.u end 

function Getv(this::FSP)
return this.v end


# ======================================================================================================================
# FMP, the Upper bound locator. 
#   * Obtains a lower bound by giving a demands that can break the feasibility for all discrete secondary 
#     decision variables. 
# ======================================================================================================================

"""
    Given the primary parameters and the secondary discrete decision variables 
    meshed into a giant feasible set of many many constraints, we are 
    searching for adversarial demands that can break our delivery system. 
        * Determines the upper bound for the feasibility slack. 
"""
mutable struct FMP<:Problem
    w::Vector{Float64}                            # Primary Generator decision variables.               (GIVEN CONSTANT)
    gamma::Float64                                # the bound for the demands.                          (GIVEN CONSTANT)
    q::Vector{Vector{Int}}                        # the secondary discrete decision variables           (GIVEN CONSTANT)
    v::Vector{Vector{JuMP.VariableRef}}           # The slack decision variables for each of the previous demands. 
    u::Vector{Vector{JuMP.VariableRef}}           # The secondary continuous decision variables. 
    k::Int                                        # The iteration number from the ccga. 
    d::Vector{JuMP.VariableRef}                   # The demand decision variable, as a giant vector. 
    eta::JuMP.VariableRef                         # The eta lower bound for all feasibility. 
    lambda::Vector{Vector{JuMP.VariableRef}}      # the dual decision variables. 
    rho_plus::Vector{Vector{JuMP.VariableRef}}    # binary decision variables for the bilinear demands
    rho_minus::Vector{Vector{JuMP.VariableRef}}   # binary decision variables for the bilinear demands
    d_hat::Vector{Float64}                        # the average testing demands vector.                (GIVEN CONSTAINT)

    model::JuMP.Model                             # The JuMP model for a certain instance. 

    function FMP()
        @warn("Construction shouldn't be called except for debugging purposes. ")
        d̂ = 100.0*ones(size(RobustOptim.H, 2))
    return FMP(zeros(size(RobustOptim.B, 2)), 10.0, d̂) end

    function FMP(w::Vector, gamma, d_hat::Vector)
        this = new()
        this.w = w
        this.gamma = gamma
        this.q = Vector{Vector}()
        this.v = Vector{Vector}()
        this.u = Vector{Vector}()
        this.rho_plus = Vector{Vector}()
        this.rho_minus = Vector{Vector}()        
        this.lambda = Vector{Vector}()
        this.k = 1
        this.model = Model(HiGHS.Optimizer)
        this.d_hat = d_hat
        PrepareVariables!(this)
        PrepareConstraints!(this)
        PrepareObjective!(this)
    return this end    

end

# The adding of the constraints API methods for FMP 

"""
    Prepare a new set of decision variables for the given instance of the problem. 
    * u[k]::The continuous decision variable, for the kth iterations of the CCGA Algorithm. 
        * u[k] will be a compositive decision variables for all the modeling decision variables in the original problem.
    * v, λ follows a similar token. 
"""
function PrepareVariables!(this::FMP, q_given::Union{Nothing, Vector{JuMP.VariableRef}})
    k = this.k
    # Preparing variable: u into JuMP model. 
    u = RobustOptim.u
    @info "Prepareing variables for current FMP model. "
    Decisionvariables = Vector()
    push!(Decisionvariables,
        @variable(this.model, [IndicesList(u[1])], lower_bound=0, base_name="c[$(k)]")...)
    push!(Decisionvariables,
        @variable(this.model, [IndicesList(u[2])], lower_bound=0, base_name="c′[$(k)]")...)
    push!(Decisionvariables,
        @variable(this.model, [IndicesList(u[3])], lower_bound=0, base_name="p[$(k)]")...)
    push!(Decisionvariables,
        @variable(this.model, [IndicesList(u[4])], lower_bound=0, base_name="p′[$(k)]")...)
    push!(Decisionvariables,
        @variable(this.model, [IndicesList(u[5])], lower_bound=0, base_name="regu[$(k)]")...)
    push!(Decisionvariables, 
        @variable(this.model, [IndicesList(u[6])], lower_bound=0, base_name="regu′[$(k)]")...)
    push!(Decisionvariables, 
        @variable(this.model, [IndicesList(u[7])], lower_bound=0, base_name="regd[$(k)]")...)
    push!(Decisionvariables, 
        @variable(this.model, [IndicesList(u[8])], lower_bound=0, base_name="regd′[$(k)]")...)
    push!(Decisionvariables, 
        @variable(this.model, [IndicesList(u[9])], lower_bound=0, base_name="sr[$(k)]")...)
    push!(Decisionvariables, 
        @variable(this.model, [IndicesList(u[10])], lower_bound=0, base_name="sr′[$(k)]")...)
    push!(Decisionvariables, 
        @variable(this.model, [IndicesList(u[11])], lower_bound=0, base_name="h[$(k)]")...)
    push!(Decisionvariables, 
        @variable(this.model, [IndicesList(u[12])], lower_bound=0, base_name="g_plus[$(k)]")...)
    push!(Decisionvariables, 
        @variable(this.model, [IndicesList(u[13])], lower_bound=0, base_name="g_minus[$(k)]")...)
    push!(Decisionvariables, 
        @variable(this.model, [IndicesList(u[14])], lower_bound=0, base_name="nsp[$(k)]")...)
    push!(Decisionvariables, 
        @variable(this.model, [IndicesList(u[15])], lower_bound=0, base_name="nsp′[$(k)]")...)
    push!(this.u, Decisionvariables)
    
    # Prepare for variable q, which is going to be a constant. 
    if q_given === nothing
        Decisionvariables = Vector()
        q = RobustOptim.q
        push!(Decisionvariables, rand((0, 1), size(q[1])...)...)
        push!(Decisionvariables, rand((0, 1), size(q[2])...)...)
        push!(Decisionvariables, rand((0, 1), size(q[3])...)...)
        push!(this.q, Decisionvariables)
        this.d = @variable(this.model, d[1:size(RobustOptim.H, 2)])[:]
        this.eta = @variable(this.model, η)
    else
        push!(this.q, (q_given.|>value)[:])
    end

    # Prepare for variable v, the slack. 
    v = @variable(this.model, [1:length(RobustOptim.h)], lower_bound=0, base_name="v[$(k)]")
    push!(this.v, v[:])

    # Parepare the variable lambda, the dual decision variables. 
    λ = @variable(this.model, [1:length(RobustOptim.h)], lower_bound=-1, upper_bound=0, base_name="λ[$(k)]")
    push!(this.lambda, λ[:])

return this end


"""
    Prepare the variables for the r=k th iterations of the CCGA algorithm. 
"""
function PrepareVariables!(this::FMP)
    @assert this.k == 1 "This function can only be called for the first creation of the FMP problem, where initial "*
    "q is not given. "
    PrepareVariables!(this, nothing)
return this end


"""
    Prepare the constraints for the FMP problem, however, it will only add for all the most recent variables with 
    r = k. 
"""
function PrepareConstraints!(this::FMP)
    k = this.k
    η = this.eta
    u = this.u[k]
    v = this.v[k]
    w = this.w
    H = RobustOptim.H
    B = RobustOptim.B 
    G = RobustOptim.G
    C = RobustOptim.C
    h = RobustOptim.h
    λ = this.lambda[k]
    d = this.d
    d̂ = this.d_hat
    γ = this.gamma
    q = this.q[k]
    model = this.model

    @constraint(model, η <= sum(v), base_name="constraint[$(k)][1]")
    @constraint(model, λ'*C .<= 0, base_name="constraint[$(k)][3]")
    @constraint(model, d .<= d̂ .+ γ, base_name="constraint[$(k)][4]")
    @constraint(model, d .>= d̂ .- γ, base_name="constraint[$(k)][5]")

    # Sparse bilinear constraints setup: 
    B̄ = size(H, 2)
    ρ⁺ = @variable(model, [1:B̄], Bin, base_name="ρ⁺[$(k)]")
    push!(this.rho_plus, ρ⁺)
    ρ⁻ = @variable(model, [1:B̄], Bin, base_name="ρ⁻[$(k)]")
    push!(this.rho_minus, ρ⁻)

    IdxTuples = findall(!=(0), H).|>Tuple
    ξ⁺ = @variable(model, [IdxTuples], base_name="ξ⁺[$(k)]")
    ξ⁻ = @variable(model, [IdxTuples], base_name="ξ⁻[$(k)]")
    
    for (j, b) in IdxTuples
        @constraint(model, -ρ⁺[b] <= ξ⁺[(j, b)])
        @constraint(model, ξ⁺[(j, b)] <= ρ⁺[b])
        @constraint(model, λ[j] - (1 - ρ⁺[b])<=ξ⁺[(j, b)])
        @constraint(model, ξ⁺[(j, b)] <= λ[j] + (1 - ρ⁺[b]))
        @constraint(model, -ρ⁻[b] <= ξ⁻[(j, b)])
        @constraint(model, ξ⁻[(j, b)] <= ρ⁻[b])
        @constraint(model, λ[j] - (1 - ρ⁻[b]) <= ξ⁻[(j, b)])
        @constraint(model, ξ⁻[(j, b)] <= λ[j] + (1 - ρ⁻[b]))
    end
    @constraint(model, [b=1:B̄], ρ⁺[b] + ρ⁻[b] == 1, base_name="constraint[$(k)][6]")
    # Binlinear constraints
    @constraint(
        model, 
        dot(λ, H*d̂) + 
        γ*dot(H[H.!=0], ξ⁺[:] .- ξ⁻[:]) + 
        dot(λ,h - B*w - G*q) == sum(v), 
        base_name="constraint[$(k)][7]"
    )
    # new added demand feasibility constraints. 
    @constraint(
        model,
        C*u - v .<= H*(d̂ + γ*ρ⁺ + γ*ρ⁻) + h - B*w - G*q, 
        base_name="constraint[$(k)][2]"
    )  

return this end

"""
    Adversarial demands are coming from vertices of the hyper cube, this will recover 
    the term (d̂ + γ.*ρ⁺ - γ.*ρ⁻), which represents the extreme demands which breaks the system. 
        * It's used for the feasibility cut for the master problem to determine primary decision variable w̄. 
        * Returns the actual value for the FSP lower bound searcher. 
"""
function GetDemandVertex(this::FMP)
    ρ⁺ = value.(this.rho_plus[end])
    ρ⁻ = value.(this.rho_minus[end])
    γ = this.gamma
    d̂ = this.d_hat

return d̂ + γ*ρ⁺ - γ*ρ⁻ end


"""
    Introduce a new secondary binary decision variable as a fixed variable for te system. 
"""
function Introduce!(this::FMP, q::Vector{JuMP.VariableRef})
    this.k += 1
    @assert length(q) == size(RobustOptim.G, 2) "Wrong size for the passed in discrete secondary decision variables: q. "
    PrepareVariables!(this, q)
    PrepareConstraints!(this)
return this end


function PrepareObjective!(this::FMP)
    @objective(this.model, Max, this.eta)
return this end



### ====================================================================================================================
### Generics Functions for abstract type of problems that contain a JuMP inner model, specific overrides too. 
### ====================================================================================================================

"""
    * NaN is returned if there is no solution to the problem, which means it's either infeasible
    or unbounded. 

"""
function objective_value(this::Problem)
    model = GetModel(this)
    status = primal_status(model)
    if status == NO_SOLUTION
        return NaN
    end
return JuMP.objective_value(model) end 

"""
    Solve the inner JuMP optimization problem by directly invoking `optimize!` in JuMP on the model. 
"""
function Solve!(this::Problem)
return optimize!(this|>GetModel) end


"""
    Get the JuMP instance from this model. 
"""
function GetModel(this::Problem) 
return this.model end

function GetModel(this::MP)
return this.M end

"""
    Print out a debug report for the given Problem object. 
"""
function DebugReport(this::Problem)
return end

### ====================================================================================================================
### Trying to experiment with the FSP, FMP instance. 
### ====================================================================================================================





