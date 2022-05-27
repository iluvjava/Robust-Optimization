include("utilities.jl")
include("matrix_construction_export.jl")
include("master_problem.jl")
using Logging

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



### ============================================================================
# FSP: Lower bound searcher! 
### ============================================================================
"""
    Given demand from FMP, test how feasible it is to determine an Lower Bound 
    for the feasibility slack variables. 
"""
mutable struct FSP
    model::JuMP.Model
    u::Vector
    v::Vector
    q::Vector
    
    # Given parameters 
    w::Vector
    d_star::Vector
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

### FSP methods clusters -------------------------------------------------------

function ObjVal(this::FSP)
return objective_value(this.model) end


function Getq(this::FSP)
return value.(this.q) end


function Getu(this::FSP)
return value.(this.u) end 


function Solve!(this::FSP)
    optimize!(this.model)
return objective_value(this.model) end

"""
    Get the model from the instance of FSP. 
"""
function JuMPModel(this::FSP) return this.model end


"""
    Given the primary parameters and the secondary discrete decision variables 
    meshed into a giant feasible set of many many constraints, we are 
    searching for adversarial demands that can break our delivery system. 
        * Determines the upper bound for the feasibility slack. 
"""
mutable struct FMP
    w::Vector                   # Primary Generator decision variables. 
    gamma::Number               # the bound for the demands. 
    q::Vector{Vector}           # the secondary discrete decision variables (GIVEN CONSTANT).  
    v::Vector{Vector}           # The slack decision variables for each of the previous demands. 
    u::Vector{Vector}           # The secondary continuous decision variables. 
    k::Number                   # The iteration number from the ccga. 
    eta::JuMP.VariableRef       # The eta lower bound for all feasibility. 
    lambda::Vector{Vector}

    model::JuMP.Model            # The JuMP model for a certain instance. 

    function FMP()
        @warn("Construction shouldn't be called except for debugging purposes. ")
    return FMP(zeros(size(RobustOptim.B, 2)), 0) end

    function FMP(w, gamma)
        this = new()
        this.w = w
        this.gamma = gamma
        this.q = Vector{Vector}()
        this.v = Vector{Vector}()
        this.u = Vector{Vector}()
        this.lambda = Vector{Vector}()
        this.k = 1
        this.model = Model(HiGHS.Optimizer) 

        PrepareVariables!(this)
    return this end
    


    

end

# ======================================================================================================================
# The adding of the constraints API methods for FMP 
# ======================================================================================================================

"""
    Prepare a new set of decision variables for the given instance of the problem. 
    * u[k]::The continuous decision variable, for the kth iterations of the CCGA Algorithm. 
        * u[k]  will be a compositive decision variables for all the modeling decision variables in the original problem.
    * v, λ follows a similar token. 
"""
function PrepareVariables!(this::FMP, q_given::Union{Nothing, JuMP.VariableRef})
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
    else
        push!(this.q, (q_given.|>value)[:])
    end

    # Prepare for variable v, the slack. 
    v = @variable(this.model, [1:length(RobustOptim.h)], lower_bound=0, base_name="v[$(k)]")
    push!(this.v, v[:])

    # Parepare the variable lambda, the dual decision variables 
    λ = @variable(this.model, [1:length(RobustOptim.h)], lower_bound=-1, upper_bound=0, base_name="λ[$(k)]")
    push!(this.lambda, λ[:])
    
    # prepare the variable η
    @variable(this.model, η)
    this.k += 1

return this end


function PrepareVariables!(this::FMP)
    @assert this.k == 1 "This function can only be called for the first creation of the FMP problem, where initial "*
    "q is not given. "
    PrepareVariables!(this, nothing)
return this end


"""
    Prepare the constraints for the FMP problem. 
"""
function PrepareConstraints!(this::FMP)
    for r in 1:this.k

    end
return this end

function Introduce!(this::FMP, q::JuMP.VariableRef)


return this end


function PrepareObjective!(this::FMP)

return this end



### ============================================================================
### Trying to experiment with the FSP instance. 
### ============================================================================

# this = FSP()
# Solve!(this)
# q = Getq(this)
# u = getu(this)

# open("constraint_print.txt", "w") do file
#     for II = 1:length(this[:c1])
#         write(file, this[:c1][II]|>repr)
#         write(file, "\n")
#     end
# end
this = FMP()