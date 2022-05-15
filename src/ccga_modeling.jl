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
            1000, 
            2*ones(RobustOptim.d|>length)
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
        # no lowerbound, negative v means more than enough feasibility 
        # for the system. Possitive means infeasible. 
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
        @constraint(this.model, C*u - v .<= H*d + h - B*w - G*q)
    return end

    function ObjVal(this::FSP)
    return objective_value(this.model) end

end

"""
    Given the primary parameters and the secondary discrete decision variables 
    meshed into a giant feasible set of many many constraints, we are 
    searching for adversarial demands that can break our delivery system. 
        * Determines the upper bound for the feasibility slack. 
"""
mutable struct FMP
    w::Vector
    gamma::Number
    q::Vector

    function FMP(w, gamma, q)
        this.w = w
        this.gamma = gamma
        this.q = q
    return end
end

this = FSP()
optimize!(this.model)

