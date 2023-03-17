### ====================================================================================================================
### Important functions for variables constructions for JuMP models.
### ====================================================================================================================


"""
Does a cartesian outter product on a list of vectors
and then return all the products in the form of a tuple packed into
a long long list.
### Arguments
- `v_listAbstractVector`: A list of ordered elements, the output is a tuple where each 
element from the tuple is from one of the list in the same order, constituting 
the cartesian product of all the vectors.
"""
function CartesianOutterProductList(v_list::AbstractVector...)
    result = Iterators.product(v_list...)|>collect
return result[:] end



"""
Given a variable coefficient holder, using its dimension to convert
it to a list of tuples representing all the possible indexing of this
multi-dimensional variable.

### Argument
- `holder::MatrixConstruct.VariableCoefficientHolder`: The variable holder that is a container 
for a variable indexed by a lot of indices. 
- `k::Union{Nothing, Int}=nothing`: An extra index, this is made so that problems such as FMP, 
MSP can introduce variables with one extra indexing dimension. 

### Returns:
All possible indices for the decision variables represented by the
instance of the coefficient holder.
"""
function IndicesList(
    holder::MatrixConstruct.VariableCoefficientHolder,
    k::Union{Nothing, Int}=nothing
)
    NdimList = [1:III for III in size(holder)]
    if k === nothing
        return CartesianOutterProductList(NdimList...)
    end
    push!(NdimList, k:k)
return CartesianOutterProductList(NdimList...) end



"""
    PrepareVariablesForTheModel!(
        model::Model,
        variable::Symbol,
        ccga_itr::Union{Nothing, Int}=0
    )

Given a JuMP model, prepare the q, u variable for that model. This methods make the `u,q` variables for all the 
instances like FMP, FSP, MSP etc. 

# Arguments
- `model::Model`: An instance of the JuMP model. 
- `variable::Symbol`: An instance of the symbol for the new variable that we wish to create. It has to be one of 
{u, v}
- `ccga_itr::Union{Nothing, Int}=0`: The number of inner/outer interations for the CCGA. Depends on the context. 

# Returns
- Returns the variable 'u' or 'q' packed into Vector{JuMP.VariableRef}.

"""
function PrepareVariablesForTheModel!(
    model::Model,
    variable::Symbol,
    ccga_itr::Union{Nothing, Int}=0
)
    ccga_itr = ccga_itr == 0 ? "" : "[$(ccga_itr)]"
    if variable == :u
        u = Vector{JuMP.VariableRef}()
        for v in MatrixConstruct.u
            push!(u, @variable(model, [IndicesList(v)], base_name="$(v.v)$(ccga_itr)", lower_bound=0)...)
        end
        return u
    end
    
    if variable == :q
        q = Vector{JuMP.VariableRef}()
        for v in MatrixConstruct.q
            push!(q, @variable(model, [IndicesList(v)], Bin, base_name="$(v.v)$(ccga_itr)")...)
        end
        return q
    end
    error("The symbol $(variable) is passed into the function `PrepareVariablesForTheModel`, and we don't"*
    " support this symbol. ")
return nothing end



# ======================================================================================================================
# The Template For a PROBLEM type in general
# ======================================================================================================================
"""
An aobstract type problem is just a super type for all the problems in this project. 
Check `@ProblemTemplate` Macro or struct `ProblemTemplate`. 

"""
abstract type Problem
    
end

"""
The problem template is created to specify 2 of the common fields for all the optimization problem in this project. 
All the optimization problem in this project requires the field of: 
# Fields:
- `model::Model`: An instance of the JuMP model. 
- `con::Vector{ConstraintRef}`: Just an vector storing all the references to the contraints we are interested, which 
are also in the JuMP model itself. 
"""
@premix mutable struct ProblemTemplate
    "An instance of the JuMP model. "
    model::Model
    "A vector storing all the references to the contraints in the same order that they are being added by our problem instances.  "
    con::Vector{ConstraintRef}
end


# ======================================================================================================================
# FMP, the Upper bound locator.
#   * Obtains a upper bound by giving a demands that can break the feasibility for all discrete secondary
#     decision variables.
# Supports:
#   * Solving the system to obtain
# ======================================================================================================================


"""
An abstract FMP should model some type of FMP, the upper bound searcher for the feasibility problems. 

# Fields
`w::Vector{Float64}`, `gamma::Vector{Float64}`, `q::Vector{Vector{Int}}`,`d_hat::Vector{Float64}`,
`v::Vector{VariableRef}`, `u::Vector{Vector{VariableRef}}`, `eta::VariableRef`, `k::Int`. 
"""
abstract type AbsFMP <: Problem
    # FMP like problem. Different way of working with solving the FMP. Check @AbsFMPTemplate. 

end


@premix mutable struct AbsFMPTemplate
    """
    Primary generator discrete decision variable flattened into a constant vector that is consistent with the `Getw`
    function implemented in the master problem and the main problem. This is a given constant for the FMP. 
    """
    w::Vector{Float64}
    """
    This is the uncertainty bound passed from the master problem, and it's the upper bound and lower bound for the demands 
    decision variables in the instance of FMPH, it's dimension and ordering should be consistent with `GetGamma` function
    for the type MSP. 
    """
    gamma::Vector{Float64}
    """
    This is the secondary discrete decision variable for the generators. It's a constant. It's created to be 
    all zeros before any cuts are introduced to the instances of FMPs, after that, for each cut, there is an 
    constaint q introduced by the FSP. A given constant. 
    """
    q::Vector{Vector{Int}}
    "The center of the uncertainty interval for the vector. A given constant. "
    d_hat::Vector{Float64}
    "The continuous, feasibility decision variable. It's a vector where each element is a scalar variable for the feasibility
    for each of the cuts introduced to the instance of FMPs. "
    v::Vector{VariableRef}
    "The continuous decision variables for the secondary generators. Each cuts is associated with an separate set of secondary decision variables. "
    u::Vector{Vector{VariableRef}}
    # d::Vector{VariableRef}
    """
    The lower bound for all of the different feasibility decision variable `v` on each of the cuts. Maximizing this quantity is the objectives of FMPs.
    """
    eta::VariableRef
    # lambda::Vector{Vector{VariableRef}}
    "An counter for the number of cuts that is already introduced to FMPs, and the initial cut with `q` being all zeros is counted as the first one."
    k::Int                                        

end
