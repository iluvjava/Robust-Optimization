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
Given a JuMP model, prepare the q, u variable for that model. This methods make the `u,q` variables for all the 
instances like FMP, FSP, MSP etc. 

### Arguments
- `model::Model`: An instance of the JuMP model. 
- `variable::Symbol`: An instance of the symbol for the new variable that we wish to create. It has to be one of 
{u, v}
- `ccga_itr::Union{Nothing, Int}=0`: The number of inner/outer interations for the CCGA. Depends on the context. 

### Returns
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
"""
abstract type Problem
    # has a JuMP model in it.
    # MUST IMPLEMENT: GetModel(::Problem)
    # has a con vector that stores the references to all the constraints created for the model. 
end

"""
The problem template is created to specify 2 of the common fields for all the optimization problem in this project. 
All the optimization problem in this project requires the field of: 
### Fields:
- `model::Model`: An instance of the JuMP model. 
- `con::Vector{ConstraintRef}`: Just an vector storing all the references to the contraints we are interested, which 
are also in the JuMP model itself. 
"""
@premix mutable struct ProblemTemplate
    model::Model
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
### Fields
- `w::Vector{Float64}`:Primary Generator decision variables, (GIVEN CONSTANT)
- `gamma::Vector{Float64}` The bound for the demands on all the buses during a specific time, (GIVEN CONSTANT)
- `q::Vector{Vector{Int}}`: The secondary discrete decision variables, (GIVEN CONSTANT)
- `d_hat::Vector{Float64}`: The average testing demands vector, (GIVEN CONSTANT)
- `v::Vector{VariableRef}`: The slack decision variables for each of the previous demands.
- `u::Vector{Vector{VariableRef}}`: The secondary continuous decision variables.
- `d::Vector{VariableRef}`: The demand decision variable, as a giant vector.
- `eta::VariableRef`: The eta lower bound for all feasibility.
- `lambda::Vector{Vector{VariableRef}}`: The dual decision variables.
- `k::Int`:The iteration number from the ccga, this instance keep tracks of it by itself, because it uses this parameter
internally frequently. 
"""
abstract type AbsFMP <: Problem
    # FMP like problem. Different way of working with solving the FMP. 

end


@premix mutable struct AbsFMPTemplate
    w::Vector{Float64}                            # Primary Generator decision variables.               (GIVEN CONSTANT)
    gamma::Vector{Float64}                        # the bound for the demands on all the buses during a specific time  (GIVEN CONSTANT)
    q::Vector{Vector{Int}}                        # the secondary discrete decision variables           (GIVEN CONSTANT)
    d_hat::Vector{Float64}                        # the average testing demands vector.                 (GIVEN CONSTANT)

    v::Vector{VariableRef}                         # The slack decision variables for each of the previous demands.
    u::Vector{Vector{VariableRef}}                # The secondary continuous decision variables.
    # d::Vector{VariableRef}                      # The demand decision variable, as a giant vector.
    eta::VariableRef                              # The eta lower bound for all feasibility.
    # lambda::Vector{Vector{VariableRef}}           # the dual decision variables.
    
    k::Int                                        # The iteration number from the ccga.

end
