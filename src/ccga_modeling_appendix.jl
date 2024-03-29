### ====================================================================================================================
### GENETRIC FUNCTIONS FOR ALL ABOVE TYPES AND SUBSETS OF THOSE TYPES. 
###     - Generics Functions for abstract type of problems that contain a JuMP inner model, specific overrides too.
### ====================================================================================================================


"""
Solve the inner JuMP optimization problem by directly invoking `optimize!` in JuMP on the model. This 
function is generic to the type `Problem`. It returns `optimize!(this|>GetModel)`

"""
function Solve!(this::Problem)
return optimize!(this|>GetModel) end


"""
Get the JuMP instance from this model.
"""
function GetModel(this::Problem)
return this.model end

"""
Get a list of references to the constraints of the model. These constraints should be 
adhere to the same order when they are being added when the model is constructed. 
"""
function GetConstraints(this::Problem)
    return this.con
end

"""
return the list of the constraints made for an instance of the FMPH1
"""
function GetConstraints(this::FMPH1)
    return hcat(this.con, this.dual_cons)
end


"""
Produce a report for the MP (master problem), which also checks feasibility of the original problem and stuff.
---
* Print out the string representation of the model.
* Print out the constraints for the model.
* Print out the objective for the model.
"""
function DebugReport(this::Problem, filename="debug_report")
    open("$(filename).txt", "w") do io
        write(io, this|>GetModel|>repr)
        write(io, "\n")
        this|>GetConstraints.|>(to_print) -> println(io, to_print)
        if !(this|>objective_value|>isnan)  # it's solved. 
            ["$(t[1]) = $(t[2])" for t in zip(this|>all_variables, this|>all_variables.|>value)].|>(to_print) -> println(io, to_print)
        end
    end
return this end


"""
Print out the constraints, the `.con`, and stream them to a file. 
"""
function PrintConstraintsGroup(con::Vector{JuMP.ConstraintRef}, filename::String="Constraints_report")
    open("$(filename).txt", "w") do io
        con.|>repr.|>(to_print) -> println(io, to_print)
    end
return end


"""
    Port out an variable for the user to write a functional to access a specific variable. 
        * Port out the decision variable from the model, if you give it the
            * Problem 
            * The symbol for the field for the problem. 
    User's Functional Should: 
        * Accept on parameters of type `JuMP.VariableRef`, and perform the wanted modificatons 
        for the decision variable. 
"""
function PortOutVariable!(port_out::Function, this::Problem, variable::Symbol)
    if variable in this|>typeof|>fieldnames
        return port_out(getfield(this, variable))
    end
    return (this |> GetModel)[variable]
return end


"""
    * Port out the model of the system as the first parameters of the context function. 
    * Port out the field of the instance as a dictionary of symbols ->  Field references as the second 
    parameters for the function. 
"""
function PortOutEverything!(port_out::Function, this::Problem)
    FieldDict = Dict([field => getfield(this, field) for field in fieldnames(this)])
return port_out(GetModel(this), FieldDict) end


# ----------------------------------------------------------------------------------------------------------------------

"""
    Get the value of the decision variable q, it's suitable for the type: Union{FSP, MP, FMP}. 
"""
function Getq(this::Union{FSP, MP})
return this.q.|>value.|>(x) -> round(x, digits=1)  end

"""
    Get the value of the u decision variables. 
"""
function Getu(this::Union{FSP, MP})
return this.u.|>value end


function Getv(this::Union{FSP, MP})
return this.v.|>value end


"""
    If the instance is an FMP, then it returns the variable v that is the most recent. 
"""
function Getv(this::FMP)
return this.v[end].|>value end



### ====================================================================================================================
### Methods Shared with JuMP.jl Model type
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
    List out the variable references for all the decision variables that are in the problem's model

"""
function all_variables(this::Problem) return this|>GetModel|>JuMP.all_variables end


"""
    Inherit the same indexer from the JuMP.Model. 
"""
function Base.getindex(this::Problem, param::Any)
return GetModel(this)[param] end

### ====================================================================================================================
### Methods that doesn't fit into any of the above category
### ====================================================================================================================

function PrintoutVariables(v::Vector{VariableRef})
    for (v, val) in zip(v, v.|>value)
        println("$v=$val")
    end
end