using JuMP, GLPK 

function AddVariables(model::JuMP.Model, var::Symbol)
    @variable(model, var[1:3, 1:3])
end

model = Model()
AddVariables(model, :x)
