using GLPK, JuMP
model = Model(GLPK.Optimizer)
x = @variable(model, x[1:5, 1:5], Bin)
y = @variable(model, y[1:5], Bin)
hcat(x..., y...)