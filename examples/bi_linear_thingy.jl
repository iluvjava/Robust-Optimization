using LinearAlgebra, JuMP, HiGHS

model = Model(HiGHS.Optimizer)
x = @variable(model, x[1:2]>=0)
y = @variable(model, y[1:2]>=0)
fix.(x, 0, force=true)

@constraint(model, dot(x, y) <= 10)
