using JuMP
using Kronecker
using LinearAlgebra

model = Model()
γ = @variable(model, γ[1:4]>=0)
Γ  = kronecker(γ[:], ones(3))