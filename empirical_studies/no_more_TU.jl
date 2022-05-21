using JuMP, HiGHS, LinearAlgebra
include("../src/matrix_construction_export.jl")
H = RobustOptim.H


function MakeModel(B::Int=6, J::Int=100; ip::Bool=false)
    model = Model(HiGHS.Optimizer)
    if ip
        ρ⁺ = @variable(model, ρ⁺[1:B], Bin)
        ρ⁻ = @variable(model, ρ⁻[1:B], Bin)
    else
        ρ⁺ = @variable(model, ρ⁺[1:B], lower_bound=0, upper_bound=1)
        ρ⁻ = @variable(model, ρ⁻[1:B], lower_bound=0, upper_bound=1)
    end
    ξ⁺ = @variable(model, ξ⁺[1:B, 1:J])
    ξ⁻ = @variable(model, ξ⁻[1:B, 1:J])
    λ = @variable(model, λ[1:J], lower_bound=-1, upper_bound=0)

    @constraint(model,[b=1:B, j=1:J], -ρ⁺[b] <= ξ⁺[b, j])
    @constraint(model,[b=1:B, j=1:J], ξ⁺[b, j] <= ρ⁺[b])
    @constraint(model,[b=1:B, j=1:J], λ[j] - (1 - ρ⁺[b])<=ξ⁺[b, j])
    @constraint(model,[b=1:B, j=1:J], ξ⁺[b, j]<=λ[j] + (1 - ρ⁺[b]))

    @constraint(model,[b=1:B, j=1:J], -ρ⁻[b] <= ξ⁻[b, j])
    @constraint(model,[b=1:B, j=1:J], ξ⁻[b, j] <= ρ⁻[b])
    @constraint(model,[b=1:B, j=1:J], λ[j] - (1 - ρ⁻[b]) <= ξ⁻[b, j])
    @constraint(model,[b=1:B, j=1:J], ξ⁻[b, j] <= λ[j] + (1 - ρ⁻[b]))
    @constraint(model,[b=1:B], ρ⁺[b] + ρ⁻[b] == 1)


return model, vcat(ξ⁺[:], ξ⁻[:], λ[:]) end


"""
    Provide specific objective for a model, returns that model
"""
function RunThis()

return end