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
    ξ⁺ = @variable(model, ξ⁺[1:B])
    ξ⁻ = @variable(model, ξ⁻[1:B])
    λ = @variable(model, λ[1:J], lower_bound=-1, upper_bound=0)

    @constraint(model,[b=1:B, j=1:J], -ρ⁺[b] <= ξ⁺[b])
    @constraint(model,[b=1:B, j=1:J], ξ⁺[b] <= ρ⁺[b])
    @constraint(model,[b=1:B, j=1:J], λ[j] - (1 - ρ⁺[b])<=ξ⁺[b])
    @constraint(model,[b=1:B, j=1:J], ξ⁺[b]<=λ[j] + (1 - ρ⁺[b]))

    @constraint(model,[b=1:B, j=1:J], -ρ⁻[b] <= ξ⁻[b])
    @constraint(model,[b=1:B, j=1:J], ξ⁻[b] <= ρ⁻[b])
    @constraint(model,[b=1:B, j=1:J], λ[j] - (1 - ρ⁻[b]) <= ξ⁻[b])
    @constraint(model,[b=1:B, j=1:J], ξ⁻[b] <= λ[j] + (1 - ρ⁻[b]))
    @constraint(model,[b=1:B], ρ⁺[b] + ρ⁻[b] == 1)


return model, vcat(ξ⁺[:], ξ⁻[:], λ[:]) end


"""
    Provide specific objective for a model, returns that model
"""
function RunThis()
    J = size(H, 1)
    B = size(H, 2)
    d = 100 .+ 20*rand(B)
    γ = 20
    I = ones(J)
    @info "LP Model Solve: "
    model, vars = MakeModel(B, J)
    λ = model[:λ]
    ξ⁺ = model[:ξ⁺]
    ξ⁻ = model[:ξ⁻]
    @objective(model, Max, dot(λ, H*d) + γ*sum(H*ξ⁺) - γ*sum(H*ξ⁻))
    optimize!(model)
    ObjectVal = objective_value(model)
    
    @info "MIP Model Solve: "
    model_mip, vars = MakeModel(B, J, ip=true)
    λ = model_mip[:λ]
    ξ⁺ = model_mip[:ξ⁺]
    ξ⁻ = model_mip[:ξ⁻]
    @objective(model_mip, Max, dot(λ, H*d) + γ*sum(H*ξ⁺) - γ*sum(H*ξ⁻))
    optimize!(model_mip)
    ObjectVal_mip = objective_value(model_mip)

    println("Obj Val LP: $(ObjectVal)")
    println("Obj Val MIP: $(ObjectVal_mip)")
return model, model_mip end

model, model_mip = RunThis()