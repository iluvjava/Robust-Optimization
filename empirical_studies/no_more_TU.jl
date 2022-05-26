using JuMP, HiGHS, LinearAlgebra
include("../src/matrix_construction_export.jl")
H = RobustOptim.H


function MakeModel(H; ip::Bool=false)
    J, B = size(H)
    model = Model(HiGHS.Optimizer)
    if ip
        ρ⁺ = @variable(model, ρ⁺[1:B], Bin)
        ρ⁻ = @variable(model, ρ⁻[1:B], Bin)
    else
        ρ⁺ = @variable(model, ρ⁺[1:B], lower_bound=0, upper_bound=1)
        ρ⁻ = @variable(model, ρ⁻[1:B], lower_bound=0, upper_bound=1)
    end
    IdxTuples = Vector{Tuple}()
    for CartIdx in findall(!=(0), H)
        j, b = Tuple(CartIdx)
        push!(IdxTuples, (j, b))
    end
    ξ⁺ = @variable(model, ξ⁺[IdxTuples])
    ξ⁻ = @variable(model, ξ⁻[IdxTuples])
    λ = @variable(model, λ[1:J], lower_bound=-1, upper_bound=0)

    for (j, b) in IdxTuples
        @constraint(model, -ρ⁺[b] <= ξ⁺[(j, b)])
        @constraint(model, ξ⁺[(j, b)] <= ρ⁺[b])
        @constraint(model, λ[j] - (1 - ρ⁺[b])<=ξ⁺[(j, b)])
        @constraint(model, ξ⁺[(j, b)]<=λ[j] + (1 - ρ⁺[b]))

        @constraint(model, -ρ⁻[b] <= ξ⁻[(j, b)])
        @constraint(model, ξ⁻[(j, b)] <= ρ⁻[b])
        @constraint(model, λ[j] - (1 - ρ⁻[b]) <= ξ⁻[(j, b)])
        @constraint(model, ξ⁻[(j, b)] <= λ[j] + (1 - ρ⁻[b]))
        @constraint(model,[b=1:B], ρ⁺[b] + ρ⁻[b] == 1)
    end

return model, H[H.!=0] end


"""
    Provide specific objective for a model, returns that model
"""
function RunThis(ip::Bool)
    H = RobustOptim.H
    d = 100 .+ 20*rand(B)
    γ = 20
    I = ones(J)
    @info "LP Model Solve: "
    model, Ĥ = MakeModel(H, ip=ip)
    λ = model[:λ]
    ξ⁺ = model[:ξ⁺]
    ξ⁻ = model[:ξ⁻]
    @objective(model, Max, dot(λ, H*d) + γ*dot(Ĥ, ξ⁺) - γ*dot(Ĥ, ξ⁻[:]))
    optimize!(model)
    ObjectVal = objective_value(model)
    println("Obj Val LP: $(ObjectVal)")
    
return model end

model = RunThis(false)
optimize!(model)
model_mip = RunThis(true)
optimize!(model_mip)