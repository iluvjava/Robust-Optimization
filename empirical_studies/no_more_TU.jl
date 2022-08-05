using JuMP, HiGHS, LinearAlgebra
include("../src/matrix_construction_export.jl")
H = MatrixConstruct.H

"""
    Make the original model with full ξ. 
"""
function MakeDenseModel(J::Int, B::Int)
    model = Model(HiGHS.Optimizer)
    ρ⁺ = @variable(model, ρ⁺[1:B], lower_bound=0, upper_bound=1)
    ρ⁻ = @variable(model, ρ⁻[1:B], lower_bound=0, upper_bound=1)
    
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

return model end

"""
Make the model using the sparse matrix's entries to construct the constraints. 
"""
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
    λ  = @variable(model, λ[1:J], lower_bound=-1, upper_bound=0)

    for (j, b) in IdxTuples
        @constraint(model, -ρ⁺[b] <= ξ⁺[(j, b)])
        @constraint(model, ξ⁺[(j, b)] <= ρ⁺[b])
        @constraint(model, λ[j] - (1 - ρ⁺[b])<=ξ⁺[(j, b)])
        @constraint(model, ξ⁺[(j, b)]<=λ[j] + (1 - ρ⁺[b]))

        @constraint(model, -ρ⁻[b] <= ξ⁻[(j, b)])
        @constraint(model, ξ⁻[(j, b)] <= ρ⁻[b])
        @constraint(model, λ[j] - (1 - ρ⁻[b]) <= ξ⁻[(j, b)])
        @constraint(model, ξ⁻[(j, b)] <= λ[j] + (1 - ρ⁻[b]))
    end
    @constraint(model, [b=1:B], ρ⁺[b] + ρ⁻[b] == 1)

return model end


"""
    Provides an objective for a given model and then it runs it. 
"""
function RunModel(model, Ĥ, d, γ)
    @info "LP Model Solve: "
    λ = model[:λ]
    ξ⁺ = model[:ξ⁺]
    ξ⁻ = model[:ξ⁻]
    @objective(model, Max, dot(λ, H*d) + γ*dot(Ĥ, ξ⁺[:]) - γ*dot(Ĥ, ξ⁻[:]))
    # @objective(model, Max, dot(λ, H*d))
    optimize!(model)
return objective_value(model) end

### Setting up the objective 
d = 100 .+ 20*rand(size(H, 2))
γ = 40

### Solve the objective value of the LP programming problem. 
ModelSparseLP = MakeModel(H)
Objective1 = RunModel(ModelSparseLP, H[H.!=0], d, γ)

## Solve the MIP model and find the objective. 
ModelSparseMIP = MakeModel(H, ip=true)
Objective2 = RunModel(ModelSparseMIP, H[H.!=0], d, γ)

ModelDenseLP = MakeDenseModel(size(H, 1), size(H, 2))
Objective3 = RunModel(ModelDenseLP, H'[:], d, γ)

λ1 = ModelSparseLP[:λ].|>value
λ2 = ModelDenseLP[:λ].|>value
sum(λ1 .== -1)|>println
sum(λ2 .== -1)|>println



