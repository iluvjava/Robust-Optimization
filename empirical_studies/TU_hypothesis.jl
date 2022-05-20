using JuMP, HiGHS

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
    λ = @variable(model, λ[1:J])
    @constraint(model,[b=1:B, j=1:J],ρ⁺[b]<=ξ⁺[b, j])
    @constraint(model,[b=1:B, j=1:J],ξ⁺[b, j]<=ρ⁺[b])
    @constraint(model,[b=1:B, j=1:J],λ[j] - (1 - ρ⁺[b])<=ξ⁺[b, j])
    @constraint(model,[b=1:B, j=1:J],ξ⁺[b, j]<=λ[j] + (1 - ρ⁺[b]))

    @constraint(model,[b=1:B, j=1:J],ρ⁻[b]<=ξ⁻[b, j])
    @constraint(model,[b=1:B, j=1:J],ξ⁻[b, j]<=ρ⁻[b])
    @constraint(model,[b=1:B, j=1:J],λ[j] - (1 - ρ⁻[b])<=ξ⁻[b, j])
    @constraint(model,[b=1:B, j=1:J],ξ⁻[b, j]<=λ[j] + (1 - ρ⁻[b]))


return model,vcat(ρ⁺[:], ρ⁻[:], ξ⁺[:], ξ⁻[:]) end

model, variables = MakeModel()

for II =1:100
    c = randn(variables|>length)
    @objective(model, Min, dot(c, variables))
    optimize!(model)
    Obj_LP = objective_value(model)

    model, variables = MakeModel(ip=true)
    @objective(model, Min, dot(c, variables))
    optimize!(model)
    Obj_IP = objective_value(model)
    @assert Obj_IP ≈ Obj_LP
end