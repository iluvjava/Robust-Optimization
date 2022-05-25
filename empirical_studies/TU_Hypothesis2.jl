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



model = nothing
model_ip = nothing
for II =1:100
    model, variables = MakeModel()
    c = randn(variables|>length)

    
    @objective(model, Min, dot(c, variables))
    optimize!(model)
    Obj_LP = objective_value(model)

    model_ip, variables = MakeModel(ip=true)
    @objective(model_ip, Min, dot(c, variables))
    optimize!(model_ip)
    Obj_IP = objective_value(model_ip)
    if abs(Obj_IP - Obj_LP) > 1e-10
        display("failed at $(II)")
        display("ip objective: $(objective_value(model_ip))")
        display("lp objective: $(objective_value(model))")
        break
    end
    

end

# Seems like this MIP problem is having total unimodularity for the variable: ρ. 

value.(model_ip[:ρ⁺])|>display
value.(model[:ρ⁻])|>display

