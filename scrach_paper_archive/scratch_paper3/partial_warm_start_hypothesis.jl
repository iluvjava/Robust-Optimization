using JuMP, Gurobi


function MakeBinPacking(N::Int=100, weight_limit::Real=5)
    l = ones(N)
    l[1:3] .= 0
#    l[-1] = -1
    model = Model(Gurobi.Optimizer)
    x = model[:x] = @variable(model, [1:length(l)], Bin, base_name="x")
    model[:con] = @constraint(model, dot(l, x) <= weight_limit)
    @objective(model, Max, dot(l, x))

    return model
end


model = MakeBinPacking()
set_start_value.(model[:x][1:3], 0)
optimize!(model)
