include("../src/ccga_modeling.jl")
import MathOptInterface as MOI
# using ProgressMeter

model = Model(
    optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false),
)
mp = MP(model)
Solve!(mp)
d = 75*ones(size(mp.d))



"""
    Try to perturbed a given demand with a given model, this is for testing the whether the hypercude cetnered 
    given demand is a good hypercube, it's done using some random sampling tricks. 
    
"""
function DemandPerturbations(
    this::MP, 
    demand::Vector{N}, 
    epsilon::M=10, 
    max_iter::Int=1000
    ) where {N<:Number, M<:Number}

    FeasibleCount = 0
    for __ in 1: max_iter
        TestDemand = demand + (rand(length(d)).- 0.5)*epsilon
        if DemandFeasible(this, TestDemand)
            FeasibleCount += 1
        end
    end
    "Demand perturbation, with ϵ=$(epsilon), feasibility percent: $(FeasibleCount/max_iter)" |> println
return FeasibleCount/max_iter end

function DemandRandomObjectiveAverage()

return end



FeasibilityPercent = Vector{Float64}()
for ϵ in LinRange(10, 50, 10)
    push!(FeasibilityPercent, DemandPerturbations(mp, d, ϵ))
end