include("../src/ccga_modeling.jl")
import MathOptInterface as MOI
using Statistics
using ProgressMeter




"""
    Try to perturbed a given demand with a given model by constructing a hypercube with the demand being the center
    and sample randomly and count the percentage of feasible demands for the system. 
    Sample region: 
        -1ϵ + d̂ + 1ϵ
         
"""
function DemandPerturbations(
    this::MP, 
    demand::Vector{N}, 
    epsilon::M=10, 
    max_iter::Int=1000
) where {M<: Number, N<:Number}

    FeasibleCount = 0
    @showprogress for __ in 1: max_iter
        TestDemand = demand + (rand(length(demand)).- 0.5)*2*epsilon
        TestDemand = max.(TestDemand, 0)
        if DemandFeasible(this, TestDemand)
            FeasibleCount += 1
        end
    end
    "Demand perturbation, with ϵ=$(epsilon), feasibility percent: $(FeasibleCount/max_iter)" |> println
return FeasibleCount/max_iter end


"""
    Given a demand, it tries to figure out the largest hypercube region centered at the given demand, probabilistically. 
    * Returns the best demand interval for the given demand. 
"""
function DemandsBestIntervalSearch(
    mp::MP, 
    demand::Vector; 
    epsilon=100, 
    delta=100, 
    max_iter::Int=20
)
    Itr = 0
    while delta > 1
        Percent = DemandPerturbations(mp, demand, epsilon)
        if Percent == 1
            if Itr == max_iter
                break 
            end
            epsilon += delta
        else
            epsilon -= delta
            delta /= 2
        end
        Itr += 1
    end

return epsilon end


"""
    Compute candicate demand by computing the average of random objectives for the main problem. 
"""
function DemandRandomObjectiveAverage(mp::MP, N=1000)
    d = mp.d
    model = mp |> GetModel
    v = Vector{Vector}()
    @showprogress for __ = 1:N
        c = randn(length(d))
        @objective(model, Max, dot(c, d))
        Solve!(mp)
        push!(v, mp.d.|>value)
    end

return mean(hcat(v...), dims = 2)  end


function RunThis()
    model = Model(
        optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false),
    )
    
    
    mp = MP(model)
    Solve!(mp)
    @info "Computing average demand from random objective. "
    
    d̂ = 30*ones(length(mp.d)) # 
    #d̂ = DemandRandomObjectiveAverage(mp)
    @info "Suggusted Average Demands: "
    
    d̂ |> display
    @info "Performing Demand best interval search on above suggusted demand vector."
    
    BestDemandInterval = DemandsBestIntervalSearch(mp, d̂[:]; epsilon=30, delta=30)
    @info "The best demand interval seems to be $(BestDemandInterval)"
return end

RunThis()

