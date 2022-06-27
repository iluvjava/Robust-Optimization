include("../src/ccga_modeling.jl")
import MathOptInterface as MOI
using Statistics
using ProgressMeter


"""
    Try to perturbed a given demand with a given model by constructing a hypercube with the demand being the center
    and sample randomly and count the percentage of feasible demands for the system. 
    Sample region: 
        max(-1ϵ, 0) + d̂ + 1ϵ         
"""
function DemandPerturbations(
    this::MP, 
    demand::Vector{N}, 
    epsilon::M=10, 
    max_iter::Int=5000
) where {M<: Number, N<:Number}

    FeasibleCount = 0

    TestDemand = max.(demand .+ epsilon, 0)
    if !DemandFeasible(this, TestDemand)
        "Demand for perturbation with ϵ=$(epsilon) failed for extreme upper demand." |> println
        return false end
        TestDemand = max.(demand .- epsilon, 0)
    if !(DemandFeasible(this, TestDemand))
        "Demand for perturbation with ϵ=$(epsilon) failed for extreme lower demand." |> println
        return false end
    
    "Sampling from vertices of the hypercube: max(-1ϵ, 0) + d̂ + 1ϵ" |> println
    @showprogress for __ in 1: max_iter
        TestDemand = demand + rand([0, 1], length(demand))*epsilon
        TestDemand = max.(TestDemand, 0)
        if DemandFeasible(this, TestDemand)
            FeasibleCount += 1
        else
            "Demand perturbation with ϵ=$(epsilon) failed before reaching $(max_iter) iterations. " |> println
            return false
        end
    end

    "Demand perturbation with ϵ=$(epsilon) successfully passed $(max_iter) iterations. " |> println
return true end


"""
    Given a demand, it tries to figure out the largest hypercube region centered at the given demand, probabilistically. 
    * Returns the best demand interval for the given demand. 
"""
function DemandsBestIntervalSearch(
    mp::MP, 
    demand::Vector; 
    epsilon=0, 
    delta=100, 
    max_iter::Int=20,
    tol=2^(-5)
)
    Itr = 1
    HistoricRecords = Dict{Float64, Bool}()
    while delta > tol && Itr <= max_iter
        if epsilon in keys(HistoricRecords)
            DemandPassed = HistoricRecords[round(epsilon, digits=10)]
        else
            DemandPassed = DemandPerturbations(mp, demand, epsilon)
            HistoricRecords[round(epsilon, digits=10)] = DemandPassed
        end
        if DemandPassed
            epsilon += delta
        else
            epsilon -= delta
            delta /= 2
        end
        Itr += 1
    end

return maximum([d for d in keys(HistoricRecords) if HistoricRecords[d]]) end


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

return mean(hcat(v...), dims=2)  end


### Results Reporter ===================================================================================================
### Reports the stability of a demand by sampling then and visualizing them, all of them will be tested on feasibility.
### * construct a plot for all the feasible demand sampled from the inside of the hypercube. 
### * construct a plot for all the infeasible demand sampled from the inside of the hypercube. 
### The plot: 
###     columns of dots denoting each elements of the vectors involved. 
### ====================================================================================================================


mutable struct Reporter
    mp::MP
    d_hat::Vector
    epsilon::N where N<:Number
    
    feasible_demands::Vector{Vector{N}} where N <: Number
    infeasible_demands::Vector{Vector{N}} where N <: Number

    function Reporter(mp::MP, epsilon::Number, average_demand::Vector{Number})
        this = new()
        this.mp = mp
        this.d_hat = average_demand
        this.epsilon = epsilon

    return this end
    

end

"""
    Generate Statistics report about this given demands with a given interval. 
        * Test extreme demand
        * sample from the vertices of the hypercube. 
        * sample from the interior of the hypercube. 
"""
function Generate(this::Reporter)
    this.feasible_demands = Vector{Vector{Float64}}()
    this.infeasible_demands = Vector{Vector{Float64}}()
    feasible = this.feasible_demands
    infeasible = this.infeasible_demands
    # Phase I, testing extreme demand value. 
    mp = this.mp
    

    
return end



# Reporter Structs END =================================================================================================

function RunThis()
    model = Model(
        optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false),
    )
    
    mp = MP(model)
    Solve!(mp)
    @info "Computing average demand from random objective. "
    
    d̂ = 40*ones(length(mp.d))  
    @info "Suggusted Average Demands: "
    
    d̂ |> display
    @info "Performing Demand best interval search on above suggusted demand vector."
    
    BestDemandInterval = DemandsBestIntervalSearch(mp, d̂[:]; epsilon=10, delta=30)
    @info "The best demand interval seems to be $(BestDemandInterval)"


return end

RunThis()

