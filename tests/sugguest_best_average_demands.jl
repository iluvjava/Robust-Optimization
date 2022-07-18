include("../src/ccga_modeling.jl")
import MathOptInterface as MOI
using Statistics
using ProgressMeter
using Plots


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
        TestDemand = demand + rand([-1, 1], length(demand))*epsilon
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

    function Reporter(mp::MP, epsilon::T1, average_demand::Vector{T2}) where {T1, T2 <: Number}
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
function Generate!(this::Reporter, sample_count=30)
    this.feasible_demands = Vector{Vector{Float64}}()
    this.infeasible_demands = Vector{Vector{Float64}}()
    feasible = this.feasible_demands
    infeasible = this.infeasible_demands
    mp = this.mp
    d̂ = this.d_hat
    ϵ = this.epsilon
    function DemandClassify(d)
        if DemandFeasible(mp, d) 
            push!(feasible, d)
        else
            push!(infeasible, d)
        end
    return end

    # Phase I: classify extreeme demands. 
    DemandClassify(d̂ .+ ϵ)
    DemandClassify(max.(d̂ .- ϵ, 0))

    
    # Phase II: Classify vertices. 
    for __ in 1:sample_count
        DemandClassify(
            max.(
                d̂ + rand([-1, 1], length(d̂))*ϵ, 0
            )
        )
    end

    # Phase II: Classify demands in the interior of the hypercube. 
    l = @. max(d̂ - ϵ, 0) 
    u = d̂ .+ ϵ
    for __ in 1:sample_count
        DemandClassify(l + rand(length(d̂)).*(u - l))
    end
    
return end


"""
    Visualize the feasible and infeasible demands and print it to the 
    through IO in the current project. 
"""
function VisualizeAll(this::Reporter)
    if !isdefined(this, :feasible_demands)
        this |> Generate!
    end
    feasible_demands = hcat(this.feasible_demands...)
    infeasible_demands = hcat(this.infeasible_demands...)
    fig = plot(
        feasible_demands,
        linewidth=0, 
        markershape=:hline, 
        # dpi=200, 
        # size=(800, 600),
        legend=nothing, 
        color=:blue, 
        title="Feasible and Infeasible Demands"
    )
    plot!(
        fig, 
        infeasible_demands,
        linewidth=0.2, 
        markershape=:x, 
        color=:red
    )
    plot!(
        fig, 
        this.d_hat, 
        linewidth=0, 
        markershape=:circle, 
        color=:black
    )

return fig end


"""
    Print out the summary statistics for the current list of 
    feasible demands and infeasible demands. 
"""
function DescriptiveStatisticsForDemandsSum(this::Reporter)
    FeasibleDemandsSum = sum(hcat(this.feasible_demands...), dims=1)[:]
    InFeasibleDemandsSum = sum(hcat(this.infeasible_demands...), dims=1)[:]
    s = ""
    s *= "Feasible Demands sum(d) for all:\n" 
    s *= "   Avg: $(mean(FeasibleDemandsSum))\n" 
    s *= "   std: $(std(FeasibleDemandsSum))\n" 
    s *= "   75% Quantile: $(quantile(FeasibleDemandsSum, 0.75))\n" 
    s *= "   25% Quantile: $(quantile(FeasibleDemandsSum, 0.25))\n" 
    s *= "   Min: $(minimum(FeasibleDemandsSum))\n" 
    s *= "   Max: $(maximum(FeasibleDemandsSum))\n" 
    s *= "   Instance count: $(FeasibleDemandsSum |> length)\n" 
    s *= "\nInfeasible Demands sum(d) for all:\n" 
    s *= "   Avg: $(mean(InFeasibleDemandsSum))\n" 
    s *= "   std: $(std(InFeasibleDemandsSum))\n" 
    s *= "   75% Quantile: $(quantile(InFeasibleDemandsSum, 0.75))\n" 
    s *= "   25% Quantile: $(quantile(InFeasibleDemandsSum, 0.25))\n" 
    s *= "   Min: $(minimum(InFeasibleDemandsSum))\n" 
    s *= "   Max: $(maximum(InFeasibleDemandsSum))\n" 
    s *= "   Instance count: $(InFeasibleDemandsSum |> length)\n" 

return s end

"""
    Save the reported results for all the demands for the system. 
"""
function SaveResults(this::Reporter, filename::String="")
    p = mkpath("average_demands_reports")
    fig = this |> VisualizeAll
    stats = this |> DescriptiveStatisticsForDemandsSum
    savefig(fig, joinpath(p, "$(filename)(Plots).png"))
    open("$(filename)(stats).txt", "w") do io
        print(io, stats)
    end

return this end


# Reporter Structs END =================================================================================================

function DemandProfile1(testname::String="1")

    model = Model(
        optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false),
    )
    mp = MP(model)

    @info "Computing average demand from random objective. "
    d̂ = 40*ones(length(mp.d))
    @info "Suggusted Average Demands: "
    d̂ |> display
    @info "Performing Demand best interval search on above suggusted demand vector."
    BestDemandInterval = DemandsBestIntervalSearch(mp, d̂[:]; epsilon=10, delta=30)
    @info "The best demand interval seems to be $(BestDemandInterval)"

    reporter = Reporter(mp, BestDemandInterval + 40, d̂[:])
    SaveResults(reporter, testname)

return reporter end


function DemandProfile2()
    model = Model(
        optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false),
    )
    mp = MP(model)
    Solve!(mp)
    @info "Computing average demand from random objective. "
    d̂ = DemandRandomObjectiveAverage(mp)./2
    @info "Suggusted Average Demands: "
    d̂ |> display
    @info "Performing Demand best interval search on above suggusted demand vector."
    BestDemandInterval = DemandsBestIntervalSearch(mp, d̂[:]; epsilon=10, delta=30)
    @info "The best demand interval seems to be $(BestDemandInterval)"

    reporter = Reporter(mp, BestDemandInterval + 20, d̂[:])
    reporter |> Generate!
    reporter |> VisualizeAll
    reporter |> DescriptiveStatisticsForDemandsSum |> print

return end


function DemandProfile3()
    model = Model(
        optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false),
    )
    mp = MP(model)
    Solve!(mp)
    @info "Computing average demand from random objective. "
    d̂ = 40*ones(length(mp.d)) + rand(length(mp.d))*20
    @info "Suggusted Average Demands: "
    d̂ |> display
    @info "Performing Demand best interval search on above suggusted demand vector."
    BestDemandInterval = DemandsBestIntervalSearch(mp, d̂[:]; epsilon=10, delta=30)
    @info "The best demand interval seems to be $(BestDemandInterval)"

    reporter = Reporter(mp, BestDemandInterval + 5, d̂[:])
    reporter |> Generate!
    reporter |> VisualizeAll
    reporter |> DescriptiveStatisticsForDemandsSum |> print

return end


reporter = DemandProfile1()




