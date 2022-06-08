include("ccga_modeling.jl")
### ====================================================================================================================
### The main problem does the following: 
###     * Test whether the configurations is feasible at all. 
###     * Test whether 

"""
    First Phase of the Main Problem: 
        * Does there exists an nontrivial demands for the system? 
"""
function MainProblemPhaseI()
    
return end


mp = MP()
Solve!(mp)
DebugReport(mp)
@info "Results"


if !isnan(objective_value(mp))
    @info "Primary discrete decision variables:"
    mp.w.|>value.|>println
    @info "Secondary Discreate decision variables:"
    mp.q.|>value.|>println;
    @info "Non-zero slack"
    findall(>(0), mp.v.|>value).|>println;
    @info "Coitinuous Decision Vars"
    zip(mp.u, mp.u.|>value).|>println;
else
    println("Problem is Infeasible. ")
end

