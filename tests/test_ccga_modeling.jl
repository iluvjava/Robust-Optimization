include("../src/ccga_modeling.jl")



### Test the existence of any feasible configurations for the system. 
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



# """
#     First iteration of the CCGA establish the following: 
#     * Using a random q for FMP, which is definitely infeasible. 
#     * Get a random setup for Q and perform FMP to get a branch cut for MP. 
#     * Cut it on MP and get a legit setup for the primary decision variable. 
# """
# # function FirstCCGA()
#     d̂ = 200.0*ones(size(RobustOptim.H, 2))
#     mp = MP(200)
#     @info "First MP solve"
#     Solve!(mp)
#     # CCGA First Iterations: Initializations
#     w, γ = Getw(mp).|>value, GetGamma(mp).|>value
#     @info "First FMP solve with given d̂, w from MP, γ from MP. "
#     fmp = FMP(w, γ, d̂)
#     Solve!(fmp)
#     d = GetDemandVertex(fmp)
#     @info "First FSP solve with w, γ from MP, d from FMP. "
#     fsp = FSP(w, γ, d.|>value)
#     Solve!(fsp)
#     q = Getq(fsp)
#     u = Getu(fsp)
#     Upper, Lower = objective_value(fmp), objective_value(fsp) # This is huge, without a doubt. 
#     @info "CCGAA: Upper bound: $(Upper), Lower Bound: $(Lower)"
#     ρ⁺ = GetRhoPlus(fmp)
#     ρ⁻ = GetRhoMinus(fmp)
#     @info "Introduce Cut to MP: "
#     IntroduceCut!(mp, u.|>value, q.|>value, d̂, ρ⁺, ρ⁻)
#     Solve!(mp)
# # return end 
# # CCGA Repeated iterations. 



