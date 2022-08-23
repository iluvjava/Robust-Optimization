include("ccga_Innerloops.jl")

### =====================================================================================================================
### Problems we are identifying: 
### Trying to executing the CCGA Inner forloops and see what could be causing the infeasibility to the cut to the master 
### problem. 

ϵ = 0.1
γ̄ = 10
d̂ = 200*(size(MatrixConstruct.H, 2)|>ones)
d̂ = d̂ + randn(size(d̂))*10
model_mp = Model(HiGHS.Optimizer); mp = MP(model_mp, γ̄)
model_msp = Model(HiGHS.Optimizer); msp = MSP(model_msp, d̂, γ̄)
PortOutVariable!(mp, :d) do d fix.(d, d̂, force=true) end
PortOutVariable!(mp, :v) do v fix.(v, 0, force=true) end
Solve!(mp)
w̄ = Getw(mp)
Solve!(msp)

Results = CCGAInnerLoop(γ̄, w̄, d̂, sparse_vee=true)

fig = plot(Results.upper_bounds, label="upper_fmp", marker=:x)
plot!(fig, Results.lower_bounds, label="lower_fsp", marker=:x)
fig|>display

IntroduceCut!(
    msp, 
    Getu(Results.fsp), 
    Getq(Results.fsp), 
    GetRhoPlus(Results.fmp), 
    GetRhoMinus(Results.fmp)
)

DebugReport(msp, "msp_with_cut")
DebugReport(mp, "main_problem_for_reference")

