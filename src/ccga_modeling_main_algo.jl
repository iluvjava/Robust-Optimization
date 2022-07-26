### READ CCGA.md For more info. 

include("utilities.jl")
include("matrix_construction_export.jl")
include("ccga_modeling.jl")

ϵ = 0.1
M = 25
d̂ = 30*(size(MatrixConstruct.H, 2)|>ones)

model_mp = Model(HiGHS.Optimizer)
mp = MP(model_mp, M)

PortOutVariable!(mp, :d) do d
    fix.(d, d̂, force=true)
end

Solve!(mp)
q0 = Getq(mp)
w0 = Getw(mp)

for II in 1:1
    # MSP Initialize Primary Decision Variables
    if II == 1
        model_msp = Model(
            optimizer_with_attributes(HiGHS.Optimizer, "output_flag" =>true),
        )

        global msp = MSP(
            model_msp, 
            d̂,
            M
        )
    end

    Solve!(msp)
    if isnan(objective_value(msp))
        error("Master Problem Infeasible. ") 
    end
    global w̄ = Getw(msp)
    global γ̄ = GetGamma(msp)

    @info "MSP, γ̄, w̄, created for the primary system. " 
    
    # FMP, Initialize a NEW FMP each time for inner CCGA. 
    model_fmp = Model(Gurobi.Optimizer)
    global fmp = FMP(w̄, γ̄, d̂, model_fmp)
    
    @info "FMP, instance is initialized. "
    
    
    for III in 1:1
        
        if III != 1
            Introduce!(fmp, q)
        end
        
        Solve!(fmp)
        U = objective_value(fmp)
        if isnan(U)
            @warn "FMP is not feasible."
        end
        global ρ⁺ = GetRhoPlus(fmp)
        global ρ⁻ = GetRhoMinus(fmp)
        global d = GetDemandVertex(fmp)
        
        model_fsp =  Model(
            optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)
        )

        global fsp = FSP(w̄, γ̄, d, model_fsp)
        Solve!(fsp)
        L = objective_value(fsp)
        global u = Getu(fsp)
        global q = Getq(fsp)
        global v = Getv(fsp)
        @info "Upper Bound U: $U, Lower Bound L: $L" 
        
        # Analyze L, U; and Make decisions 
        # if L > 0
        #     IntroduceCut!(msp, u, q, ρ⁺, ρ⁻)
        #     "Cut Introduced to the MSP, restarts inner CCGA. " |>println
        #     break
        # elseif U <= ϵ
        #     "A good robust solution has been found it seems. "|>println
        # else
        #     # Nothing, it's here just for logical closure. 
        # end
        
        Constraints = IntroduceCut!(msp, u, q, ρ⁺, ρ⁻)
        DebugReport(msp, "msp_after_first_cut")
        DebugReport(mp, "main_problem_referece")
    end

end
