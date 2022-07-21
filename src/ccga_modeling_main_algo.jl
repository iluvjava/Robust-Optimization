### READ CCGA.md For more info. 

include("utilities.jl")
include("matrix_construction_export.jl")
include("ccga_modeling.jl")

ϵ = 0.1
M = 30
d̂ = 30*(size(MatrixConstruct.H, 2)|>ones)
model_mp = Model(HiGHS.Optimizer)
mp = MP(model_mp, M)
PortOutVariable!(mp, :d) do d
    fix.(d, d̂, force=true)
end

Solve!(mp)
q0 = Getq(mp)
w0 = Getw(mp)

# for II in 1:1
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
    w̄ = Getw(msp)
    γ̄ = GetGamma(msp)

    @info "MSP, γ̄, w̄, created for the primary system. " 
    
    # FMP, Initialize a NEW FMP each time for inner CCGA. 
    model_fmp = Model(Gurobi.Optimizer)
    global fmp = FMP(w̄, γ̄, d̂, model_fmp)
    
    @info "FMP, instance is initialized. "
    
    let q = u = ρ⁺ = ρ⁻ = d = nothing; for III in 1:1
        
        if III != 1
            IntroduceVariables!(fmp, q)
            PrepareConstraints!(fmp)
        end
        
        Solve!(fmp)
        U = objective_value(fmp)
        ρ⁺ = GetRhoPlus(fmp)
        ρ⁻ = GetRhoMinus(fmp)
        d = GetDemandVertex(fmp)
        
        model_fsp =  Model(
            optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)
        )

        global fsp = FSP(w̄, γ̄, d, model_fsp)
        Solve!(fsp)
        L = objective_value(fsp)
        u = Getu(fsp)
        q = Getq(fsp)
        v = Getv(fsp)
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
        DebugReport(msp, "msp_all_constraints")
        DebugReport(mp, "main_problem_constraints_sets")

    end end

# end
