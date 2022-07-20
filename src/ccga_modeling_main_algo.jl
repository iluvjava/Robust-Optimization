### Package the CCGA Structs into a new module. 

include("utilities.jl")
include("matrix_construction_export.jl")
include("ccga_modeling.jl")

# Read the reports generated from the average demand suggest file.
# We use the conent of that file to hardcode some of the parameters
# into the CCGA initialization parameters. 

ϵ = 0.1
M = 5
d̂ = 40*(size(MatrixConstruct.H, 2)|>ones)
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
    w̄ = w0
    γ̄ = GetGamma(msp)

    "MSP, γ̄, w̄, created for the primary system. " |> println

    # FMP, Initialize a NEW FMP each time for inner CCGA. 
    model_fmp = Model(
        optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)
    )
    global fmp = FMP(w̄, γ̄, d̂, model_fmp)
    "FMP, instance is initialized. " |> println
    let q = u = ρ⁺ = ρ⁻ = d = nothing; for III in 1:3
        if III != 1
            IntroduceVariables!(fmp, q)
            PrepareConstraints!(fmp)
        end
        
        Solve!(fmp)
        U = objective_value(fmp)
        ρ⁺ = GetRhoPlus(fmp)
        ρ⁻ = GetRhoMinus(fmp)
        d = GetDemandVertex(fmp)
        "Upper Bound U: $U" |> println
        model_fsp =  Model(
            optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)
        )
        global fsp = FSP(w̄, γ̄, d, model_fsp)
        Solve!(fsp)
        L = objective_value(fsp)
        u = Getu(fsp)
        q = Getq(fsp)
        "Lower Bound L: $L" |> println

        # Analyze L, U; and Make decisions 
        if L > 0
            IntroduceCut!(msp, u, q, ρ⁺, ρ⁻)
            "Cut Introduced to the MSP, restarts inner CCGA. " |>println
            break
        elseif U <= ϵ
            "A good robust solution has been found it seems. "|>println
        else
            # Nothing, it's here just for logical closure. 
        end
    end end

end
