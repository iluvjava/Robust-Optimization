include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")

# using Logging

M = 30
d̂ = 40*(size(MatrixConstruct.H, 2)|>ones)

# model_mp =  Model(
#     optimizer_with_attributes(
#             HiGHS.Optimizer, 
#             "output_flag" =>true, 
#             "mip_feasibility_tolerance"=>1e-8, 
#             "ipm_optimality_tolerance" => 1e-6, 
#             "dual_feasibility_tolerance" => 1e-8, 
#             "primal_feasibility_tolerance" => 1e-8, 
#             "mip_rel_gap" => 0.0001,
#             "mip_abs_gap" => 1e-04
#         ) 
#     )

model_mp = Model(
    optimizer_with_attributes(
       Gurobi.Optimizer, "BarConvTol" => 0
    )
)

mp = MP(model_mp, 30)
PortOutVariable!(mp, :d) do d
    fix.(d, d̂; force=true)
end
Solve!(mp)
u = Getu(mp)
q = Getq(mp)
ρ⁺ = ones(d̂ |> length)
ρ⁻ = zeros(d̂ |> length)

# model_msp =  Model(optimizer_with_attributes(HiGHS.Optimizer, "output_flag" =>true))
model_msp = Model(Gurobi.Optimizer)
msp = MSP(model_msp, d̂, M)
IntroduceCut!(msp, u, q, ρ⁺, ρ⁻)


function TestConstraints()
    B = MatrixConstruct.B
    h = MatrixConstruct.h
    G = MatrixConstruct.G
    C = MatrixConstruct.C
    H = MatrixConstruct.H
    w = Getw(mp)
    d = mp.d.|>value
    v = mp.v.|>value
    res = h - (B*w + C*u + G*q + H*d - v)
    res = map((x) -> (x < 0 ? x : 0), res)

return res end

res = TestConstraints()
