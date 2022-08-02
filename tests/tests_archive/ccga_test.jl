include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")

# ======================================================================================================================
# A basic cut performed by using the main problem. 
# ======================================================================================================================

M = 10
d̂ = 40*(size(MatrixConstruct.H, 2)|>ones)


model_mp = Model(
    optimizer_with_attributes(
       Gurobi.Optimizer, "BarConvTol" => 0
    )
)

mp = MP(model_mp, M)
PortOutVariable!(mp, :d) do d
    fix.(d, d̂; force=true)
end

Solve!(mp)
u = Getu(mp)
q = Getq(mp)
ρ⁺ = ones(d̂ |> length)
ρ⁻ = zeros(d̂ |> length)
v = Getv(mp)

model_msp = Model(Gurobi.Optimizer)
msp = MSP(model_msp, d̂, M)
IntroduceCut!(msp, u, q, ρ⁺, ρ⁻, v)

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

# ======================================================================================================================
# All zeros cut. 
# ======================================================================================================================
model_msp = Model(Gurobi.Optimizer)
msp = MSP(model_msp, d̂, M)
IntroduceCut!(msp, 
    zeros(size(MatrixConstruct.C, 2)), 
    zeros(size(MatrixConstruct.G, 2)), 
    zeros(size(MatrixConstruct.H, 2)), 
    zeros(size(MatrixConstruct.H, 2))
)
