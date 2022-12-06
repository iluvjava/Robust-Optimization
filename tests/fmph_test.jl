using Infiltrator, ProgressMeter, Dates
include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")

N = 3*4*6
fmph1 = FMPH1(zeros(N), 10*ones(24), 100*ones(24))
Solve!(fmph1)
# IntroduceCut!(fmph1, zeros(size(MatrixConstruct.G, 2)))
Solve!(fmph1)
fmph2 = FMPH2(zeros(N), 10*ones(24), 100*ones(24), fmph1.lambda[1].|>value)
Solve!(fmph2)

