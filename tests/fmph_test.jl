using Infiltrator, ProgressMeter, Dates
include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")

N = 3*4*6
fmph = FMPH1(zeros(N), 20*ones(24), 20*ones(24))

