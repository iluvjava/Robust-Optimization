using Test

include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")


@testset "Testing FMPH related functionalities" begin 
    
    function FMPHSetup()

        return true
    end

    @test FMPHSetup()
end
