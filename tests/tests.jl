using Test

@testset "SANITY TEST" begin 

    function sanity_test()
        return true
    end


    @test sanity_test()
end