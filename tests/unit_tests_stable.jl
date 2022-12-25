### The point of unit tests is not to catch any bug, but a preliminary verification that we 
### didn't break the code while adding more code to existing code. 
using Test

include("../src/utilities.jl")
include("../src/matrix_construction_export.jl")
include("../src/ccga_modeling.jl")


@testset "Testing Basics" begin 
"""
The minimal amount of code to test to make sure the newly introduce code are not breaking the old 
code. 
"""
    C = MatrixConstruct.C
    B = MatrixConstruct.B
    H = MatrixConstruct.H
    G = MatrixConstruct.G
    h = MatrixConstruct.h
    u = MatrixConstruct.u
    u_len = u.|>length|>sum
    w = MatrixConstruct.w
    w_len = w.|>length|>sum
    q = MatrixConstruct.q
    q_len = q.|>length|>sum

    """
    Verifying that the matrices: C,B,H,G and the vector u, w, q are all having 
    compatible dimensons with each other. 
    """
    function VerifyMatrices()
        @assert w_len == size(B, 2) "vector w and matrix B has incompatible size. "
        @assert u_len == size(C, 2) "Vector u and matrix C has incompatible size. "
        @assert q_len == size(G, 2) "Vector q andmatrix G has incompatible size."
        @assert size(C, 1) == size(B, 1) == size(G, 1) == size(H, 1) "The number of rows for the matrices: C, B, G, H are not equal. "*
        "Currently they are: C=$(size(C, 1)), B=$(size(B, 1)), G=$(size(G, 1)), H=$(size(H, 1)). "
        return true 
    end

    """
    Using the main problem to solve for a configurations to establish the variable 
    w̄, γ̄ for the system. 
    """
    function MainProblem()
        mp = MainProblem()
        return true
    end

    # actually running these tests. 
    @test VerifyMatrices()
end
