### Package the CCGA Structs into a new module. 

module CCGA

    include("utilities.jl")
    include("matrix_construction_export.jl")
    # using Logging
    include("ccga_modeling.jl")

    # Read the reports generated from the average demand suggest file.
    # We use the conent of that file to hardcode some of the parameters
    # into the CCGA initialization parameters. 


    ϵ = 0.1
    model_msp = Model(
        optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false),
    )
    msp = MSP(
        model_msp, 
        40*(size(CCGA.MatrixConstruct.H, 2)|>ones), 
        40
    )
    Solve!(msp)
    w̄ = Getw(msp)
    
    for II in 1:10

    end

end