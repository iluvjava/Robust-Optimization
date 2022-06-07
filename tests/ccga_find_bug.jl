include("../src/ccga_modeling.jl")
mp = MP()

Solve!(mp)
DebugReport(mp)
@info "Results"


if !isnan(objective_value(mp))
    @info "Primary discrete decision variables:"
    mp.w.|>value.|>println
    @info "Secondary Discreate decision variables:"
    mp.q.|>value.|>println;
    @info "Non-zero slack"
    findall(>(0), mp.v.|>value)
else
    println("Problem is Infeasible. ")
end

