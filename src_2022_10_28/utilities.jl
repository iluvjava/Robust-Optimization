@static if Sys.islinux()
    @error "I don't know what to do if it's linux"
elseif Sys.isapple()
    ENV["GUROBI_HOME"] = "/Library/gurobi950/mac64"
elseif Sys.iswindows()
    ENV["GUROBI_HOME"] = "C:\\gurobi952\\win64"
end
using SparseArrays, LinearAlgebra, JuMP, GLPK, HiGHS, Gurobi, Plots
