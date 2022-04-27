### This file focuses on the constructions of matrices. 

# include("coef_mgr_v2.jl")


include("coef_mgr_v2.jl")
include("problem_parameters.jl")

# Enumeration set Cardinality. 
T = CONST_PROBLEM_PARAMETERS.HORIZON
N = PRIMARY_GENERATORS|>length
M = SECONDARY_GENERATORS|>length
S = STORAGE_SYSTEM|>length

x = VariableCoefficientHolder(:x, N, T)
y = VariableCoefficientHolder(:y, N, T)
z = VariableCoefficientHolder(:z, N, T)
x′ = VariableCoefficientHolder(:x′, N, T)
y′ = VariableCoefficientHolder(:y′, N, T)
z′ = VariableCoefficientHolder(:z′, N, T)

c = VariableCoefficientHolder(:c, N, T)
c′ = VariableCoefficientHolder(:c′, M, T)
p = VariableCoefficientHolder(:p, N, T)
p′ = VariableCoefficientHolder(:p′, M, T)
regu = VariableCoefficientHolder(:regu, N, T)
regu′ = VariableCoefficientHolder(:regu′, M, T)
regd = VariableCoefficientHolder(:regd, N, T)
regd′ = VariableCoefficientHolder(:regd′, M, T)
sr = VariableCoefficientHolder(:sr, N, T)
sr′ = VariableCoefficientHolder(:sr′, M, T)
h = VariableCoefficientHolder(:h, S, T)
g_plus = VariableCoefficientHolder(:g_plus, S, T)
g_minus = VariableCoefficientHolder(:g_minus, S, T)
nsp = VariableCoefficientHolder(:nsp, N, T)
nsp′ = VariableCoefficientHolder(:nps′, M, T)

function MakeCoefMatrices()
    B = CoefficientMatrix()
    C = CoefficientMatrix()
    G = CoefficientMatrix()
    B(x, y, z)
    C(
        c, 
        c′, 
        p, 
        p′, 
        regu, 
        regu′, 
        regd, 
        regd′, 
        sr, 
        sr′, 
        h, 
        g_plus, 
        g_minus, 
        nsp, 
        nsp′
    )
    G(x′, y′, z′)
return B, C, G end

B, C, G = MakeCoefMatrices()

## All global variable, let's not worry about the bad programming for now, 
# improve them later. 


"""
    Adding the first row, the fuel constraints. RHS is returned by the function, as 
    a vector. 
"""
function FuelConstraints()
    pg = PRIMARY_GENERATORS
    sg = SECONDARY_GENERATORS
    rhs = Vector{Number}()
    for m in 1:M, t in 1:T
        x′[m, t] = sg.SU[m]
        z′[m, t] = sg.SD[m]
        c′[m, t] = 1
    end
    for n in 1:N, t in 1:T
        x[n, t] = pg.SU[n]
        z[n, t] = pg.SD[n]
        c[n, t] = 1
    end
    G(x′, z′); G()
    C(c′, c); C()
    B(x, z); B()
    
    push!(rhs, CONST_PROBLEM_PARAMETERS.Φ)

    for n in 1:N, t in 1:T, k in 1:(pg.alphas|>length)
        p[n, t] = pg.alphas[k]
        y[n, t] = pg.betas[k]
        c[n, t] = -1
        C(p, c); B(y)
        C(); B()
        push!(rhs, 0)
    end

    for m in 1:M, t in 1:T, k in 1:(sg.alphas|>length)
        p′[m, t] = sg.alphas[k]
        y′[m, t] = sg.betas[k]
        c′[m, t] = -1
        C(p′, c′); B(y′)
        C(); B()
        push!(rhs, 0)
    end
    # SYNC MATRICES
    
    SyncRow(B, C, G)
return rhs end

function QuickStartConstraints()
    rhs = Vector{Number}()
    

return rhs end



# ------------------------------------------------------------------------------
# CALL these construction functions in the correct oder 
FuelRHS = FuelConstraints()
