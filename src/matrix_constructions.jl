### This file focuses on the constructions of matrices. 

# Enumeration set Cardinality. 
T = CONST_PROBLEM_PARAMETERS.HORIZON
N = PRIMARY_GENERATORS|>length
M = SECONDARY_GENERATORS|>length
S = STORAGE_SYSTEM|>length
B̄ = (SIGMAS|>size)[1]
L = (SIGMAS|>size)[2]

# ==============================================================================
# set up all decisions variables and their dimensions using the enumeration sets
# cardinality. 
# ==============================================================================
# For matrix B
x = VariableCoefficientHolder(:x, N, T)
y = VariableCoefficientHolder(:y, N, T)
z = VariableCoefficientHolder(:z, N, T)

# For matrix G
x′ = VariableCoefficientHolder(:x′, N, T)
y′ = VariableCoefficientHolder(:y′, N, T)
z′ = VariableCoefficientHolder(:z′, N, T)

# For Matrix C
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

# For matrix F
d = VariableCoefficientHolder(:d, B̄, T)

function MakeCoefMatrices()
    B = CoefficientMatrix()
    C = CoefficientMatrix()
    G = CoefficientMatrix()
    F = CoefficientMatrix()
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
    F(d)
    G(x′, y′, z′)
return B, C, G, F end


# All global variables, let's not worry about the bad programming for now, 
# improve them later.

# ==============================================================================
# Constraints adding functions
# ==============================================================================

"""
    Adding the first row, the fuel constraints. RHS is returned by the function, as 
    a vector. 
    constraints (6, 7, 8)
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
        y[m, t] = sg.betas[k]
        c′[m, t] = -1
        C(p′, c′); B(y)
        C(); B()
        push!(rhs, 0)
    end

    
return rhs end

"""
    Adding the quick start binary constraints 
    Constraints (9, ..., 12)
"""
function QuickStartConstraints()
    sg = SECONDARY_GENERATORS
    rhs = Vector{Number}()
    for t = 1:T, m = 1:M
        # (9)
        x′[m, t] = 1
        y′[m, t] = 1
        G(x′, y′)
        G()
        push!(rhs, 1)
    end

    for t = 2:T, m = 1:M  # note, t starts with 2. 
        # (10)
        y′[m, t] = 1; y′[m ,t - 1] = -1
        x′[m ,t] = -1; z′[m, t] = 1
        G(x′, y′, z′); G()
        push!(rhs, 0)
        y′[m, t] = -1; y′[m ,t - 1] = 1
        x′[m ,t] = 1; z′[m, t] = -1
        G(x′, y′, z′); G()
        push!(rhs, 0)
    end

    for t = 1:T, m = 1:M
        # (11)
        for tau = t - sg.Tminu[m] + 1
            x′[m, tau] = 1
        end
        y′[m, t] = -1
        G(x′, y′); G()
        push!(rhs, 0)
    end

    for t = 1:T, m = 1:M
        # (11)
        for tau = t - sg.Tmind[m] + 1
            z′[m, tau] = 1
        end
        y′[m, t] = 1
        G(x′, y′); G()
        push!(rhs, 1)
    end

return rhs end


"""
    The capacity constraints of the first stage generator. 
    constraints (13, ..., 20)
    * Pass in the generator instance, matrix C, or Matrix G to specify 
    whether these sets of constraints are for primary, or secondary generators. 
"""
function CapacityConstraints(
    gen::Generators, 
    K::CoefficientMatrix,
    p::VariableCoefficientHolder, 
    x::VariableCoefficientHolder, 
    y::VariableCoefficientHolder, 
    z::VariableCoefficientHolder,
    sr::VariableCoefficientHolder, 
    regd::VariableCoefficientHolder, 
    regu::VariableCoefficientHolder,
    nsp::VariableCoefficientHolder
)
    rhs = Vector{Number}()

    for t = 2:T, n = 1:(gen|>length)
        # (13)
        p[n, t] = 1
        p[n, t - 1] = -1
        y[n, t - 1] = gen.RU[n]
        x[n, t] = - gen.RU_bar[n]
        C(p); C()
        K(y, x); K()
        push!(rhs, 0)
    end
    for t = 2:T, n = 1:(gen|>length)
        # (14)
        p[n, t - 1] = 1
        p[n, t] = -1
        y[n, t] = -gen.RD[n]
        z[n, t] = - gen.RD_bar[n]
        C(p); C()
        K(y, z); K()
        push!(rhs, 0)
    end
    for t = 1:T, n = 1:(gen|>length)
        # (15)
        sr[n, t] = 1
        y[n, t] = -gen.SR[n]
        C(sr); C(); K(y); K();
        push!(rhs, 0)
    end
    for t = 1:T, n = 1:(gen|>length)
        # (16)
        p[n, t] = 1; sr[n, t]=1; regu[n, t] = 1; y[n, t] = -gen.Pmax[n]
        C(p, sr, regu); C(); K(y); K();
        push!(rhs, 0)
    end
    for t = 1:T, n = 1:(gen|>length)
        # (17)
        p[n, t] = -1; regd[n, t] = 1; y[n, t] = gen.Pmin[n]
        C(p, regd); C(); K(y); K();
        push!(rhs, 0)
    end
    for t = 1:T, n = 1:(gen|>length)
        # (18)
        nsp[n, t] = 1; y[n, t] = gen.NSP[n]
        C(nsp); C(); K(y); K();
        push!(rhs, gen.NSP[n])
    end
    for t = 1:T, n = 1:(gen|>length)
        # (19)
        regu[n, t] = 1; y[n, t] = - gen.REGU[n]
        C(regu); C(); K(y); K();
        push!(rhs, 0)
    end
    for t = 1:T, n = 1:(gen|>length)
        # (20)
        regd[n, t] = 1; y[n, t] = - gen.REGD[n]
        C(regd); C(); K(y); K();
        push!(rhs, 0)
    end
    
return rhs end

"""
    Constraints 30 to 33. 
    
"""
function MinimumRequirement()
    rhs = Vector{Number}()
    glb = CONST_PROBLEM_PARAMETERS
    for t = 1: T
        for n = 1:N
            regu[n, t] = -1 
        end
        for m = 1:M
            regu′[m, t] = -1
        end
        C(regu, regu′); C();
        push!(rhs, -glb.RREGU[t])  
    end
    for t = 1: T
        for n = 1:N
            regd[n, t] = -1 
        end
        for m = 1:M
            regd′[m, t]  = -1
        end
        C(regd, regd′); C();
        push!(rhs, -glb.RREGD[t]) 
    end
    for t = 1:T
        for n = 1:N
            nsp[n, t] = -1
        end
        for m = 1:M
            nsp′[m, t] = -1
        end
        C(nsp, nsp′); C();
        push!(rhs, -glb.RNSP[t])
    end

return rhs end


"""
    Batter constraints, constraints 34 ... 38
"""
function BatteryConstraints()
    rhs = Vector{Number}()
    v = STORAGE_SYSTEM.Efficiency
    H̄ = STORAGE_SYSTEM.Capacity
    G_plus = STORAGE_SYSTEM.CharingLim
    G_minus = STORAGE_SYSTEM.DischargingLim
    # (34)
    for t = 1:T - 1, s = 1:S
        h[s, t + 1] = 1
        h[s, t] = -1
        g_plus[s, t] = -v[s]
        g_minus[s, t] = v[s]
        C(h, g_plus, g_minus); C()
        push!(rhs, 0)
    end
    # (35)
    for t = 1:T - 1, s = 1:S
        h[s, t + 1] = -1
        h[s, t] = 1
        g_plus[s, t] = v[s]
        g_minus[s, t] = -v[s]
        C(h, g_plus, g_minus); C()
        push!(rhs, 0)
    end
    # (36)
    for t = 1:T, s = 1:S
        h[s, t] = 1
        C(h); C()
        push!(rhs, H̄[s])
    end
    # (37)
    for t = 1:T, s = 1:S
        g_plus[s, t] = 1
        C(g_plus); C()
        push!(rhs, G_plus[s])
    end
    # (38)
    for t = 1:T, s = 1:S
        g_minus[s, t] = 1
        C(g_minus); C()
        push!(rhs, G_minus[s])
    end

return rhs end


"""
    Constraints (39) to (41)
"""
function DemandBalanceConstraints()
    rhs = Vector{Number}()
    μ = DISFACTORS.the_matrix
    σ = SIGMAS.the_matrix
    f = TRANSMISSION_SYSTEM.Limit
    # (39)
    for t = 1:T
        p[:, t] .= 1
        p′[:, t] .= 1
        g_minus[:, t] .= 1
        g_plus[:, t] .= -1
        d[:, t] .= -1
        C(p, p′, g_minus, g_plus); F(d)
        C();F();
        push!(rhs, 0)
    end
    # (40)
    for t = 1:T
        p[:, t] .= -1
        p′[:, t] .= -1
        g_minus[:, t] .= -1
        g_plus[:, t] .= 1
        d[:, t] .= 1
        C(p, p′, g_minus, g_plus); F(d)
        C();F();
        push!(rhs, 0)
    end
    # (41)
    # currently each bus has one primary generator, and one secondary generator
    # currently all transmission line has the same storage system. 

    for l=1:L, t=1:T
        p[:, t] .= sum(σ[:, l])
        p′[:, t] .= sum(σ[:, l])
        g_minus[:, t] .= sum(μ[:])
        g_plus[:, t] .= -sum(μ[:])
        d[:, t] = -σ[:, l]
        C(p, p′, g_minus, g_plus); C();
        F(d); F();
        push!(rhs, f[l])
    end    

    # (42)
    for l=1:L, t=1:T
        p[:, t] .= -sum(σ[:, l])
        p′[:, t] .= -sum(σ[:, l])
        g_minus[:, t] .= -sum(μ[:])
        g_plus[:, t] .=  sum(μ[:])
        d[:, t] = σ[:, l]
        C(p, p′, g_minus, g_plus); C();
        F(d); F();
        push!(rhs, f[l])
    end

return rhs end



# ------------------------------------------------------------------------------
# CALL these construction functions in the correct oder 
# Visualize the matrices for debugging purpose. 


B, C, G, F = MakeCoefMatrices()

# Fuel Constraints
# RHS = Vector{Float64}()
RHS = FuelConstraints()
SyncRow(B, C, G, F)

# Quick start Constraints
RHS = vcat(RHS, QuickStartConstraints())
SyncRow(B, C, G, F)

# Capacity constraints primary
RHS = vcat(
    RHS, 
    CapacityConstraints(
        PRIMARY_GENERATORS,
        B,
        p,
        x,
        y,
        z,
        sr,
        regd,
        regu,
        nsp
    )
)
SyncRow(B, C, G, F)

# Capacity constraints secondary
RHS = vcat(
    RHS, 
    CapacityConstraints(
        SECONDARY_GENERATORS,
        G,
        p′,
        x′,
        y′,
        z′,
        sr′,
        regd′,
        regu′,
        nsp′
    )
)
SyncRow(B, C, G, F)

# minimum requirement constraints
RHS = vcat(RHS, MinimumRequirement())
SyncRow(B, C, G, F)

# battery constraints 
RHS = vcat(RHS, BatteryConstraints())
SyncRow(B, C, G, F)

# demand balance
RHS = vcat(RHS, DemandBalanceConstraints())
SyncRow(B, C, G, F)

B = B|>GetMatrix
C = C|>GetMatrix
G = G|>GetMatrix
H = F|>GetMatrix

## Start grouping the cofficient holder for each of the groupped decision variables

w = [x, y, z]
u = [c, c′, p, p′, regu, regu′, regd, regd′, sr, sr′, h, g_plus, g_minus, nsp, nsp′]
q = [x′, y′, z′]
# d is it's own vector. 


