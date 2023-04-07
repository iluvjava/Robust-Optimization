# Enumeration set Cardinality.
T = CONST_PROBLEM_PARAMETERS.HORIZON
N = PRIMARY_GENERATORS|>length
S = STORAGE_SYSTEM|>length
L = DEMAND_RESPONSE.level

# For matrix C
c = VariableCoefficientHolder(:c, N, T)
p = VariableCoefficientHolder(:p, N, T)
h = VariableCoefficientHolder(:h, S, T)
g_plus = VariableCoefficientHolder(:g_plus, S, T)
g_minus = VariableCoefficientHolder(:g_minus, S, T)
dr = VariableCoefficientHolder(:dr, T)
μ = VariableCoefficientHolder(:μ, L, T)

# For matrix B
x = VariableCoefficientHolder(:x, N, T)
y = VariableCoefficientHolder(:y, N, T)
z = VariableCoefficientHolder(:z, N, T)

# For matrix G
σ⁺ = VariableCoefficientHolder(:σ⁺, S, T)
σ⁻ = VariableCoefficientHolder(:σ⁻, S, T)

# For matrix H
d = VariableCoefficientHolder(:d, T)

# All Matrices: 
B = CoefficientMatrix()
C = CoefficientMatrix()
G = CoefficientMatrix()
"""
Actually the negative `H` matrix in the paper, it's now on the left hand side and takes positive sign.
It's a skinny matrix with mostly zeros on the top and negative identity matrix 
on the bottom few rows for the coefficients of d, for demand balance constraints. 
"""
H = CoefficientMatrix()
B(x, y, z)
C(
    c,
    p,
    h,
    g_plus,
    g_minus,
    dr,
    μ
)
G(σ⁺, σ⁻)
H(d)

Sync() = SyncRow(B, C, G, H)



"""
Fuel constraints. 
"""
function FuelConstraints()
    pg = PRIMARY_GENERATORS
    ρ = DEMAND_RESPONSE.rho
    rhs = Vector{Number}()
    
    for n in 1:N, t in 1:T
        x[n, t] = pg.SU[n]
        z[n, t] = pg.SD[n]
        c[n, t] = 1
    end
    for l in 1:L, t in 1:T
        μ[l, t] = ρ[l]
    end
    C(c, μ); C()
    B(x, z); B()
    push!(rhs, CONST_PROBLEM_PARAMETERS.Φ)
    Sync()

    for n in 1:N, t in 1:T, k in 1:(size(pg.alphas, 2))
        p[n, t] = pg.alphas[n, k]
        y[n, t] = pg.betas[n, k]
        c[n, t] = -1
        C(p, c); B(y)
        C(); B()
        push!(rhs, 0)
    end
    Sync()
    
return rhs end


function CapacityContraints()
    gen = PRIMARY_GENERATORS
    rhs = Vector{Number}()
    for t = 2:T, n = 1:(gen|>length)
        # (13)
        p[n, t] = 1
        p[n, t - 1] = -1
        y[n, t - 1] = -gen.RU[n]
        x[n, t] = -gen.RU_bar[n]
        C(p); C()
        B(y, x); B()
        push!(rhs, 0)
    end
    Sync()
    for t = 2:T, n = 1:(gen|>length)
        # (14)
        p[n, t - 1] = 1
        p[n, t] = -1
        y[n, t] = -gen.RD[n]
        z[n, t] = - gen.RD_bar[n]
        C(p); C()
        B(y, z); B()
        push!(rhs, 0)
    end
    Sync()
    for t = 1:T, n = 1:(gen|>length)
        # (16)
        p[n, t] = 1
        y[n, t] = -gen.Pmax[n]
        C(p); C(); B(y); B();
        push!(rhs, 0)
    end
    Sync()
    for t = 1:T, n = 1:(gen|>length)
        # (17)
        p[n, t] = -1
        y[n, t] = gen.Pmin[n]
        C(p); C(); B(y); B();
        push!(rhs, 0)
    end
    Sync()
    return rhs
end

function DemandResponseConstraints()
    rhs = Vector()
    R = DEMAND_RESPONSE.R
    dr_max = DEMAND_RESPONSE.dr_max

    for l in 1:L, t in 1:T
        μ[l, t] = -1
        push!(rhs, 0)
        C(μ); C()
    end
    Sync()

    for l in 1:L, t in 1:T
        μ[l, t] = 1
        push!(rhs, R[l] - R[l - 1])
        C(μ); C()
    end
    Sync()

    for t in 1: T
        for l in 1:L
            μ[l, t] = -1
        end
        dr[t] = 1
        push!(rhs, 0)
        C(μ, dr); C() 
    end
    for t in 1: T
        for l in 1:L
            μ[l, t] = 1
        end
        dr[t] = -1
        push!(rhs, 0)
        C(μ, dr); C() 
    end
    Sync()

    for t in 1:T
        dr[t] = 1
        push!(rhs, dr_max)
        C(dr); C()
    end
    Sync()

    return rhs
end


"""
Battery constraints, constraints 34 ... 38
"""
function BatteryConstraints()
    rhs = Vector{Number}()
    v = STORAGE_SYSTEM.Efficiency
    H̄ = STORAGE_SYSTEM.Capacity
    G_plus = STORAGE_SYSTEM.CharingLim
    G_minus = STORAGE_SYSTEM.DischargingLim

    for s in 1: S
        # Note, the h[0, s] has been set to a constant value: 0, manually. 
        h[s, 1] = 1
        C(h); C()
        push!(rhs, 0)

        h[s, 1] = -1
        C(h); C()
        push!(rhs, 0)
        
    end
    
    for t in 2:T, s in 1:S
        h[s, t] = 1
        h[s, t - 1] = -1
        g_plus[s, t - 1] = -v[s]
        g_minus[s, t - 1] = v[s]
        C(h, g_plus, g_minus); C()
        push!(rhs, 0)
    end
    
    for t in 2:T, s in 1:S
        h[s, t] = -1
        h[s, t - 1] = 1
        g_plus[s, t - 1] = v[s]
        g_minus[s, t - 1] = -v[s]
        C(h, g_plus, g_minus); C()
        push!(rhs, 0)
    end

    for t in 1:T, s in 1:S
        h[s, t] = 1
        C(h); C()
        push!(rhs, H̄[s])
    end

    Sync()

    for t in 1:T, s in 1:S
        g_plus[s, t] = 1
        σ⁺[s, t] = -G_plus[s]
        C(g_plus)
        G(σ⁺)
        C()
        G()
        push!(rhs, 0)
    end
    Sync()

    for t in 1:T, s in 1:S
        g_minus[s, t] = 1
        σ⁻[s, t] = -G_minus[s]
        C(g_plus)
        G(σ⁻)
        C()
        G()
        push!(rhs, 0)
    end
    Sync()

    for t in 1:T, s in 1:S
        σ⁺[s, t] = 1
        σ⁻[s, t] = 1
        G(σ⁺, σ⁻)
        G()
        push!(rhs, 1)
        σ⁺[s, t] = -1
        σ⁻[s, t] = -1
        G(σ⁺, σ⁻)
        G()
        push!(rhs, -1)
    end
    Sync()

return rhs end

function DemandsBalanceConstraints() 
    rhs = Vector{Number}()
    for t in 1:T
        for n in 1: N
            p[n, t] = 1
        end
        for s in 1:S
            g_plus[s, t] = -1
            g_minus[s, t] = 1
        end
        dr[t] = 1
        d[t] = -1
        C(p, g_plus, g_minus, dr); C()
        H(d); H()
        push!(rhs, 0)
    end

    for t in 1:T
        for n in 1: N
            p[n, t] = -1
        end
        for s in 1:S
            g_plus[s, t] = 1
            g_minus[s, t] = -1
        end
        dr[t] = -1
        d[t] = 1
        C(p, g_plus, g_minus, dr); C()
        H(d); H()
        push!(rhs, 0)
    end

    Sync()

    return rhs
end


"""
Given a symbol and the matrix that is supposed to responsbible for it, it returns 
what range of columns of the matrix corresponds to the given variable. 

The columns will start indexing at one. 
"""
function ColumnRegimeFor(matrix::CoefficientMatrix, var::VariableCoefficientHolder)
    if !(var.v in matrix.var_posi|>keys)
        return -1 
    end
    starting = matrix.var_posi[var.v]
return (starting + 1, starting + length(var)) end

rhs1 = FuelConstraints()
rhs2 = CapacityContraints()
rhs3 = DemandResponseConstraints()
rhs4 = BatteryConstraints()
rhs5 = DemandsBalanceConstraints()

B_ = B
C_ = C
G_ = G
H_ = G

rhs = vcat(rhs1, rhs2, rhs3, rhs4, rhs5)
B, C, G, H = [B, C, G, H].|>GetMatrix

w = [x, y, z]
u = [c, p, h, g_plus, g_minus, dr, μ]
q = [σ⁺, σ⁻]
# d is it's own vector.

@info "Matrices for robust optimizations successfully constructed. "