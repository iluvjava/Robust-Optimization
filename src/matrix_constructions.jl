# CONTEXT OF THE CODE ===================================================================================================
# Requires, problem_parameters.jl, coeff-mgr.jl, util.jl, prior to loading. 
#   The code that generated the first stage constraints problem for the MSP problem is legacy code hanled in ccga_modeling_mp.jl
#   separately. 
# The code is simple and it does this: 
#   - Assign coefficient for a varibale. 
#   - Call it with the matrix to add the coefficientto a specific row, i.e: `C(x)`, 
#   - Make a new line on the coefficient matrix, i.e: `C()`, crates a new row on the sparse matrix. 
#   - Sync the row for all matrices, so that the current row is the same for all matrix. 
# This code creates the coefficient matrix for the abstract form in the paper. 
# ===================================================================================================



# Enumeration set Cardinality.
global T = CONST_PROBLEM_PARAMETERS.HORIZON
global N = PRIMARY_GENERATORS|>length
global S = STORAGE_SYSTEM|>length
global L = DEMAND_RESPONSE.level

function EstablishMatrices()

    # For matrix C
    global c = VariableCoefficientHolder(:c, N, T)
    global p = VariableCoefficientHolder(:p, N, T)
    global h = VariableCoefficientHolder(:h, S, T)
    global g_plus = VariableCoefficientHolder(:g_plus, S, T)
    global g_minus = VariableCoefficientHolder(:g_minus, S, T)
    global dr = VariableCoefficientHolder(:dr, T)
    global μ = VariableCoefficientHolder(:μ, L, T)

    # For matrix B
    global x = VariableCoefficientHolder(:x, N, T)
    global y = VariableCoefficientHolder(:y, N, T)
    global z = VariableCoefficientHolder(:z, N, T)

    # For matrix G
    global σ⁺ = VariableCoefficientHolder(:σ⁺, S, T)
    global σ⁻ = VariableCoefficientHolder(:σ⁻, S, T)

    # For matrix H
    global d = VariableCoefficientHolder(:d, T)

    # All Matrices: 
    global B = CoefficientMatrix()
    global C = CoefficientMatrix()
    global G = CoefficientMatrix()
    """
    Actually the negative `H` matrix in the paper, it's now on the left hand side and takes positive sign.
    It's a skinny matrix with mostly zeros on the top and negative identity matrix 
    on the bottom few rows for the coefficients of d, for demand balance constraints. 
    """
    global H = CoefficientMatrix()
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
        # [x]: Add initial conditions here. 
        for t = 1:1, n = 1:(gen|>length)
            p[n, t] = 1
            x[n, t] = -gen.RU_bar[n]
            C(p); C()
            B(y, x); B()
            push!(rhs, -gen.RU[n]*gen.initial_status[n] + gen.initial_pg[n])
        end
        Sync()
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
        # [x]: Change here, on the value of p, the initial condition. 
        for t = 1:1, n = 1:(gen|>length)
            p[n, t] = -1
            y[n, t] = -gen.RD[n]
            z[n, t] = -gen.RD_bar[n]
            C(p); C()
            B(y, z); B()
            push!(rhs, gen.RD[n]*gen.initial_status[n] - gen.initial_pg[n])
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

    rhs1 = FuelConstraints()
    rhs2 = CapacityContraints()
    rhs3 = DemandResponseConstraints()
    rhs4 = BatteryConstraints()
    rhs5 = DemandsBalanceConstraints()

    global B_ = B
    global C_ = C
    global G_ = G
    global H_ = G

    global rhs = vcat(rhs1, rhs2, rhs3, rhs4, rhs5)
    global B, C, G, H = [B, C, G, H].|>GetMatrix

    global w = [x, y, z]
    global u = [c, p, h, g_plus, g_minus, dr, μ]
    global q = [σ⁺, σ⁻]
    # d is it's own vector.


end

EstablishMatrices()


"""
Given a symbol and the matrix that is supposed to responsbible for it, it returns 
what range of columns of the matrix corresponds to the given variable. 

The columns will start indexing at one. 

WARN: i think this function should be in `coef_mgr_v2.jl` instead. 
"""
function ColumnRegimeFor(matrix::CoefficientMatrix, var::VariableCoefficientHolder)
    if !(var.v in matrix.var_posi|>keys)
        return -1 
    end
    starting = matrix.var_posi[var.v]
return (starting + 1, starting + length(var)) end

@info "Matrices for robust optimizations successfully constructed. "


