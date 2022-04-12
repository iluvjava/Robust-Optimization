include("coef_mgr_v2.jl")

x = VariableCoefficientHolder(:x, 3,3)
y = VariableCoefficientHolder(:y, 3,3)
z = VariableCoefficientHolder(:z, 2,3)
cm = CoefficientMatrix()

x[1, 1] = 1
y[1, 1] = -1
cm(x, y)  # add coefficient to current row
!cm       # next row on the matrix!

x[2,2] = 1
y[2,2] = 1
cm(x, y)  # add coefficint of x, y to that row 
!cm       # NEXT ROW!

z[:, 1].= 1
z[:, 2].= -1
cm(z)

A = cm |> GetMatrix

# Now we make the jump model variables in the same order of how we constructed the 
# matrix 
jmodel= Model()
x = @variable(jmodel, x[1:3,1:3])
y = @variable(jmodel, y[1:3,1:3])
z = @variable(jmodel, z[1:2,1:3])
@constraint(jmodel, A*vcat(x[:], y[:], z[:]) .== 0)