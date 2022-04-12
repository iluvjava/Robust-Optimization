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


cm |> GetMatrix|> Matrix
cm |> GetMatrix|> Matrix