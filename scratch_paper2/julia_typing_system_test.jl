function f(a::T1, b::T2) where {T1, T2 <: Number}
return a + b end

f(1, 1.2) |> println
f(2, Complex(1, 1)) |> println 

