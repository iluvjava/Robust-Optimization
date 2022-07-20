function GenericPortOut(f::Function, a, b)
return f(a, b) end

GenericPortOut(1, 2) do a, b
return a + b end