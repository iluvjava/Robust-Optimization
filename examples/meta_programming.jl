using Mixers

@premix mutable struct YourPoint
    x::Int
    y::Int
end

@premix mutable struct HerPoint
    z::Int
end

@macroexpand @premix mutable struct HerPoint
    z::Int
end

@HerPoint @YourPoint mutable struct MyPoint end

macro show_value(variable::Symbol)
    quote
        println("The ", $(string(variable)), " you passed is ", $(esc(variable)))
    end
end
