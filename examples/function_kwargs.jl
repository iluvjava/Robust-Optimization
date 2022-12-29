function fn1(a, b; kwargs...)
    fn2(a, b; kwargs...)
    fn3(a, b; kwargs...)
end

function fn2(a, b;c=0, kwargs...)
    println(a, b, c)
end

function fn3(a, b; d=0, kwargs...)
    println(a, b, d)
end

fn1(1, 1, a=1, b=1, c=1)


