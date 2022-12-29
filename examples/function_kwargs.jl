function fn(a, b; kwargs...)
    return fn2(a, b; kwargs...)
end

function fn2(a, b; c=0)
    return a, b, c
end

fn(1, 1, c=1)