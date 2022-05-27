
mutable struct Foo
    x::Int
    
    function Foo(x)
    return new(x) end
    
    function Get(this::Foo)
    return this.x end

    () -> (Get)
end


foo = Foo(10) 