function Foo()
    x::Int = -1
    setx(x_) = (x = x_) 
    getx() = x
    () -> (getx;setx)
end

f = Foo()
f.getx()
f.setx(22)
f.getx()
f.setx(1)
f.getx()