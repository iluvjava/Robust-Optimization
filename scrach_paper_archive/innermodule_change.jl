module SomeModule

    function ChangeIt(t)
        global VAR = t
    end
    global G, H = [1, 2]
end

SomeModule.ChangeIt(2)

println("SomeModule.VAR: $(SomeModule.VAR)")
println("SomeModule.G: $(SomeModule.VAR)")
println("SomeModule.H: $(SomeModule.VAR)")