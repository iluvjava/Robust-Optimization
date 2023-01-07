function Subroutine1()
    return println("Subroutine one")
end

function Subroutine2()
    return println("Subroutine two")
end

function Dispatcher(routine::typeof(Subroutine1))
    routine()
    println("$(routine|>Symbol|>String) dispatched. ")
end

function Dispatcher(routine::typeof(Subroutine2))
    routine()
    println("$(routine|>Symbol|>String) dispatched. ")
end

Dispatcher(Subroutine1)
Dispatcher(Subroutine2)