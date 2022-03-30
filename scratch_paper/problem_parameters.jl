struct CONST
    HORIZON 
    GENERATOR
    QS
    BUSES
    ARCS 
    BREAKPOINTS
    Φ
    function CONST()
        this = new()
        this.Φ = 1000000 #budget
        this.HORIZON = 3 #time HORIZON 
        this.HORIZON #same thing as horizon but more consise notation
        this.GENERATORS = 6 #number of generators
        this.QS = 6 #number of quick start generators
        this.BUSES = 6 #number of BUSES
        this.ARCS= 7
        this.BREAKPOINTS = 5 
    return this end
end

CONST()