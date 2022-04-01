struct ConstParameters
    HORIZON 
    GENERATOR_PRIMARY
    QS
    BUSES
    ARCS 
    BREAKPOINTS
    Φ
    function ConstParameters()
        this = new(
            3, # HORIZON
            6, # GENERATOR_PRIMARY
            6, 
            6,
            7,
            5,
            1000000
        )
    return this end
end

mutable struct Parameters
    Pmin; Pmax;
    REGU; REGD
    RMT
    SU
    SD
    Tminu; Tmind
    RU; RD
    RU_bar; RD_bar
    NSP
    H
    S
    Pmin_prime; Pmax_prime
    REGU_prime; REGD_prime
    RMT_prime; SU_prime
    SD_prime
    Tminu_prime; Tmind_prime
    RU_prime; RD_prime
    RU_bar_prime; RD_bar_prime; NSP_prime
    H_prime
    S_prime
    req_data
    RREGU; RREGD
    RNSP
    function Parameters()
        this = new()
        this.Tminu = 1
        this.Tmind = 1
    return this end
    
end

if !(@isdefined(CONST_PROBLEM_PARAMETERS))
    const CONST_PROBLEM_PARAMETERS = ConstParameters();
    const PROBLEM_PARAMETERS = Parameters()
end