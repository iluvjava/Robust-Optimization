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
            6, # QS: Quick start 
            6, # BUSES
            7, # ARCS
            5, # BREAKPOINTS
            1000000 # Budget
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
    return this end
    
end


CONST_PROBLEM_PARAMETERS = ConstParameters();
PROBLEM_PARAMETERS = Parameters()
PROBLEM_PARAMETERS.Tmind = 2
PROBLEM_PARAMETERS.Tminu = 2
