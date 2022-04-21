using CSV
ALPHAS = CSV.File(open("data/alpha.csv"))
BETAS = CSV.File(open("data/beta.csv"))
P_GEN = CSV.File(open("data/ro_gen_data.csv"))
STORAGE = CSV.File(open("data/storage_data.csv"))
TRANS_LIMIT = CSV.File(open("data/transmission_limits.csv"))


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

mutable struct Generators
    generator_count::Int
    Pmin::Vector{Number}  
    Pmax::Vector{Number}
    SU::Vector{Number}
    SD::Vector{Number}
    Tminu::Vector{Number}
    Tmind::Vector{Number}
    RD_bar::Vector{Number}
    RD::Vector{Number}
    RU_bar::Vector{Number}
    RU::Vector{Number}
    RREGU::Vector{Number}
    RREGD::Vector{Number}
    REGU::Vector{Number}
    REGD::Vector{Number}
    SR::Vector{Number}
    RNSP::Vector{Number}
    NSP::Vector{Number}

    alphas::Matrix{Number}
    betas::Matrix{Number}
    function Generators(file::CSV.File)
        this = new()
        this.generator_count = P_GEN|>length
        this.Pmin = file["Pmin"]
        this.Pmax = file["Pmax"]
        this.SU = file["SU"]
        this.SD = file["SD"]
        this.Tminu = file["Tminu"]
        this.Tmind = file["Tmind"]
        this.RD_bar = file["RD_bar"]
        this.RU_bar = file["RU_bar"]
        this.RU = file["RU"]
        this.RREGU = file["RREGU"]
        this.RREGD = file["RREGD"]
        this.REGU = file["REGU"]
        this.REGD = file["REGD"]
        this.SR = P_GEN["SR"]
        this.RNSP = P_GEN["RNSP"]
        this.NSP = P_GEN["NSP"]
        this.alphas = ConvertCSV(ALPHAS)
        this.betas = ConvertCSV(BETAS)
    return this end
    
end

mutable struct StorageSystem


end


mutable struct Transmission
    
end

function ConvertCSV(data::CSV.File)
    properties = data|>propertynames
    m = zeros(data|>length, (properties|>length) - 1)
    for II in 1:size(m, 1), III in 1:size(m, 2)
        m[II, III] = data[II][properties[III + 1]]
    end
return m end


CONST_PROBLEM_PARAMETERS = ConstParameters();
PRIMARY_GENERATORS = Generators(P_GEN);
