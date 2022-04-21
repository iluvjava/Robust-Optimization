using CSV
CSV_ALPHAS = CSV.File(open("data/alpha.csv"))
CSV_BETAS = CSV.File(open("data/beta.csv"))
CSV_P_GEN = CSV.File(open("data/ro_gen_data.csv"))
CSV_STORAGE = CSV.File(open("data/storage_data.csv"))
CSV_TRANS_LIMIT = CSV.File(open("data/transmission_limits.csv"))
CSV_POWERFLOW_BUS = CSV.File(open("data/sigma.csv"))

struct ConstParameters
    HORIZON 
    QS
    BUSES
    ARCS 
    BREAKPOINTS
    Φ
    function ConstParameters()
        this = new(
            3, # HORIZON
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
        this.generator_count = CSV_P_GEN|>length
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
        this.SR = CSV_P_GEN["SR"]
        this.RNSP = CSV_P_GEN["RNSP"]
        this.NSP = CSV_P_GEN["NSP"]
        this.alphas = ConvertCSV(CSV_ALPHAS)
        this.betas = ConvertCSV(CSV_BETAS)
    return this end
    
    function Base.length(this::Generators)
        return this.generator_count
    end
end

# !!! currently storage system models a single instance. 
mutable struct StorageSystem
    Efficiency::Number  # nu
    Distfactor::Number  # mu
    Capacity::Number
    CharingLim::Number
    DischargingLim::Number

    function StorageSystem()
        this = new()
        properties = CSV_STORAGE|>propertynames
        this.Efficiency = CSV_STORAGE[1][properties[2]]
        this.Distfactor = CSV_STORAGE[1][properties[3]]
        this.Capacity = CSV_STORAGE[1][properties[4]]
        this.CharingLim = CSV_STORAGE[1][properties[5]]
        this.DischargingLim = CSV_STORAGE[1][properties[6]]
    return this end
    
end


mutable struct Transmission
    Limit::Vector{Number}
    
    function Transmission()
        this = new()
        this.Limit = CSV_TRANS_LIMIT["Limit"]
    return this end
end


mutable struct Sigmas
    SigmaMatrix::Matrix{Number}
    function Sigmas()
        this = new()
        this.SigmaMatrix = ConvertCSV(CSV_POWERFLOW_BUS)
    return this end

end


function ConvertCSV(data::CSV.File)
    properties = data|>propertynames
    m = zeros(data|>length, (properties|>length) - 1)
    for II in 1:size(m, 1), III in 1:size(m, 2)
        m[II, III] = data[II][properties[III + 1]]
    end
return m end


CONST_PROBLEM_PARAMETERS = ConstParameters();
PRIMARY_GENERATORS = Generators(CSV_P_GEN);
STORAGE_SYSTEM = StorageSystem();
TRANSMISSION_SYSTEM = Transmission();
SIGMAS = Sigmas()
