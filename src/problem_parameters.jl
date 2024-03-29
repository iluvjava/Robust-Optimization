using CSV

"Primary generators fuel costs parameters. "
const CSV_P_ALPHAS = CSV.File(open("$DATA_DIR/alpha.csv"))

"Primary generators fuel costs. "
const CSV_P_BETAS = CSV.File(open("$DATA_DIR/beta.csv"))

"Primary generator generation limit. "
const CSV_P_GEN = CSV.File(open("$DATA_DIR/ro_gen_data.csv"))

"The data for the storage devices. "
const CSV_STORAGE = CSV.File(open("$DATA_DIR/storage_data.csv"))

const CSV_DEMAND_RESPONSE = CSV.File(open("$DATA_DIR/dr_data.csv"))


mutable struct ConstParameters
    HORIZON::Int
    Φ::Number

    function ConstParameters()
        this = new()
        this.HORIZON = 24; 
        this.Φ = 1e8
    return this end
end

"""
- The alpha file, 
- beta file. 
- generator file. 
"""
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
    initial_status::Vector{Number}
    initial_pg::Vector{Number}
    
    function Generators(
        gen_file::CSV.File, 
        alpha_file::CSV.File, 
        beta_file::CSV.File
    )
        this = new()
        # [x]: Change here, adapt to the new data format. 
        this.generator_count = gen_file|>length
        this.Pmin = gen_file["Pmin"]
        this.Pmax = gen_file["Pmax"]
        this.SU = gen_file["SU"]
        this.SD = gen_file["SD"]
        this.Tminu = gen_file["Tminu"]
        this.Tmind = gen_file["Tmind"]
        this.RD_bar = gen_file["RD_bar"]
        this.RD = gen_file["RD"]
        this.RU_bar = gen_file["RU_bar"]
        this.RU = gen_file["RU"]
        # this.RREGU = gen_file["RREGU"]
        # this.RREGD = gen_file["RREGD"]
        # this.REGU = gen_file["REGU"]
        # this.REGD = gen_file["REGD"]
        # this.SR = CSV_P_GEN["SR"]
        # this.RNSP = CSV_P_GEN["RNSP"]
        # this.NSP = CSV_P_GEN["NSP"]
        this.alphas = ConvertCSV(alpha_file)
        this.betas = ConvertCSV(beta_file)
        this.initial_status = gen_file["status"]
        this.initial_pg = gen_file["Pg"]

    return this end
    
    function Base.length(this::Generators)
        return this.generator_count
    end
end


mutable struct StorageSystem
    "ν, greek nu"
    Efficiency::Array{Number}  # nu 
    # Distfactor::Array{Number}  # mu
    "H̄"
    Capacity::Array{Number}    # H̅
    "Ḡ⁺"
    CharingLim::Array{Number}  # G̅+
    "Ḡ⁻"
    DischargingLim::Array{Number}  # G̅-
    "total number of storage system"
    s  

    function StorageSystem()
        this = new()
        properties = CSV_STORAGE|>propertynames
        this.Efficiency = CSV_STORAGE[properties[2]]
        this.Capacity = CSV_STORAGE[properties[3]]
        this.CharingLim = CSV_STORAGE[properties[4]]    
        this.DischargingLim = CSV_STORAGE[properties[5]]
        this.s = CSV_STORAGE |> length
    return this end

    function Base.length(this::StorageSystem)
    return this.s end
    
end


mutable struct DemandResponse
    level::Int
    rho::Vector
    R::Dict{Int, Number}

    function DemandResponse(file::CSV.File)
        this = new()
        this.level = length(file)
        this.rho = file["rho"]
        this.R = zip(1:length(file["R"]), file["R"])|>collect|>Dict
        this.R[0] = 0
        return this
    end

end


function ConvertCSV(data::CSV.File)
    properties = data|>propertynames
    m = zeros(data|>length, (properties|>length) - 1)
    for II in 1:size(m, 1), III in 1:size(m, 2)
        m[II, III] = data[II][properties[III + 1]]
    end
return m end


"The constant parameters for the problem. "
CONST_PROBLEM_PARAMETERS = ConstParameters()

"All the data related to the primary generators. "
PRIMARY_GENERATORS = Generators(CSV_P_GEN, CSV_P_ALPHAS, CSV_P_BETAS)

"All the data related to the storage system. "
STORAGE_SYSTEM = StorageSystem()


DEMAND_RESPONSE = DemandResponse(CSV_DEMAND_RESPONSE)


# "Stores all the structures of the busses in different transmission line. "
# const BUSES = Buses()

@info "Problem Parameters Successfully Loaded"

