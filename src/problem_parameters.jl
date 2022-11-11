using CSV
# Primary Generators
CSV_P_ALPHAS = CSV.File(open("data/alpha.csv"))
CSV_P_BETAS = CSV.File(open("data/beta.csv"))
CSV_P_GEN = CSV.File(open("data/ro_gen_data.csv"))
CSV_S_GEN = CSV.File("data/quick_start_data.csv")
CSV_S_ALPHAS = CSV.File("data/alpha_prime.csv")
CSV_S_BETAS = CSV.File("data/beta_prime.csv")
CSV_DISFACTOR = CSV.File("data/mu.csv")

# Others
CSV_STORAGE = CSV.File(open("data/storage_data.csv"))
CSV_TRANS_LIMIT = CSV.File(open("data/transmission_limits.csv"))
CSV_POWERFLOW_BUS = CSV.File(open("data/sigma.csv"))

mutable struct ConstParameters
    HORIZON 
    Φ
    RREGD::Vector{Number}
    RREGU::Vector{Number}
    RNSP::Vector{Number}

    function ConstParameters()
        this = new()
        this.HORIZON = 4; 
        this.Φ = 1e8; 
        this.RREGD = zeros(this.HORIZON)
        this.RREGU = zeros(this.HORIZON)
        this.RNSP = zeros(this.HORIZON)
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
    RREGU::Vector{Number}  # wrong parameters 
    RREGD::Vector{Number}  # wrong parameters 
    REGU::Vector{Number}
    REGD::Vector{Number}
    SR::Vector{Number}
    RNSP::Vector{Number}   # move to global
    NSP::Vector{Number}

    alphas::Matrix{Number}
    betas::Matrix{Number}
    function Generators(
        gen_file::CSV.File, 
        alpha_file::CSV.File, 
        beta_file::CSV.File
    )
        this = new()
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
        this.RREGU = gen_file["RREGU"]
        this.RREGD = gen_file["RREGD"]
        this.REGU = gen_file["REGU"]
        this.REGD = gen_file["REGD"]
        this.SR = CSV_P_GEN["SR"]
        this.RNSP = CSV_P_GEN["RNSP"]
        this.NSP = CSV_P_GEN["NSP"]
        this.alphas = ConvertCSV(alpha_file)
        this.betas = ConvertCSV(beta_file)
    return this end
    
    function Base.length(this::Generators)
        return this.generator_count
    end
end


mutable struct StorageSystem
    Efficiency::Array{Number}  # nu 
    # Distfactor::Array{Number}  # mu
    Capacity::Array{Number}    # H̅
    CharingLim::Array{Number}  # G̅+
    DischargingLim::Array{Number}  # G̅-
    
    s  # total number of storage system 
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


mutable struct Transmission
    Limit::Vector{Number}
    
    function Transmission()
        this = new()
        this.Limit = CSV_TRANS_LIMIT["Limit"]
    return this end
end


"""
    Contains info for both buses and transmission line interactions!!!
        [Busses Index, Transimission Line]
"""
mutable struct DataMatrix
    the_matrix::Matrix{Number}
    
    function DataMatrix(input_file::CSV.File)
        this = new()
        this.the_matrix = ConvertCSV(input_file)'
    return this end
    
    function Base.size(this::DataMatrix)
    return size(this.the_matrix) end 
end


function ConvertCSV(data::CSV.File)
    properties = data|>propertynames
    m = zeros(data|>length, (properties|>length) - 1)
    for II in 1:size(m, 1), III in 1:size(m, 2)
        m[II, III] = data[II][properties[III + 1]]
    end
return m end


"""
Models a group of generators that are in the same busses. 
It has primary/secondary generator in it, and it maps the buses to a set of 
indices for generators in that bus. 

* The array could be empty. 
* CURRENTLY NO DATA IS PREPARE FOR THIS YET, SO IT'S HARD CODED IN. 
"""
mutable struct Buses
    primary::Dict{Int, Set{Int}}
    secondary::Dict{Int, Set{Int}}
    
    """
        Generator info are HARD CODED IN!!! 
    """
    function Buses(num_of_busses=6, num_primary_gen=6, num_secondary_gen=1)
        this = new()
        this.primary = Dict()
        this.secondary = Dict()
        for II in 1:num_of_busses
            this.primary[II] = Set{Int}()
            push!(this.primary[II], II)
            this.secondary[II] = Set{Int}()
        end
        push!(this.secondary[2], 1)
    return this end
end

"The constant parameters for the problem. "
CONST_PROBLEM_PARAMETERS = ConstParameters()
"All the data related to the primary generators. "
PRIMARY_GENERATORS = Generators(CSV_P_GEN, CSV_P_ALPHAS, CSV_P_BETAS)
"All the data related to the secondary generators. "
SECONDARY_GENERATORS = Generators(CSV_S_GEN, CSV_S_ALPHAS, CSV_S_BETAS)
"All the data related to the storage system. "
STORAGE_SYSTEM = StorageSystem()
"All the data related to the transmission system for the model. "
TRANSMISSION_SYSTEM = Transmission()
"The sigma transmission coefficients between the powerline and the generators. "
SIGMAS = DataMatrix(CSV_POWERFLOW_BUS)
"The disfactor parameters. "
DISFACTORS = DataMatrix(CSV_DISFACTOR)
"Stores all the structures of the busses in different transmission line. "
BUSES = Buses()

@info "Problem Parameters Successfully Loaded"