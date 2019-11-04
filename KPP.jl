# Configure here:
const mechanism = "small_strato"
const mechanism_safe = mechanism

# NO USER CONFIGURABLE PARAMETERS BELOW
include("parser.jl")
include("generator.jl")

# Change to current directory
cd(dirname(@__FILE__))

# Code below
if !isfile("./inputs/$mechanism.spc") || !isfile("./inputs/$mechanism.eqn") || !isfile("./inputs/$mechanism.def")
    error("Cannot find mechanism files: verify if spc, eqn and def exist for $mechanism")
end

# Notes:
# (1) Regexes are written with assumption of line by line parsing
# (2) Rate law functions ARR, FALL, EP3... are hard coded with shims for now

generatemodel()
