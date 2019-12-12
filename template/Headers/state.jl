mutable struct ModelState
    # current model time [seconds since 0000z model start day]
    t::Float64

    # time stepping
    dt::Float64
    tstart::Float64
    tend::Float64

    # grid dimensions
    NX::Int
    NY::Int
    NZ::Int

    # time for history output (set to 0.0 if not desired)
    dthistory::Float64
    lasthistory::Float64

    # Constructor function boilerplate
    # This should no longer be necessary in the future https://github.com/JuliaLang/julia/issues/10146
    function ModelState(;dt=900.0,tstart=0.0,tend=3600.0,NX=1,NY=1,NZ=1,dthistory=0.0,lasthistory=0.0)
        new(0.0,dt,tstart,tend,NX,NY,NZ,dthistory,lasthistory)
    end
end

mutable struct ChmState
    # number of stored chemical species
    nspecies::Int

    # species dictionary map from String -> N
    # this should be a O(1) lookup
    idx::Dict{String, Int}

    # species concentrations array
    # (N, Z, Y, X) for a generic 3-D eulerian grid
    #
    # However, Sundials solver is a Fortran interface and does not support @view
    # which leads to the decision of a Array{Array{Float64, 1}, 3} structure
    # which is essentially (Z, Y, X)(N), loop over ZYX to get N
    # species::Array{Float64, 4}
    species::Array{Array{Float64, 1}, 3}

    # Constructor function
    function ChmState(spclist::Dict{String, Int},NX=1,NY=1,NZ=1)
        nspecies = length(spclist)
        # species = zeros(Float64,nspecies,NX,NY,NZ) # 4-D species array or
        species = [zeros(Float64, nspecies) for k=1:NZ, j=1:NY, i=1:NX] # memorder ZYX
        new(nspecies,spclist,species)
    end
end

mutable struct MetState
    # SUN       Sunlight intensity for chemistry [1]
    SUN::Array{Float64, 3}

    # TEMP      Temperature for chemistry [K]
    TEMP::Array{Float64, 3}

    # Constructor function
    function MetState(NX=1,NY=1,NZ=1)
        new(zeros(Float64,NX,NY,NZ),fill(300.0,(NX,NY,NZ)))
    end
end
