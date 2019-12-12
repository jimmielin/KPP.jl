using DifferentialEquations
# $SETUP_EXTRA_DECLARES$ #

include("Headers/registry.jl")
include("Core/kpp.jl")

"""
    jlkpp_Model()
this function runs the KPP.jl generated atmospheric chemistry model.
this is a simplified version of the model where data is self-contained and returned
at the end of simulation.
    KPP.jl, version 1912.08
    (c) 2019 Haipeng Lin <hplin@seas.harvard.edu>
    Licensed under the GNU General Public License, version 2 (excluding later versions)
"""
function jlkpp_Model()
    # $INITIAL_CONDITIONS$ #
    # $SETUP_TIMESTEPS$ #

    # tstart, tend, dt
    # $SPECIAL_DEFINES_BEFOREMODEL$ #

    # Initialize model environment
    spclist = jlkpp_registry

    # # of levels for 3-D array
    IM = 1
    JM = 1
    LM = 1

    # chemistry state information
    chem_nspecies = length(spclist)
    chem_idx = spclist
    # chem_species = zeros(Float64, chem_nspecies, IM, JM, LM) # ordering: N X Y Z
    chem_species = [zeros(Float64, chem_nspecies) for k=1:LM, j=1:JM, i=1:IM] # ordering (ZYX)(N)

    # "Compile" the mechanism for one run, so the ODE solver internals
    # are ready and type-specialized. See
    # https://stackoverflow.com/questions/47501844/julia-differentialequations-jl-speed
    #
    # This returns an ODE problem
    # $GENERATE_OPROB$ #

    # Write initial conditions to model (to be read by restart file)
    # For now only background concentrations are read
    println("Reading initial species concentrations from jlKPP defaults")
    for k in 1:LM
    for j in 1:JM
    for i in 1:IM
        jlkpp_Initialize_Defaults(chem_species[k,j,i])
    end
    end
    end

    # $SPECIAL_DEFINES_AFTERMODEL$ #

    # println("Begin time stepping!")
    @inbounds for t in tstart:dt:tend
        for k in 1:LM
        for j in 1:JM
        for i in 1:IM
            # $IN_GRID_LOOP_BEFORE$ #
            jlkpp_Timestep!(chem_species[k,j,i], jlkpp_oprob,
                            t, dt,
                            300.0 # temperature for Chemistry
                            )
            # $KLUDGE_FIX_FIXED_SPECIES$ #

            # $IN_GRID_LOOP_AFTER$ #
        end
        end
        end # end grid loop

        # this may be a good place for history output

        # println("---> X-HRS: ", round(t/3600, digits=6))
    end

    # $SPECIAL_DEFINES_CLEANUP$ #

    # To actually get species concs you need to divide by CFACTOR

    # println("Time-stepping complete.")
end

@time jlkpp_Model()
