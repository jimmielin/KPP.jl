using DifferentialEquations
using DiffEqBiological

include("Headers/state.jl")
include("Headers/registry.jl")
include("Core/kpp.jl")

"""
    jlkpp_Model()
this function runs the KPP.jl generated atmospheric chemistry model.
this automatic model has been generated by the following pre-processor:
    KPP.jl, version 1910.20 "August's Rhapsody"
    (c) 2019 Haipeng Lin <hplin@seas.harvard.edu>
    Licensed under the GNU General Public License, version 2 (excluding later versions)
"""
function jlkpp_Model()
    # $INITIAL_CONDITIONS$ #
    # $SETUP_TIMESTEPS$ #

    # $SPECIAL_DEFINES_BEFOREMODEL$ #

    # Initialize model environment
    spclist = jlkpp_registry
    state_model = ModelState(dt=dt,tstart=tstart,tend=tend,NX=1,NY=1,NZ=1,dthistory=0.0,lasthistory=0.0)

    # Useful shorthands
    IM, JM, LM  = state_model.NX, state_model.NY, state_model.NZ

    state_chm   = ChmState(spclist, IM, JM, LM)
    state_met   = MetState(IM, JM, LM)

    # "Compile" the mechanism for one run, so the ODE solver internals
    # are ready and type-specialized. See
    # https://stackoverflow.com/questions/47501844/julia-differentialequations-jl-speed
    @time jlkpp_Compile(jlkpp_mechanism, @view(state_chm.species[:,1,1,1]))

    # Write initial conditions to model (to be read by restart file)
    # For now only background concentrations are read
    println("Reading initial species concentrations from jlKPP defaults")
    for k in 1:LM
    for j in 1:JM
    for i in 1:IM
        jlkpp_Initialize_Defaults(@view state_chm.species[:,i,j,k])
    end
    end
    end

    # $SPECIAL_DEFINES_AFTERMODEL$ #

    println("* B e g i n   T i m e   S t e p p i n g !! *")
    @time @inbounds for t in state_model.tstart:state_model.dt:state_model.tend
        for k in 1:LM
        for j in 1:JM
        for i in 1:IM
            # $IN_GRID_LOOP_BEFORE$ #
            onestep = jlkpp_Timestep(rs=jlkpp_mechanism,
                                     u=@view(state_chm.species[:,i,j,k]),
                                     t=state_model.t, dt=state_model.dt,
                                     t_chm=state_met.TEMP[i,j,k])

            for n in 1:state_chm.nspecies
                state_chm.species[n,i,j,k] = onestep.u[end][n]
            end
            # $KLUDGE_FIX_FIXED_SPECIES$ #

            # $IN_GRID_LOOP_AFTER$ #
        end
        end
        end # end grid loop

        # this may be a good place for history output

        println("---> X-HRS: ", round(t/3600, digits=6))
        state_model.t += state_model.dt
    end

    # $SPECIAL_DEFINES_CLEANUP$ #

    # To actually get species concs you need to divide by CFACTOR

    println("******** SIMULATION ENDED ********")
end

@time jlkpp_Model()
