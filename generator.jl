"""
    generatemodel(spclist, spcfixID, spcIC, modelparams)
generate the jlkpp model with specified information, based on the dummy model code in template/
creating the generated model in models/$mechanism_safe, e.g. models/saprc99/
"""
function generatemodel()
    # Allocate global state information
    # Note that read functions should mutate them instead of direct global access.
    # Even though performance is not a major concern in pre-processing, stick to conventions.
    spclist = String[] # Species list ["NO", "O2", "NO2" ...] (from spc)
    spcfixID = Int[] # IDs of species that are fixed (from spc)
    spcIC = Float64[] # Default initial condition (IC) for each species (from model def)
    modelparams = Dict{String, Any}(
        # Model properties
        "species" => nothing,
        "fixed_species" => nothing,
        "reactions" => nothing,

        # Model parameters read from .def
        "CFACTOR" => nothing,
        "compat_integrator_hint" => nothing, # hint for integrator from .def file (not respected in KPP.jl)

        # Model timestepping parameters [s]
        "tstart" => nothing,
        "dt" => nothing,
        "tend" => nothing,

        # Model environmental parameters
        "TEMP" => 300.0,
    ) # Dictionary holding model parameters for generation

    # Generated Julia code blurbs to be populated by each of the routines
    # Insertion code points (marked with *** are not implemented)
    # MODEL LAYER: main.jl
    #   SETUP_EXTRA_DECLARES
    #   INITIAL_CONDITIONS***
    #   SETUP_TIMESTEPS
    #   GENERATE_OPROB
    #   SPECIAL_DEFINES_BEFOREMODEL***
    #   SPECIAL_DEFINES_AFTERMODEL (with state_* initialized)***
    #   IN_GRID_LOOP_BEFORE***
    #   IN_GRID_LOOP_AFTER***
    #   KLUDGE_FIX_FIXED_SPECIES
    #   SPECIAL_DEFINES_CLEANUP***
    # KPP LAYER: Core/kpp.jl
    #   EXTRA_LIBS_KPP (setup_extra_declares_kpp)
    #   INITIALIZE_KPP --> Core/mechanism.jl
    #   KPP_COMPILE
    #   TIMESTEP_PARAM (in jlkpp_Timestep)
    #   TIMESTEP_AFTER_PARAM***
    #   TIMESTEP_AFTER_SOLVE
    # KPP LAYER: Headers/registry.jl
    #   REGISTRY_SPCLIST
    #   REGISTRY_DEFAULT_IC (this is passed as jlkpp_ic to INITIAL_CONDITIONS)
    setup_extra_declares = ""
    setup_extra_declares_kpp = ""
    setup_timesteps = ""

    initialize_kpp = ""
    initialize_kpp_spc = ""
    initialize_kpp_eqn = ""
    timestep_param = ""

    registry_spclist = ""
    registry_default_ic = ""

    kludge_fixed_species = ""

    # Generate species information
    readfile_spc(spclist, spcfixID)
    if length(spclist) < 1
        error("KPP.jl: there must be at least 1 species in the list. Please check .spc definition")
    end
    modelparams["species"] = length(spclist)
    modelparams["fixed_species"] = length(spcfixID)

    if driver == "DiffBioEq"
        # DiffBioEq.jl driver, generate a reaction_network
        initialize_kpp_spc, registry_spclist = gen_spc_declaration(spclist)
        initialize_kpp_eqn, prodloss_diag_kpp_eqn = readfile_eqn(spclist, spcfixID, modelparams)
        initialize_kpp = """
jlkpp_mechanism = @empty_reaction_network KppReactionNetwork
""" * initialize_kpp_spc * initialize_kpp_eqn * """
addodes!(jlkpp_mechanism)
        """
    elseif driver == "Native"
        # Native driver, generates the KPP mechanism and problem
        initialize_kpp_spc, registry_spclist = gen_spc_declaration(spclist)
        initialize_kpp_eqn, prodloss_diag_kpp_eqn = readfile_eqn(spclist, spcfixID, modelparams)

        # Adaptive reduction (Santillana et al., 2010)
        # Add extra code to zero out changes for the "fixed" slow species
        if incl_optimize_adapt == true
            initialize_kpp_eqn *= """
    # Adaptive reduction: slow species set rate of change to 0; concentrations
    # will be manually updated later.
    for spc in 1:length(u)
        jlkpp_adapt_slow[spc] && (du[spc] = 0)
    end
"""
        end

        # Build the final KPP block code
        initialize_kpp = """
function jlkpp_mechanism_t!(du::Array{Float64, 1}, u::Array{Float64, 1}, p, t)
    param = p(t)
""" * initialize_kpp_spc * initialize_kpp_eqn * """
    nothing
end

function jlkpp_mechanism_prodloss!(prod, loss, u, p, t)
    param = p(t)
""" * initialize_kpp_spc * prodloss_diag_kpp_eqn * """
end
        """
    end

    readfile_def(spclist, spcfixID, spcIC, modelparams)
    registry_default_ic, kludge_fix_fixed_species = gen_ic(spclist, spcfixID, spcIC, modelparams)

    # Model setup declarations depending on the driver you use.
    # For DiffBioEq its just an include.
    # For Native usually there are optimization options, etc. you can use
    setup_extra_declares = ""
    setup_extra_declares_kpp = ""
    generate_oprob = ""
    if driver == "DiffBioEq"
        setup_extra_declares *= "using DiffEqBiological\n"

        setup_extra_declares_kpp *= "using DiffEqBase: AbstractReactionNetwork\n"

        # generate_oprob *= "jlkpp_oprob = jlkpp_Compile(jlkpp_mechanism, @view(chem_species[:,1,1,1]))"
        # @view is incompatble with SUNDIALS solver
        generate_oprob *= "jlkpp_oprob = jlkpp_Compile(jlkpp_mechanism_t!, chem_species[1,1,1])"
    elseif driver == "Native"
        setup_extra_declares *= "# KPP.jl Native Driver: Use precomputed Jacobian?\n"
        setup_extra_declares *= "const jlkpp_prejac = true\n"

        generate_oprob *= "jlkpp_oprob = jlkpp_Compile(jlkpp_mechanism_t!, chem_species[1,1,1])"
    end

    setup_extra_declares *= "# KPP.jl Solver: By default we use SUNDIALS as it provides exceptional performance. You can change the solver in kpp.jl.\n"
    setup_extra_declares *= "using Sundials\n"

    setup_extra_declares *= "# KPP.jl: Time stepping method. 'remake' remakes ODE problem, 'interfaced' uses the DiffEq Integrator interface.\n"
    setup_extra_declares *= "const jlkpp_stepmethod = 'remake'\n"

    # Model parameters write into main.jl as necessary
    setup_timesteps *= string("# These timesteps are set here by KPP.jl for convenience. If you are building a model, you should be able to read these from a configuration file.\n")
    setup_timesteps *= string("    tstart = ", modelparams["tstart"], "\n")
    setup_timesteps *= string("    dt = ", modelparams["dt"], "\n")
    setup_timesteps *= string("    tend = ", modelparams["tend"], "\n")

    # KPP ODEProblem compilation setup.
    # This generates the ODEProblem and does one solve iteration
    if driver == "DiffBioEq"
    kpp_compile = """
    function jlkpp_Compile(rs::AbstractReactionNetwork, u_scratch::AbstractArray{Float64, 1})::ODEProblem
        # println("KPP.jl: Running ODE solver for type specialization")
        ODEProblem(rs, u_scratch, (0.0, 600.0), t -> (0, 300.0, 2.4476e13))
        # solve(oprob, alg_hints=[:stiff], dense=false, save_on=false, calck=false)
        # println("KPP.jl: Finished")
        # oprob # Return the ODEProblem
    end
    """
    elseif driver == "Native"
    kpp_compile = """
    function jlkpp_Compile(f, u_scratch::AbstractArray{Float64, 1})::ODEProblem
        ODEProblem(f, u_scratch, (0.0, 600.0), t -> (0, 300.0, 2.4476e13))
        # solve(oprob, alg_hints=[:stiff], dense=false, save_on=false, calck=false)
        # println("KPP.jl: Finished")
    end
    """
    end


    # KPP parameters write into kpp.jl - jlkpp_Timestep (TIMESTEP_PARAM)
    # In the format of p::Tuple{Float64, ...} for parameters. This is returned by readfile_eqn
    # FIXME: Hardcoded for now; need to coordinate with readfile_eqn's addparam!s
    timestep_param = "\n"
    timestep_param_extra_code = ""
    timestep_after_solve = "\n"
    kpp_code_extra = ""
    if incl_optimize_adapt == true && driver == "Native"
        # Adaptive reduction mechanism diagnose and update slow species
        # (Native driver only)
        timestep_param *= "    # Adaptive reduction algorithm (Santillana et al., 2010)\n"
        timestep_param *= "    # Updates global scratchpads jlkpp_scratch_* and jlkpp_adapt_slow\n"
        timestep_param *= "    _optim_slow_count = jlkpp_adapt_diagn(u, t)\n"
        # timestep_param_extra_code *= ", _optim_slow_spc"
        # now uses scratchpads, no longer through parameters. I hate global variables
        # as much as you do, but this saves allocations
        timestep_after_solve *= "    # Adaptive reduction: update slow spc\n"
        timestep_after_solve *= "    jlkpp_adapt_slow!(u, dt)\n"

        # Add scratchpads for jlkpp_adapt to work with and save production/losses
        # and current slow species state
        setup_extra_declares_kpp *= "const jlkpp_scratch_p = zeros(Float64, " * string(length(spclist)) * ")\n"
        setup_extra_declares_kpp *= "const jlkpp_scratch_l = zeros(Float64, " * string(length(spclist)) * ")\n"
        setup_extra_declares_kpp *= "const jlkpp_adapt_slow = zeros(Bool, " * string(length(spclist)) * ")\n"

        setup_extra_declares_kpp *= "const jlkpp_CFACTOR = " * string(modelparams["CFACTOR"])

        # Inject the adaptive reduction mechanism code implementing
        # jlkpp_adapt_diagn and jlkpp_adapt_slow!
        open("./optimizers/santillana2010.jl", "r") do io
            adapt_diagn_code = read(io, String)
            kpp_code_extra *= adapt_diagn_code
        end
    end

    timestep_param *= "    # SUN, TEMP, CFACTOR\n"
    timestep_param *= string("    p = t -> (jlkpp_SUN(t), T_chm, ", modelparams["CFACTOR"], timestep_param_extra_code, ")\n")

    # ----------------------- END PREP ------------------------- #

    # debugging code for emergency escape
    # return

    # Copy into a template folder. Prompt if destination is existing
    if isdir("models/$mechanism_safe")
        error("KPP.jl: Destination model folder models/$mechanism_safe exists. Check if you want to overwrite and delete as needed.")
    end

    # Act on IO!
    println("Creating target model ./models/$mechanism_safe")
    cp("./template/", "./models/$mechanism_safe")

    if !isfile("./models/$mechanism_safe/main.jl")
        error("KPP.jl: Could not create destination folder models/$mechanism_safe. Check permissions.")
    end

    # Modify files as necessary
    open("./models/$mechanism_safe/main.jl", "r") do io
        file = read(io, String)
    end
    file = replace(file, "# \$SETUP_TIMESTEPS\$ #" => setup_timesteps)
    file = replace(file, "# \$SETUP_EXTRA_DECLARES\$ #" => setup_extra_declares)
    file = replace(file, "# \$GENERATE_OPROB\$ #" => generate_oprob)
    file = replace(file, "# \$KLUDGE_FIX_FIXED_SPECIES\$ #" => kludge_fix_fixed_species)
    open("./models/$mechanism_safe/main.jl", "w") do io
        write(io, file)
    end

    open("./models/$mechanism_safe/Core/kpp.jl", "r") do io
        file = read(io, String)
    end
    file = replace(file, "# \$EXTRA_LIBS_KPP\$ #" => setup_extra_declares_kpp)
    file = replace(file, "# \$TIMESTEP_PARAM\$ #" => timestep_param)
    file = replace(file, "# \$TIMESTEP_AFTER_SOLVE\$ #" => timestep_after_solve)
    file = replace(file, "# \$KPP_COMPILE\$ #" => kpp_compile)
    file = replace(file, "# \$KPP_CODE_EXTRA\$ #" => kpp_code_extra)
    file = replace(file, "# \$MECHANISM\$ #" => initialize_kpp)
    open("./models/$mechanism_safe/Core/kpp.jl", "w") do io
        write(io, file)
    end

    open("./models/$mechanism_safe/Headers/registry.jl", "r") do io
        file = read(io, String)
    end
    file = replace(file, "# \$REGISTRY_SPCLIST\$ #" => registry_spclist)
    file = replace(file, "# \$REGISTRY_DEFAULT_IC\$ #" => registry_default_ic)
    open("./models/$mechanism_safe/Headers/registry.jl", "w") do io
        write(io, file)
    end

    println("Finished generating model $mechanism_safe")
end
