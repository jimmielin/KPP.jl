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
    #   INITIAL_CONDITIONS***
    #   SETUP_TIMESTEPS
    #   SPECIAL_DEFINES_BEFOREMODEL***
    #   SPECIAL_DEFINES_AFTERMODEL (with state_* initialized)***
    #   IN_GRID_LOOP_BEFORE***
    #   IN_GRID_LOOP_AFTER***
    #   KLUDGE_FIX_FIXED_SPECIES
    #   SPECIAL_DEFINES_CLEANUP***
    # KPP LAYER: Core/kpp.jl
    #   INITIALIZE_KPP --> Core/mechanism.jl
    #   TIMESTEP_PARAM (in jlkpp_Timestep)
    #   TIMESTEP_AFTER_PARAM***
    # KPP LAYER: Headers/registry.jl
    #   REGISTRY_SPCLIST
    #   REGISTRY_DEFAULT_IC (this is passed as jlkpp_ic to INITIAL_CONDITIONS)
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

    initialize_kpp_spc, registry_spclist = gen_spc_declaration(spclist)
    initialize_kpp_eqn = readfile_eqn(spclist, modelparams)
    initialize_kpp = """
jlkpp_mechanism = @empty_reaction_network KppReactionNetwork
""" * initialize_kpp_spc * initialize_kpp_eqn * """
addodes!(jlkpp_mechanism)
"""

    readfile_def(spclist, spcfixID, spcIC, modelparams)
    registry_default_ic, kludge_fix_fixed_species = gen_ic(spclist, spcfixID, spcIC, modelparams)

    # Model parameters write into kpp.jl as necessary
    setup_timesteps *= string("# These timesteps are set here by KPP.jl for convenience. If you are building a model, you should be able to read these from a configuration file.\n")
    setup_timesteps *= string("tstart = ", modelparams["tstart"], "\n")
    setup_timesteps *= string("dt = ", modelparams["dt"], "\n")
    setup_timesteps *= string("tend = ", modelparams["tend"], "\n")

    # KPP parameters write into kpp.jl - jlkpp_Timestep (TIMESTEP_PARAM)
    # In the format of p::Tuple{Float64, ...} for parameters. This is returned by readfile_eqn
    # FIXME: Hardcoded for now; need to coordinate with readfile_eqn's addparam!s
    timestep_param = "# SUN, TEMP, CFACTOR\n"
    timestep_param *= string("p = (jlkpp_SUN(t), T_chm, ", modelparams["CFACTOR"], ")")

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
    file = replace(file, "# \$KLUDGE_FIX_FIXED_SPECIES\$ #" => kludge_fix_fixed_species)
    open("./models/$mechanism_safe/main.jl", "w") do io
        write(io, file)
    end

    open("./models/$mechanism_safe/Core/kpp.jl", "r") do io
        file = read(io, String)
    end
    file = replace(file, "# \$TIMESTEP_PARAM\$ #" => timestep_param)
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
