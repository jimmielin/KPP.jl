function readfile_spc(spclist::Array{String, 1}, spcfixID::Array{Int, 1})
    println("Begin parsing $mechanism.spc...")
    parserstate = "none" # none, defvar, deffix,

    open("./inputs/$mechanism.spc") do f
        for (i, line) in enumerate(eachline(f))
            # #DEFVAR, #DEFFIX, #include ...
            # CH3OOH = C + 4H + 2O ; {CH4O2      methylperoxy alcohol}
            # O3  	= 3O;
            # ACET =  IGNORE;

            # Match commands.
            if (m = match(r"\s*#([a-zA-Z]+)\s*(.*)", line)) != nothing
                cmd = lowercase(m[1])
                if cmd == "include"
                        println("# NotImplemented: include command in .eqn file.")
                elseif cmd == "deffix"
                    parserstate = "deffix"
                elseif cmd == "defvar"
                    parserstate = "defvar"
                else
                    println("# Unsupported command ", cmd, " parsed. Skipping. You might want to report this.")
                end
            # Match species. Warning: This regex DOES NOT MATCH atoms for now.
            # In the future implementing atom & IGNORE check (2H + C + O) or (IGNORE); {}
            # will be necessary.
            elseif (m = match(r"^\s*([a-zA-Z0-9_]+)\s*=\s*", line)) != nothing
                spc = m[1] # species are case-sensitive
                if spc in spclist
                    print("# ERROR: Duplicate species $spc in list found. Skipping.")
                end
                push!(spclist, spc)
                if parserstate == "deffix"
                    push!(spcfixID, length(spclist))
                end
            end

            m = nothing
        end
    end

    println("Finished parsing species file for a total of ", length(spclist), " species")
end

"""
    gen_spc_declaration(spclist)
generates the specs for species declaration, both for the KPP layer and the registry
specification within Headers/registry.jl
"""
function gen_spc_declaration(spclist::Array{String, 1})
    registry_spc = "const jlkpp_registry = Dict("
    registry_meta_spc = "const jlkpp_registry_meta = Dict{String, Any}["
    initialize_kpp_spc = ""

    # FIXME: There really should be escaping for the name... right now it does not match
    # anything other than [a-zA-Z0-9_], but we should never trust user input (hplin, 10/20/19)
    for (index, name) in enumerate(spclist)
        # registry_spc expects a Dictionary of type Dict{String, Int}
        registry_spc *= "\"$name\" => $index, \n"

        registry_meta_spc *= "Dict(\"index\" => $index, \"name\" => \"$name\"), \n"

        # initialize_kpp expects a species list for DiffBioEq driver
        if driver == "DiffBioEq"
            initialize_kpp_spc *= "addspecies!(jlkpp_mechanism, :$name)\n"
        elseif driver == "Native"
            # Instead, map species to a list like so
            # O, O1D, O3, NO, NO2, M, O2 = u  # species
            initialize_kpp_spc *= ", $name"
        end
    end

    # Chop and process for Native driver
    if driver == "Native"
        initialize_kpp_spc *= " = u\r\n"
        initialize_kpp_spc  = initialize_kpp_spc[3:end] # chop off the initial ", " (2c)
    end

    registry_spc *= ")\n"
    registry_meta_spc *= "]"
    registry_spc *= registry_meta_spc

    initialize_kpp_spc, registry_spc
end

function parsekcoeff(str::AbstractString)
    # Hard-coded reaction coefficient functions, shimmed for compatibility with
    # ancient F90 code. Replace with tokens (:A0, :B0, :C0). Enclose tokens with ()
    # for safety
    # Brought to you by hplin 10/18/19

    # Arrhenius equation
    # Params: TEMP
    rgx_ARR = r"ARR\(\s*(?<A0>[\d\.e\-\+]+)\s*,\s*(?<B0>[\d\.e\-\+]+)\s*,\s*(?<C0>[\d\.e\-\+]+)\s*\)"
    tpl_ARR = s"jlkpp_ARR(TEMP, \g<A0>, \g<B0>, \g<C0>)"

    # Simplified Arrhenius, with two arguments
    # Params: TEMP
    rgx_ARR2 = r"ARR2\(\s*(?<A0>[\d\.e\-\+]+)\s*,\s*(?<B0>[\d\.e\-\+]+)\)"
    tpl_ARR2 = s"jlkpp_ARR2(TEMP, \g<A0>, \g<B0>)"

    # EP2
    # K0 = (:A0) * exp(-(:C0)/TEMP)
    # K2 = (:A2) * exp(-(:C2)/TEMP)
    # K3 = (:A3) * exp(-(:C3)/TEMP) * 1.0e6 * CFACTOR
    # Params: TEMP, CFACTOR
    rgx_EP2 = r"EP2\(\s*(?<A0>[\d\.e\-\+]+)\s*,\s*(?<C0>[\d\.e\-\+]+)\s*,\s*(?<A2>[\d\.e\-\+]+)\s*,\s*(?<C2>[\d\.e\-\+]+)\s*,\s*(?<A3>[\d\.e\-\+]+)\s*,\s*(?<C3>[\d\.e\-\+]+)\s*\)"
    tpl_EP2 = s"jlkpp_EP2(TEMP, CFACTOR, \g<A0>, \g<C0>, \g<A2>, \g<C2>, \g<A3>, \g<C3>)"

    # EP3
    # K1 = (:A1) * exp(-(:C1)/TEMP)
    # K2 = (:A2) * exp(-(:C2)/TEMP)
    rgx_EP3 = r"EP3\(\s*(?<A1>[\d\.e\-\+]+)\s*,\s*(?<C1>[\d\.e\-\+]+)\s*,\s*(?<A2>[\d\.e\-\+]+)\s*,\s*(?<C2>[\d\.e\-\+]+)\s*\)"
    tpl_EP3 = s"jlkpp_EP3(TEMP, CFACTOR, \g<A1>, \g<C1>, \g<A2>, \g<C2>)"

    # FALL
    # K0 = (:A0) * exp(-(:B0)/TEMP) * (TEMP/300.0)^(:C0) * CFACTOR * 1.0e6
    # K1 = (:A1) * exp(-(:B1)/TEMP) * (TEMP/300.0)^(:C1)
    # K1' = ((:A0) * exp(-(:B0)/TEMP) * (TEMP/300.0)^(:C0) * CFACTOR * 1.0e6)/((:A1) * exp(-(:B1)/TEMP) * (TEMP/300.0)^(:C1))
    rgx_FALL = r"FALL\(\s*(?<A0>[\d\.e\-\+]+)\s*,\s*(?<B0>[\d\.e\-\+]+)\s*,\s*(?<C0>[\d\.e\-\+]+)\s*,\s*(?<A1>[\d\.e\-\+]+)\s*,\s*(?<B1>[\d\.e\-\+]+)\s*,\s*(?<C1>[\d\.e\-\+]+)\s*,\s*(?<CF>[\d\.e\-\+]+)\s*\)"
    tpl_FALL = s"jlkpp_FALL(TEMP, CFACTOR, \g<A0>, \g<B0>, \g<C0>, \g<A1>, \g<B1>, \g<C1>, \g<CF>)"

    if driver == "DiffBioEq"
        tpl_ARR = s"((\g<A0>) * exp(-(\g<B0>)/TEMP) * (TEMP/300.0)^(\g<C0>))"
        tpl_ARR2 = s"((\g<A0>) * exp((\g<B0>)/TEMP))"
        tpl_EP2 = s"((\g<A0>) * exp(-(\g<C0>)/TEMP) + ((\g<A3>) * exp(-(\g<C3>)/TEMP) * 1.0e6 * CFACTOR)/(1.0 + ((\g<A3>) * exp(-(\g<C3>)/TEMP) * 1.0e6 * CFACTOR)/((\g<A2>) * exp(-(\g<C2>)/TEMP))))"
        tpl_EP3 = s"((\g<A1>) * exp(-(\g<C1>)/TEMP) + ((\g<A2>) * exp(-(\g<C2>)/TEMP)) * 1.0e6 * CFACTOR)"
        tpl_FALL = s"((((\g<A0>) * exp(-(\g<B0>)/TEMP) * (TEMP/300.0)^(\g<C0>) * CFACTOR * 1.0e6)/(1.0 + ((\g<A0>) * exp(-(\g<B0>)/TEMP) * (TEMP/300.0)^(\g<C0>) * CFACTOR * 1.0e6)/((\g<A1>) * exp(-(\g<B1>)/TEMP) * (TEMP/300.0)^(\g<C1>)))) * (\g<CF>)^(1.0/(1.0 + log10((\g<A0>) * exp(-(\g<B0>)/TEMP) * (TEMP/300.0)^(\g<C0>) * CFACTOR * 1.0e6)/((\g<A1>) * exp(-(\g<B1>)/TEMP) * (TEMP/300.0)^(\g<C1>)))^2))"
        # FIXME HACK: DiffBioEq has problems interpolating functions so we have to use the old way...
    end

    # replace away, then tokenize
    str = replace(str, rgx_ARR => tpl_ARR)
    str = replace(str, rgx_ARR2 => tpl_ARR2)
    str = replace(str, rgx_EP2 => tpl_EP2)
    str = replace(str, rgx_EP3 => tpl_EP3)
    str = replace(str, rgx_FALL => tpl_FALL)

    if driver == "DiffBioEq"
        string(":(", str, ")")
    elseif driver == "Native"
        string("(", str, ")")
    end
end

function readfile_eqn(spclist::Array{String, 1}, spcfixID::Array{Int, 1}, modelparams::Dict{String, Any})
    println("Begin parsing $mechanism.eqn...")

    nspc = length(spclist)

    codeeqn = ""
    prodloss_diag_kpp_eqn = "" # P/L diagnostics

    # TODO: Hardcoded parameters - fix later (TEMP, SUN, CFACTOR)
    if driver == "DiffBioEq"
        codeeqn *= "addparam!(jlkpp_mechanism, :SUN)\r\n"
        codeeqn *= "addparam!(jlkpp_mechanism, :TEMP)\r\n"
        codeeqn *= "addparam!(jlkpp_mechanism, :CFACTOR)\r\n"
    elseif driver == "Native"
        codeeqn *= "SUN, TEMP, CFACTOR = param\r\n"
        prodloss_diag_kpp_eqn *= "SUN, TEMP, CFACTOR = param\r\n"
    end

    rxn_count = 0

    # If we are in the "Native" driver, we need to keep track of the prod, loss and rates manually
    # This is through parsing equations, reading stoiochrometric coeffs, and stuff...
    # We store string representations of: net change, prod and loss, by species index,
    # in three arrays organized by SpcID => String
    if driver == "Native"
        netrates = ["du[$i] = 0" for i in 1:nspc]
        prod     = ["prod[$i] = 0" for i in 1:nspc]
        loss     = ["loss[$i] = 0" for i in 1:nspc]
        # append netrates, prod, loss by space, sign, ratecoef, spc, myspc
        @show typeof(netrates)
    end

    # Since there may be line breaks in between equations this file is opened differently
    # and not read by line (.eqn)

    # {\d+}\s*(?<reaction>[a-zA-Z0\d][a-zA-Z0\d\s\+\=\.\_]*)\s*:\s*(?<rate>[^;]*);   version 1, hplin, 10/18/19
    rgx_equation = r"[<{][rR\d]+[}>]\s*(?<reaction>[a-zA-Z0\d][a-zA-Z0\d\s\+\=\.\_]*)\s*:\s*(?<rate>[^;]*);"
    rgx_coef     = r"(?<coef>[\d\.]+)?\s*(?<chm>[a-zA-Z\_\d]+)"

    eqnio  = open("./inputs/$mechanism.eqn", "r")
    eqnstr = read(eqnio, String)

    for m in eachmatch(rgx_equation, eqnstr)
        rxn_count += 1
        # m.reaction, m.rate
        # for log purposes, verify if is photolytic reaction
        is_photol = match(r"\+\s*hv", m[:reaction]) != nothing
        rxn = replace(strip(m[:reaction]), r"\+\s*hv" => "") # photol reactions are handled through rate constants
        rxn = replace(rxn, r"[\t\n\r]" => "") # remove outstanding whitespace
        rxn = replace(rxn, r"  +" => " ") # remove more than 2 spaces
        rxn = replace(rxn, "=" => "-->")
        if driver == "DiffBioEq"
            rxn = ":($rxn)" # enclose in statement
        end

        # parse the rate coefficients
        coef = parsekcoeff(m[:rate])

        # write the code...
        # println("Eqn#$rxn_count : $rxn @ $coef")
        # addreaction!(small_strato, :(2.643e-10*SUN*SUN*SUN), :(O2 --> 2O)) # photol
        if driver == "DiffBioEq"
            codeeqn *= "addreaction!(jlkpp_mechanism, $coef, $rxn)\r\n"
        elseif driver == "Native"
            # Further diagnosis of the system is necessary...
            # append netrates, prod, loss by space, sign, ratecoef, spc, myspc

            # Parse this equation further...
            # rxn: doesn't get more complicated than
            # MEK + OH --> 0.37RO2_R + 0.042RO2_N + 0.616R2O2+ 0.492CCO_O2 + 0.096RCO_O2 + 0.115HCHO + 0.482CCHO + 0.37RCHO
            reactants, products = strip.(split(rxn, "-->"))
            reactants = strip.(split(reactants, "+"))
            products = strip.(split(products, "+"))
            rate_expr = "(" * coef * ")"
            for x in reactants
                # Build the rate. match using
                # (?<coef>[\d\.]+)?\s*(?<chm>[a-zA-Z\_\d]+)        version 1, hplin, 11/9/2019
                rm = match(rgx_coef, x)
                rcoef = (rm[:coef] == nothing ? 1.0 : parse(Float64, rm[:coef]))
                # Seriously, this is a very stupid construct of indexof or anything for that matter.
                # CartesianIndex? Why reinvent the Tuple?
                ridx = findfirst(isequal(rm[:chm]), spclist)
                if ridx == nothing
                    error("Fatal parser error at equation ", rxn, " -- could not find species ", rm[:chm])
                end
                # ridx = ridx[2]
                rspc = rm[:chm]

                # Check passed, can build rate expression
                if rcoef == 1
                    rate_expr *= "*" * rspc
                else
                    rate_expr *= "*(" * rspc * "^$rcoef)"
                end
            end

            # Now with the complete rate expression build the loss and prod...
            for x in reactants
                # Build the loss and rxn rate.
                rm = match(rgx_coef, x)
                rcoef = (rm[:coef] == nothing ? 1 : parse(Float64, rm[:coef]))
                ridx = findfirst(isequal(rm[:chm]), spclist)
                if ridx == nothing
                    error("Fatal parser error at equation ", rxn, " -- could not find species ", rm[:chm])
                end
                # ridx = ridx[2]
                rspc = rm[:chm]

                # Check if fixed species... in spcfixID
                if findfirst(isequal(ridx), spcfixID) != nothing
                    continue
                end

                # Build the loss expression and du[]
                netrates[ridx] *= " - $rcoef * $rate_expr"
                loss[ridx]     *= " + $rcoef * $rate_expr"
            end

            # Build prod
            for x in products
                # Build the loss and rxn rate.
                rm = match(rgx_coef, x)
                rcoef = (rm[:coef] == nothing ? 1 : parse(Float64, rm[:coef]))
                ridx = findfirst(isequal(rm[:chm]), spclist)
                if ridx == nothing
                    error("Fatal parser error at equation ", rxn, " -- could not find species ", rm[:chm])
                end
                # ridx = ridx[2]
                rspc = rm[:chm]

                # Check if fixed species... in spcfixID
                if findfirst(isequal(ridx), spcfixID) != nothing
                    continue
                end

                # Build the loss expression and du[]
                netrates[ridx] *= " + $rcoef * $rate_expr"
                prod[ridx]     *= " + $rcoef * $rate_expr"
            end
        end # driver == "Native"
    end # for each match

    if driver == "Native"
        # Build the ODE and Prod/Loss Diagn for the Native driver...
        codeeqn *= join(netrates, "\r\n") * "\r\n"
        prodloss_diag_kpp_eqn *= join(prod, "\r\n") * "\r\n\r\n"
        prodloss_diag_kpp_eqn *= join(loss, "\r\n") * "\r\n"

        # For cleaner code we can remove the initial "= 0 +" or "= 0 -" if they exist
        codeeqn = replace(codeeqn, "= 0 +" => "=")
        codeeqn = replace(codeeqn, "= 0 -" => "= -")
        prodloss_diag_kpp_eqn = replace(prodloss_diag_kpp_eqn, "= 0 +" => "=")
        prodloss_diag_kpp_eqn = replace(prodloss_diag_kpp_eqn, "= 0 -" => "= -")
    end

    modelparams["reactions"] = rxn_count

    println("Finished parsing equations file for a total of ", rxn_count, " reactions")

    # print(eqnstr)
    close(eqnio)
    codeeqn, prodloss_diag_kpp_eqn
end

"""
    readfile_def(spclist, spcfixID, spcIC, modelparams)
reads the .def KPP-file and extracts model parameters and species initial concentrations into respective array/dicts.
format of the .def file is IGNORED; very rough regex matching is performed for now.
all concs without ics will initialize with ZERO values
"""
function readfile_def(spclist::Array{String, 1}, spcfixID::Array{Int, 1}, spcIC::Array{Float64, 1}, modelparams::Dict{String, Any})
    println("Begin parsing $mechanism.def...")
    # Default values
    default_ic = 0.0
    cfactor = 1.0e0

    # Preallocate the spcIC array
    resize!(spcIC, length(spclist))

    defio = open("./inputs/$mechanism.def", "r")
    defstr = read(defio, String)

    # CFACTOR = 2.4476e+13;   ALL_SPEC = 0.0e0;   NO = 1.0e-1;   NO2 = 5.0e-2;
    # First, match hardcoded parameters CFACTOR, ALL_SPEC
    if((m = match(r"CFACTOR\s*=\s*(?<num>[0-9e\-\+\.]*)", defstr)) != nothing)
        cfactor = parse(Float64, m[:num])
    end

    if((m = match(r"ALL_SPEC\s*=\s*(?<num>[0-9e\-\+\.]*)", defstr)) != nothing)
        default_ic = parse(Float64, m[:num])
    end

    # Loop through species and retrieve their initial concentrations by index
    for (index, name) in enumerate(spclist)
        m = match(Regex("\\s+" * name * "\\s*=\\s*(?<num>[0-9e\\-\\+\\.]*)"), defstr)
        spcIC[index] = (m == nothing ? default_ic : parse(Float64, m[:num])) * cfactor
        # have to multiply by cfactor!

        # if spcIC[index] != default_ic
        #     println("Index $index ($name) matched ", m[:num], " parsed ", spcIC[index])
        #     println("Regex is ", Regex(name * "\\s*=\\s*(?<num>[0-9e\\-\\+\\.]*)"))
        # end
    end

    # Read model recommended parameters
    modelparams["CFACTOR"] = cfactor

    # The below acts are potentially dangerous...
    # match the C variants, which are semicolon-delimited
    if((m = match(r"TSTART\s*=[\s\(]*(?<num>[0-9e\-\+\.\\*\/]*)[\s\)]*;"i, defstr)) != nothing)
        modelparams["tstart"] = convert(Float64, eval(Meta.parse(m[:num])))
    end
    if((m = match(r"TEND\s*=\s*(?<tstart>TSTART\s*\+\s*)?[\s\(]*(?<num>[0-9e\-\+\.\\*\/]*)[\s\)]*;"i, defstr)) != nothing)
        modelparams["tend"] = (m[:tstart] == nothing ? 0 : modelparams["tstart"]) +
                              convert(Float64, eval(Meta.parse(m[:num])))
    end
    if((m = match(r"DT\s*=\s*(?<num>[0-9e\-\+\.\\*\/]*)\s*;"i, defstr)) != nothing)
        modelparams["dt"] = convert(Float64, eval(Meta.parse(m[:num])))
    end

    close(defio)
    println("Finished parsing $mechanism.def with model parameters:")
    display(modelparams)
end

function gen_ic(spclist::Array{String, 1}, spcfixID::Array{Int, 1}, spcIC::Array{Float64, 1}, modelparams::Dict{String, Any})
    registry_default_ic = "const jlkpp_bg_ic = Float64["
    kludge_fix_fixed_species = "\n"

    # Loop through species
    for (index, name) in enumerate(spclist)
        registry_default_ic *= string(spcIC[index], ", ")

        if driver == "DiffBioEq"
            # The "Native" driver does not need this hack but DiffBioEq does because
            # it cannot "fix" a species in time
            if index in spcfixID
                kludge_fix_fixed_species *= string("            chem_species[k,j,i][$index] = ", spcIC[index], "\n")
            end
        end
    end

    registry_default_ic *= "]\n"
    registry_default_ic, kludge_fix_fixed_species
end
