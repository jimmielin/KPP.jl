# KPP.jl optimizer
# Adaptive reduction mechanism by Santillana et al., 2010
# Reference: Santillana M., P. Le Sager, D. J. Jacob, and M. P. Brenner, An adaptive reduction algorithm for efficient chemical calculations in global atmospheric chemistry models. Atmos. Environ., 44, 4426-4431, 2010

# Implemented by Haipeng Lin <hplin@seas.harvard.edu> for KPP.jl
# Atmospheric Chemistry Modeling Group, Harvard University
# 11 December 2019

function jlkpp_adapt_diagn(u::Array{Float64, 1}, t::Float64)
    # Set net rate production/loss delta threshold parameter
    # Need to vary this parameter depending on units system.
    δ = 100

    # Compute production and loss diagnostics (p is an approximation)
    # FIXME hardcoded for now
    p = t -> (jlkpp_SUN(t), 300.0, jlkpp_CFACTOR)
    jlkpp_mechanism_prodloss!(jlkpp_scratch_p, jlkpp_scratch_l, u, p, t)

    # Outs
    _optim_slow_count = 0
    for spc in 1:jlkpp_nspecies
        jlkpp_adapt_slow[spc] = (jlkpp_scratch_p[spc] < δ && jlkpp_scratch_l[spc] < δ)
        jlkpp_adapt_slow[spc] && (_optim_slow_count += 1)
    end

    _optim_slow_count
end

function jlkpp_adapt_slow!(u::Array{Float64, 1}, dt::Float64)
    # Update in-place species using explicit solution
    #          Pi*           Pi*    -ki*Δt
    # C[n+1] = --- + (C[n] - ---) e
    #          ki*           ki*
    # Where we assume first-order loss Pi-kiCi. We don't have that kind of loss information
    # so we just assume ki* = Li/Ci, which is fair enough as most equations are like this.

    # We use the global scratchpad
    for spc in 1:jlkpp_nspecies
        if jlkpp_adapt_slow[spc] # is slow spc
            if jlkpp_scratch_l[spc] != 0
                # print("* jlkpp_adapt_slow! updating slow spc ", spc, " from ", u[spc])
                # print(" p/l=", jlkpp_scratch_p[spc]/jlkpp_scratch_l[spc], " exp=", exp(-jlkpp_scratch_l[spc]*dt/u[spc]))
                u[spc] = (jlkpp_scratch_p[spc]/jlkpp_scratch_l[spc] + (1.0-jlkpp_scratch_p[spc]/jlkpp_scratch_l[spc])*exp(-jlkpp_scratch_l[spc]*dt/u[spc]))*u[spc]
                # println("-> to ", u[spc])
            else
                # print("* jlkpp_adapt_slow! updating slow spc ", spc, " from ", u[spc])
                u[spc] = jlkpp_scratch_p[spc]*dt + u[spc]
                # println("-> using linear method to ", u[spc])
            end

            # Never allow negative species conc
            u[spc] < 0 && (u[spc] = 0)
        end
    end
end
