# KPP.jl

KPP.jl is a KPP (Kinetics PreProcessor)-compatible chemical kinetics model code generator for Julia. It interfaces with `DifferentialEquations.jl` and the `DiffEqBiological.jl` to import mass-action kinetics equation descriptions to generate time-stepping models solved by the broader numerical computing infrastructure provided by Julia.

This package aims to provide KPP-like compatibility, in the sense that *pre-processing* for species and equations are provided and a model code abstracting the ODE generation, solving and time-stepping is generated according to descriptions in a format compatible with the original Kinetics PreProcessor (KPP), which generates FORTRAN, C and MATLAB code. This project is a **clean re-implementation** of KPP with no original code and no affiliation with the original authors, who may be cited below:

> Damian, V., Sandu, A., Damian, M., Potra, F. and Carmichael, G.R., 2002. The kinetic preprocessor KPP-a software environment for solving chemical kinetics. Computers & Chemical Engineering, 26(11), pp.1567-1579.

Not all features by KPP are supported and model code generated adheres to Julia standards, given KPP is more than a dozen years old. Please see the "Deprecated" section for details.

`KPP.jl` is developed by Haipeng Lin at `hplin as seas dot harvard dot edu` originally for the MIT 18.337 final project.

## Feature description
TBD

## Unsupported features (caveats)
The following features are **currently unsupported** in `KPP.jl` and will be implemented as infrastructure is ready.

* Fixed species (e.g. `O2` or generic species like `M`). Their concentrations are not fixed in internal solver time-steps and are currently re-set with a kludge at every model time-step.

## Deprecated features
The following features are **deprecated** and will remain indefinitely unsupported.
