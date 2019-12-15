# KPP.jl

KPP.jl is a KPP (Kinetics PreProcessor)-compatible chemical kinetics model code generator for Julia. It interfaces with `DifferentialEquations.jl` and *optionally* the `DiffEqBiological.jl` to import mass-action kinetics equation descriptions to generate time-stepping models solved by the broader numerical computing infrastructure provided by Julia.

This package aims to provide KPP-like compatibility, in the sense that *pre-processing* for species and equations are provided and a model code abstracting the ODE generation, solving and time-stepping is generated according to descriptions in a format compatible with the original Kinetics PreProcessor (KPP), which generates FORTRAN, C and MATLAB code. This project is a **clean re-implementation** of KPP with no original code and no affiliation with the original authors, who may be cited below:

> Damian, V., Sandu, A., Damian, M., Potra, F. and Carmichael, G.R., 2002. The kinetic preprocessor KPP-a software environment for solving chemical kinetics. Computers & Chemical Engineering, 26(11), pp.1567-1579.

Not all features by KPP are supported and model code generated adheres to Julia standards, given KPP is more than a dozen years old. Please see the "Deprecated" section for details.

`KPP.jl` is developed by Haipeng Lin at `hplin as seas dot harvard dot edu` originally for the MIT 18.337 final project.

## Drivers
`KPP.jl` provides multiple "driver" packages enabling differing kinds of functionality. By default, it provides a simple `DiffEqBiological.jl`-based driver, `DiffBioEq`, which creates the reaction network using said package. KPP.jl is essentially a parser in this case.

For high performance applications `KPP.jl` incorporates the `Native` driver, which builds the species array manually and aims to support advanced developments in solving chemical equations specific to atmospheric chemistry:
* **Diagnostics:** Production/Loss for species at each time-step
* Adaptive approximations for reducing the complexity of the chemical mechanism, based on:
> Santillana M., P. Le Sager, D. J. Jacob, and M. P. Brenner, An adaptive reduction algorithm for efficient chemical calculations in global atmospheric chemistry models. Atmos. Environ., 44, 4426-4431, 2010

and possibly in the future an optimization as described in:

> Shen, L., D.J. Jacob, M. Santillana, X. Wang, and W. Chen, An adaptive method for speeding up the numerical integration of chemical mechanisms in atmospheric chemistry models: application to GEOS-Chem version 12.0.0, Geophys. Model Dev. Discuss., https://doi.org/10.5194/gmd-2019-2792019, in review, 2019.

Experimental features will also leverage `ModelingToolkit.generate_jacobian` to take the symbolic Jacobian, which can greatly improve performance. This is a work-in-progress.

## Feature description
TBD

## Differences from legacy KPP-generated code
Generated code is a full time-stepped model `jlkpp_MECHANISM_NAME` which performs time-stepping from start to end solving chemical kinetic equations step-by-step until completion.

* Variables are not handled in a global scope like in FORTRAN (COMMON blocks).

* Legacy parameters such as `TEMP`, `CFACTOR` and `SUN` have their replacement Julia routines / hard-coded shims for compatibility.

* Legacy rate functions such as `FALL`, `ARR`, `EP2` have hard-coded Julia version shims for compatibility. Currently functions for rate coefficients require work directly modifying `KPP.jl`. Submit an issue if you want to work on this.

* Extensibility is made possible by providing clear code points to implement other processes and a main non-KPP specific time stepping loop.

## Future features (wishlist)

* Debug IO: Using `JLD.jl`, a Julia-native format variant of HDF5 is used to save out debug output at specific time intervals **defined at pre-processor generation time**. It is fixed to discourage production use and only used for easy debugging of outputs.

## Unsupported features (caveats)
The following features are **currently unsupported** in `KPP.jl` and will be implemented as infrastructure is ready.

* Fixed species (e.g. `O2` or generic species like `M`). Their concentrations are not fixed in internal solver time-steps and are currently re-set with a kludge at every model time-step.

* As a result of the above, `#DEFVAR`, `#DEFFIX` may not behave as expected.

* Mass balance checks via chemical atoms (when defined via `.spc` files, e.g. `PAN		= 2C + 3H + 5O + N;`) is unsupported.

* The integrator option e.g. `#INTEGRATOR rodas4` is ignored; `DifferentialEquations.jl` is used with algorithm hinting `alg_hints = [:stiff]` by default.

## Deprecated features
The following features are **deprecated** and will remain indefinitely unsupported.

* The `#DRIVER` option specifying the language version is obviously unsupported. All code will be Julia code.

* The `#MONITOR` option to choose output species. The generated model does NOT contain production IO functionality. Users are encouraged to implement their own.

## Dependencies
We acknowledge the authors of the following packages (and Julia itself!) without which `KPP.jl` is not possible:

* `DifferentialEquations.jl`
* `DiffEqBiological.jl`

Additionally the demo `input` files are directly from `KPP-2.1` by original authors (Damien et al., 2002; Sandu et al., 2006).
