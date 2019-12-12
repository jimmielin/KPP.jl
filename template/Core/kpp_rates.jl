# KPP.jl: Implement reaction rate equations
# These were adapted from originals in KPP 2.0 (Damien et al., 2002) saprc99_Rates.f90

function jlkpp_ARR(TEMP::Float64, A0::Float64, B0::Float64, C0::Float64)
   A0 * exp(-B0/TEMP) * (TEMP/300.0)^C0
end

# Simplified Arrhenius, with two arguments
# Note: The argument B0 has a changed sign when compared to ARR
function jlkpp_ARR2(TEMP::Float64, A0::Float64, B0::Float64)
   A0 * exp(B0/TEMP)
end

function jlkpp_EP2(TEMP::Float64, CFACTOR::Float64,
                   A0::Float64, C0::Float64,
                   A2::Float64, C2::Float64,
                   A3::Float64, C3::Float64)
   K0 = A0 * exp(-C0/TEMP)
   K2 = A2 * exp(-C2/TEMP)
   K3 = A3 * exp(-C3/TEMP) * CFACTOR * 1.0e6
   K0 + K3/(1.0+K3/K2)
end

function jlkpp_EP3(TEMP::Float64, CFACTOR::Float64,
                   A1::Float64, C1::Float64,
                   A2::Float64, C2::Float64)
   K1 = A1 * exp(-C1/TEMP)
   K2 = A2 * exp(-C2/TEMP)
   K1 + K2*(1.0e6*CFACTOR)
end

function jlkpp_FALL(TEMP::Float64, CFACTOR::Float64,
                    A0::Float64, B0::Float64, C0::Float64,
                    A1::Float64, B1::Float64, C1::Float64,
                    CF::Float64)
   K0 = A0 * exp(-B0/TEMP) * (TEMP/300.0)^C0 * CFACTOR * 1.0e6
   K1 = K0/(A1 * exp(-B1/TEMP) * (TEMP/300.0)^C1)
   (K0/(1.0+K1))*(CF^(1.0/(1.0+(log10(K1))^2)))
end
