module constants
# import Unitful as u
# import UnitfulAstro as ua

const G = 4.30091727003628e-6 # Gravitational constant in kpc * km^2 / Msun / s^2, CODATA2018
const G_CGS = 6.6743e-8 # Gravitational constant in cm^3 / g / s^2, CODATA2018
const KPC_CM = 3.0856775814913673e21 # kiloparsec in cm, CODATA2018
const KPC_KM = 3.0856775814913673e16 # kiloparsec in cm, CODATA2018
const DELTA_COLLAPSE = 1.686470199841145 # linear overdensity threshold for halo collapse for the spherical top-hat collapse model, 3/5*(3π/2)^(2/3)
const RHO_C_Z0_CGS = 1.878341616932337e-29 # * u.g / u.cm^3 # critical density at z=0 in gram / cm^3
const RHO_C_Z0_KPC = 2.77536627245708E2 # * ua.Msun / ua.kpc^3# critical density at z=0 in Msun / kpc^3
const RHO_C_Z0_MPC = 2.7753662724570807e11 # * ua.Msun / ua.kpc^3# critical density at z=0 in Msun / Mpc^3


end # module
