module constants
# import Unitful as u
# import UnitfulAstro as ua

""" Newton's gravitational constant in units of kpc km^2 / Msun / s^2. """
const G = 4.30091727003628e-6 # Gravitational constant in kpc * km^2 / Msun / s^2, CODATA2018
""" Newton's gravitational constant in units of cm^3 / g / s^2. """
const G_CGS = 6.6743e-8 # Gravitational constant in cm^3 / g / s^2, CODATA2018
""" Number of centimeters in a kiloparsec. """
const KPC_CM = 3.0856775814913673e21 # kiloparsec in cm, CODATA2018
""" Number of kilometers in a kiloparsec. """
const KPC_KM = 3.0856775814913673e16 # kiloparsec in km, CODATA2018
""" Linear overdensity threshold for halo collapse for the spherical top-hat collapse model, ``3/5*(3π/2)^{2/3}``. """
const DELTA_COLLAPSE = 1.686470199841145 # linear overdensity threshold for halo collapse for the spherical top-hat collapse model, 3/5*(3π/2)^(2/3)
""" Prefactor for the critical density of the Universe, in units of g / cm^3,
```math
\\frac{3}{8πG} \\ \\left(100 \\ \\text{km} \\ \\text{s}^{-1} \\ \\text{Mpc}^{-1} \\right)^2
```
"""
const RHO_C_Z0_CGS = 1.878341616932337e-29 # * u.g / u.cm^3 # critical density at z=0 in gram / cm^3
""" Prefactor for the critical density of the Universe, in units of Msun / kpc^3,
```math
\\frac{3}{8πG} \\ \\left(100 \\ \\text{km} \\ \\text{s}^{-1} \\ \\text{Mpc}^{-1} \\right)^2
```
"""
const RHO_C_Z0_KPC = 2.77536627245708E2 # * ua.Msun / ua.kpc^3# critical density at z=0 in Msun / kpc^3
""" Prefactor for the critical density of the Universe, in units of Msun / Mpc^3,
```math
\\frac{3}{8πG} \\ \\left(100 \\ \\text{km} \\ \\text{s}^{-1} \\ \\text{Mpc}^{-1} \\right)^2
```
"""
const RHO_C_Z0_MPC = 2.7753662724570807e11 # * ua.Msun / ua.kpc^3# critical density at z=0 in Msun / Mpc^3
""" Conversion factor between CMB temperature and neutrino temperature, `` \\frac{4}{11}^{1/3}``. """
const TNU_PREFAC = 0.7137658555036082 # Equivalent to (4//11)^(1/3); for neutrino temperature
# const TNU_PREFAC = big"0.7137658555036081706718999917626612475907965475890380669156267520845831470677149" # Equivalent to (4//11)^(1/3); for neutrino temperature
const ONU_PREFACC = 7//8 * big(4//11)^big(4//3) # Prefactor for neutrino density Ω_ν (Omega_nu)

end # module
