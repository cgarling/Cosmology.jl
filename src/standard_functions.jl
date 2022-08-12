#############################################################################################
# Type handling and parameter retrieval
#############################################################################################
"""
    h(c::T) where T<:AbstractCosmology
Return "little h" for the cosmology, defined as the Hubble constant at present day (H0) divided by 100 km / s / Mpc.
"""
h(c::T) where T<:AbstractCosmology = c.h
"""
    partype(c::T) where T<:AbstractCosmology
Return the type of the scalar fields in the `AbstractCosmology` object `c`. 
"""
partype(c::T) where T<:AbstractCosmology = typeof(h(c))
"""
    m_nu(c::T) where T<:AbstractCosmology
Return the masses of the neutrino species in the cosmology `c` in eV.
"""
m_nu(c::T) where T<:AbstractCosmology = c.m_nu
"""
    Neff(c::T) where T<:AbstractCosmology
Return the effective number of neutrino species in the cosmology `c`.
"""
Neff(c::T) where T<:AbstractCosmology = c.Neff
"""
    w0(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM})
Return the present-day value of `w` for the CDL dark energy equation of state.
"""
w0(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM}) = c.w0
w0(c::Union{FlatLCDM,OpenLCDM,ClosedLCDM}) = -1
"""
    wa(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM})
Return the evolution parameter for the CDL dark energy equation of state.
"""
wa(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM}) = c.wa
wa(c::Union{FlatLCDM,OpenLCDM,ClosedLCDM}) = 0

#############################################################################################
# Neutrinos and temperatures
#############################################################################################
"""
    T_cmb([u::Unitful.Unitlike,] c::AbstractCosmology)
    T_cmb([u::Unitful.Unitlike,] c::AbstractCosmology, z)
    T_cmb([u::Unitful.Unitlike,] z, Tcmb0)
Return the temperature of the CMB at redshift `z`, in Kelvin. Will convert to compatible unit `u` if provided.

```math
\\text{T}_\\text{CMB}(z) = \\text{T}_\\text{CMB,0} \\times \\left( 1+z \\right)
```
"""
T_cmb(c::AbstractCosmology) = c.Tcmb0 * u.K
T_cmb(c::AbstractCosmology, z) = T_cmb(c) * (1 + z) # (T=T_cmb(c); muladd(T,z,T)) # This muladd is not any faster as far as I can tell
T_cmb(z,Tcmb0::u.Temperature) = u.uconvert(u.K,Tcmb0) * (1 + z)
T_cmb(z,Tcmb0) = Tcmb0 * (1 + z) * u.K

"""
    n_nu(c::AbstractCosmology)
    n_nu(Neff::Real)
Return the number of discrete neutrinos in the cosmology, defined as `Int(floor(Neff(c)))`; see [`Neff`](@ref Cosmology.Neff).
"""
n_nu(c::AbstractCosmology) = Int(floor(Neff(c))) #number of neutrinos
n_nu(Neff::Real) = Int(floor(Neff))


"""
    T_nu([u::Unitful.Unitlike,] c::AbstractCosmology [,z::Number])
    T_nu(Tcmb0::Number [,z::Number])
    T_nu(Tcmb0::u.Temperature [,z::Number])
Return the neutrino temperature of the cosmology at redshift `z` in Kelvin. If `z` is not provided, will return neutrino temperature at `z=0`. Will convert to compatible unit `u` if provided.

```math
\\begin{aligned}
\\text{T}_{\\nu,0} &= \\frac{4}{11}^{\\left(1/3\\right)} \\times \\text{T}_\\text{CMB,0} \\newline
\\text{T}_\\nu(z) &= \\text{T}_{\\nu,0} \\times \\left( 1+z \\right) 
\\end{aligned}
```

See, e.g., Equation 3.1.20 on page 154 of Weinberg's "Cosmology" for an explanation of the prefactor.
"""
T_nu(Tcmb0::T) where T<:Number = T(constants.TNU_PREFAC) * Tcmb0 * u.K
T_nu(Tcmb0::T, z::T) where T<:Number = T_nu(Tcmb0) * (1 + z)
T_nu(Tcmb0::Number, z::Number) = T_nu(promote(Tcmb0,z)...)
T_nu(Tcmb0::u.Temperature, z::Number) = T_nu(u.ustrip(u.K,Tcmb0),z)
T_nu(Tcmb0::u.Temperature) = T_nu(u.ustrip(u.K,Tcmb0))
T_nu(c::AbstractCosmology) = T_nu(T_cmb(c))
T_nu(c::AbstractCosmology,z) = T_nu(T_cmb(c,z))

""" 
    nu_relative_density(m_nu::Number, Neff::Number, nu_temp::Number N_nu::Number=n_nu(Neff))
    nu_relative_density(m_nu, Neff::Number, nu_temp::Number, N_nu::Union{Nothing,Number}=nothing)
    nu_relative_density(c::AbstractCosmology, z)
    nu_relative_density(c::AbstractCosmology)
Return the neutrino density function relative to the energy density in photons. If `!(m_nu isa Number)`, then `m_nu` should be iterable and indexable. When called with an `AbstractCosmology` but without a redshift, returns the `z=0` value. 

# Arguments
 - `m_nu::Any`; either a `Number` or an iterable (like an `Array` or `Tuple`) that contains the neutrino masses in eV.
 - `Neff`; effective number of neutrino species; see [`Neff`](@ref Cosmology.Neff).
 - `N_nu`; number of neutrino species; see [`n_nu`](@ref Cosmology.n_nu).
 - `nu_temp`; temperature of neutrino background in Kelvin; see [`T_nu`](@ref Cosmology.T_nu). This is the argument that carries the redshift dependence. """
@inline function nu_relative_density(m_nu::T, Neff::T, nu_temp::T, N_nu::S) where {T<:Number,S<:Integer}
    (iszero(Neff) || iszero(nu_temp)) && (return zero(T))
    prefac = T(0.22710731766023898) #7/8 * (4/11)^(4/3)
    iszero(m_nu) && (return prefac * Neff)
    p = T(1.83)
    invp = T(0.54644808743)  # 1.0 / p
    k = T(0.3173)
    nu_y = m_nu / (T(8.617333262145179e-5) * nu_temp )
    rel_mass = (1 + (k * nu_y)^p)^invp
    return prefac * m_nu * Neff / N_nu
end
nu_relative_density(m_nu::Number, Neff::Number, nu_temp::Number, N_nu::Number=n_nu(Neff)) = nu_relative_density(promote(m_nu,Neff,nu_temp)...,Int(N_nu))
# Now for if m_nu isa Vector, Tuple, etc. 
@inline function nu_relative_density(m_nu, Neff::T, nu_temp::T, N_nu::S) where {T<:Number,S<:Integer}
    (iszero(Neff) || iszero(nu_temp)) && (return zero(T))
    prefac = T(0.22710731766023898) #7/8 * (4/11)^(4/3)
    p = T(1.83)
    invp = T(0.54644808743)  # 1 / p
    k = T(0.3173)
    nu_y_fac = T(11604.518121550082) # 1 / 8.617333262145179e-5
    inv_nu_temp = 1 / nu_temp
    rel_mass_per_sum = zero(T)
    massless=0
    for i in eachindex(m_nu)
        if iszero(m_nu[i])
            massless+=1
        else
            # nu_y = m_nu[i] / (8.617333262145179e-5 * nu_temp )
            nu_y = m_nu[i] * nu_y_fac * inv_nu_temp
            rel_mass_per_sum += (1 + (k * nu_y)^p)^invp
        end
    end
    return prefac * ( rel_mass_per_sum + massless ) * Neff / N_nu
end
nu_relative_density(m_nu, Neff::Number, nu_temp::Number, N_nu::Number=n_nu(Neff)) = nu_relative_density(m_nu,promote(Neff,nu_temp)...,Int(N_nu))
nu_relative_density(c::AbstractCosmology, z) = nu_relative_density(m_nu(c), Neff(c), u.ustrip(u.K,T_nu(c,z)))
nu_relative_density(c::AbstractCosmology) = nu_relative_density(m_nu(c), Neff(c), u.ustrip(u.K,T_nu(c)))
##############################################################################
"""
    a2E(c::AbstractCosmology,a)
    a2E(a,OmegaM,OmegaK,OmegaL,OmegaG,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Return the cosmological [`E`](@ref) factor times the square of the [scale factor](@ref Cosmology.scale_factor) `a`.
"""
function a2E(c::FlatLCDM, a)
    z = 1/a-1
    Or = Ω_γ(c) * (1 + nu_relative_density(c,z))
    a2 = a * a
    return sqrt(Or + Ω_m(c) * a + Ω_Λ(c) * a2 * a2) 
end

function a2E(c::Union{ClosedLCDM,OpenLCDM}, a)
    z = 1/a-1
    Or = Ω_γ(c) * (1 + nu_relative_density(c,z))
    a2 = a * a
    sqrt(Or + Ω_m(c) * a + (Ω_k(c) + Ω_Λ(c) * a2) * a2)
end

function a2E(c::Union{FlatWCDM,ClosedWCDM,OpenWCDM}, a)
    z = 1/a-1
    Or = Ω_γ(c) * (1 + nu_relative_density(c,z))
    ade = exp((1 - 3 * (w0(c) + wa(c))) * log(a) + 3 * wa(c) * (a - 1))
    sqrt(Or + (Ω_m(c) + Ω_k(c) * a) * a + Ω_Λ(c) * ade)
end

function a2E(a,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
    z = 1/a-1
    OmegaG = 4.481620089297254e-7 * Tcmb0^4 / h^2
    nu_temp = u.ustrip(u.K,T_nu(Tcmb0, z))
    Or = OmegaG * (1+nu_relative_density(m_nu, Neff, nu_temp))
    ade = exp((1 - 3 * (w0 + wa)) * log(a) + 3 * wa * (a - 1))
    sqrt(Or + (OmegaM + OmegaK * a) * a + OmegaL * ade)
end

###############################################################################
# Hubble rate
###############################################################################
""" 
    scale_factor(z::Real) or scale_factor(c::AbstractCosmology, z::Real)

Calculate the scale factor at redshift `z`. The scale factor is defined as ``a=\\frac{1}{1+z}``.
The derivative with respect to z is available as [`∇scale_factor`](@ref). 
"""
scale_factor(z::Real) = 1 / (1 + z)
scale_factor(c::AbstractCosmology, z::Real) = scale_factor(z) #for compatibility with z_at_value
""" 
    ∇scale_factor(z::Real) or ∇scale_factor(c::AbstractCosmology, z::Real)
Calculate the derivative of the scale factor at redshift `z` with respect to `z`; ``\\frac{da}{dz} = -\\frac{1}{\\left(1+z\\right)^2}``. """
∇scale_factor(z::Real) = -1 / (1 + z)^2
∇scale_factor(c::AbstractCosmology,z) = ∇scale_factor(z)

"""
    E(c::AbstractCosmology, z::Real)
    E(z,h,OmegaM,OmegaK,OmegaL,OmegaG,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Return the Hubble parameter as a function of redshift, in units of H₀. See also [`a2E`](@ref Cosmology.a2E), [`H`](@ref Cosmology.H). 
"""
E(c::AbstractCosmology, z::Real) = (a = scale_factor(z); a2E(c, a) / a^2)
E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = (a=scale_factor(z); a2E(a,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) / a^2)
"""
    H(c::AbstractCosmology, z::Real)
Return the Hubble parameter as a function of redshift, in km / s / Mpc. From the Friedmann equation,

```math
H^2 = \\left( \\frac{\\dot{a}}{a} \\right) = \\frac{8 \\pi G \\rho}{3} - \\frac{kc^2}{a^2}
```

where `a` is the cosmological [`scale factor`](@ref Cosmology.scale_factor), ``\\dot{a}`` is the time derivative of `a`, `G` is Newton's gravitational constant, `c` is the speed of light, ``\\rho`` is a mass density, and `k` is the curvature parameter. `k` is typically -1, 0, or 1; `k=-1` corresponds to an open Universe (forever expanding), `k=0` corresponds to a "flat" Universe with [`ρ==ρ_c`](@ref Cosmology.ρ_c), and `k=1` corresponds to a closed Universe (eventual contraction).
"""
H(c::AbstractCosmology, z::Real) = 100 * h(c) * E(c, z) * u.km / u.s / ua.Mpc
"""
    hubble_dist0(c::AbstractCosmology)
Calculate the Hubble distance at present-day in Mpc, defined as the speed of light times the Hubble time.
```math
    D_0 = \\frac{c}{H_0}
```
"""
hubble_dist0(c::AbstractCosmology) = partype(c)(2997.92458) / h(c) * ua.Mpc # constant is speed of light in km/s divided by 100
"""
    hubble_dist(c::AbstractCosmology, z::Real)
Return the speed of light times the Hubble time at redshift `z` in Mpc.
```math
    D(z) = \\frac{c}{H(z)}
```
"""
hubble_dist(c::AbstractCosmology, z::Real) = hubble_dist0(c) / E(c, z)

"""
    hubble_time(c::AbstractCosmology, z::Real)
Return ``\\frac{1}{\\text{H}_0}`` in units of Gyr.
"""
hubble_time0(c::AbstractCosmology) =  partype(c)(9.777922216807893) / h(c) * u.Gyr # 9.77814 # 9.777922216807893 / h(c) * u.Gyr
"""
    hubble_time(c::AbstractCosmology, z::Real)
Return ``\\frac{1}{\\text{H}\\left(z\\right)}`` in units of Gyr.
"""
hubble_time(c::AbstractCosmology, z::Real) = hubble_time0(c) / E(c, z)

###############################################################################
# Distances
###############################################################################

Z(c::AbstractCosmology, z::Real, ::Nothing; kws...) = quadgk(a->1 / a2E(c, a), scale_factor(z), 1; kws...)[1]
Z(c::AbstractCosmology, z₁::Real, z₂::Real; kws...) = quadgk(a->1 / a2E(c, a), scale_factor(z₂), scale_factor(z₁); kws...)[1]

"""
    comoving_radial_dist([u::Unitlike,] c::AbstractCosmology, [z₁,] z₂; kws...)

Calculate the [comoving radial distance](https://ned.ipac.caltech.edu/level5/Hogg/Hogg4.html) (sometimes also called comoving line-of-sight distance) in Mpc at redshift `z₂` as seen by an observer at `z₁`.  Redshift `z₁` defaults to 0 if omitted.  Will convert to compatible unit `u` if provided. `kws` are integration options passed to `QuadGK.quadgk`.
```math
\\begin{aligned}
D_R(z_1,z_2) = \\frac{c}{H_0} \\int_{z_1}^{z_2} \\frac{1}{E(z)} \\ dz = \\frac{c}{H_0} \\int_{a(z_2)}^{a(z_1)} \\frac{1}{E(a^{\\prime}) \\, a^{\\prime \\, 2}} \\ da^{\\prime}
\\end{aligned}
```
"""
comoving_radial_dist(c::AbstractCosmology, z₁, z₂ = nothing; kws...) = hubble_dist0(c) * Z(c, z₁, z₂; kws...)

"""
    comoving_transverse_dist(c::T<:AbstractCosmology, z₁, z₂ = nothing; kws...)

Returns the [comoving transverse distance](https://ned.ipac.caltech.edu/level5/Hogg/Hogg5.html) between two points with an angular separation of 1 radian at redshift `z₂` as measured by an observer at redshift `z₁`.
This is the same as the comoving distance if Ω_k is zero (as in the current concordance ΛCDM model). Will convert to compatible unit `u` if provided. `kws` are integration options passed to quadgk.
```math
D_T(z_1,z_2) = 
\\begin{cases}
\\frac{c}{H_0 \\ \\sqrt{\\Omega_{k,0}}} \\ \\text{sinh} \\left( \\frac{H_0 \\ \\sqrt{\\Omega_{k,0}}}{c} \\ D_R(z_1,z_2) \\right) \\ &\\forall \\ \\Omega_k > 0 \\newline
D_R(z_1,z_2) \\ &\\forall \\ \\Omega_k = 0 \\newline
\\frac{c}{H_0 \\ \\sqrt{-\\Omega_{k,0}}} \\ \\text{sin} \\left( \\frac{H_0 \\ \\sqrt{-\\Omega_{k,0}}}{c} \\ D_R(z_1,z_2) \\right) \\ &\\forall \\ \\Omega_k < 0 \\newline
\\end{cases}
```
where ``D_R(z_1,z_2)`` is the [`comoving radial distance`](@ref Cosmology.comoving_radial_dist).
"""
comoving_transverse_dist(c::AbstractFlatCosmology, z₁, z₂ = nothing; kws...) = comoving_radial_dist(c, z₁, z₂; kws...)
function comoving_transverse_dist(c::AbstractOpenCosmology, z₁, z₂ = nothing; kws...)
    sqrtok = sqrt(Ω_k(c))
    hubble_dist0(c) * sinh(sqrtok * Z(c, z₁, z₂; kws...)) / sqrtok
end
function comoving_transverse_dist(c::AbstractClosedCosmology, z₁, z₂ = nothing; kws...)
    sqrtok = sqrt(abs(Ω_k(c)))
    hubble_dist0(c) * sin(sqrtok * Z(c, z₁, z₂; kws...)) / sqrtok
end

"""
    angular_diameter_dist([u::Unitlike,] c::AbstractCosmology, [z₁,] z₂; kws...)

Ratio of the proper transverse size in Mpc of an object at redshift `z₂` to its angular size in radians, as seen by an observer at `z₁`.  Redshift `z₁` defaults to 0 if omitted.  Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
angular_diameter_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...) / (1 + z)
angular_diameter_dist(c::AbstractCosmology, z₁, z₂; kws...) =
    comoving_transverse_dist(c, z₁, z₂; kws...) / (1 + z₂)

"""
    luminosity_dist([u::Unitlike,] c::AbstractCosmology, z; kws...)

Bolometric luminosity distance in Mpc at redshift `z`. Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
luminosity_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...) * (1 + z)

"""
    distmod(c::AbstractCosmology, z; kws...)

Distance modulus in magnitudes at redshift `z`. kws are integration options passed to quadgk.
"""
distmod(c::AbstractCosmology, z; kws...) =
    5 * log10(luminosity_dist(c, z; kws...) / ua.Mpc) + 25

#######################################################################################
# Volumes
#######################################################################################

"""
    comoving_volume([u::Unitlike,] c::AbstractCosmology, z; kws...)

Comoving volume in cubic Gpc out to redshift `z`. Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
comoving_volume(c::AbstractFlatCosmology, z; kws...) =
    (4π / 3) * (comoving_radial_dist(ua.Gpc, c, z; kws...))^3
function comoving_volume(c::AbstractOpenCosmology, z; kws...)
    DH = hubble_dist0(ua.Gpc, c)
    x = comoving_transverse_dist(ua.Gpc, c, z; kws...) / DH
    Ok = Ω_k(c)
    sqrtok = sqrt(Ok)
    2π * DH^3 * (x * sqrt(1 + Ok * x^2) - asinh(sqrtok * x) / sqrtok) / Ok
end
function comoving_volume(c::AbstractClosedCosmology, z; kws...)
    DH = hubble_dist0(ua.Gpc, c)
    x = comoving_transverse_dist(ua.Gpc, c, z; kws...) / DH
    Ok = Ω_k(c)
    sqrtok = sqrt(abs(Ok))
    2π * DH^3 * (x * sqrt(1 + Ok * x^2) - asin(sqrtok * x) / sqrtok) / Ok
end

"""
    comoving_volume_element([u::Unitlike,] c::AbstractCosmology, z; kws...)

Comoving volume element in Gpc out to redshift `z`. Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
comoving_volume_element(c::AbstractCosmology, z; kws...) =
    hubble_dist0(ua.Gpc, c) * angular_diameter_dist(ua.Gpc, c, z; kws...)^2 / a2E(c, scale_factor(z))

#############################################################################################
# Times
#############################################################################################

T(c::AbstractCosmology, a0, a1; kws...) = quadgk(x->x / a2E(c, x), a0, a1; kws...)[1]

"""
    age([u::Unitlike,] c::AbstractCosmology, z; kws...)

Age of the universe in Gyr at redshift `z`. Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
age(c::AbstractCosmology, z; kws...) = hubble_time0(c) * T(c, 0, scale_factor(z); kws...)

"""
    lookback_time([u::Unitlike,] c::AbstractCosmology, z; kws...)

Difference between age at redshift 0 and age at redshift `z` in Gyr. Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
lookback_time(c::AbstractCosmology, z; kws...) = hubble_time0(c) * T(c, scale_factor(z), 1; kws...)

#############################################################################################
# Dark energy
#############################################################################################

""" 
    w(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM,FlatLCDM,OpenLCDM,ClosedLCDM},z)
    w(z::Real,w0::Real,wa::Real)
Evaluates the redshift dependence of the dark energy density. The scaling factor, ``I(z)``, is defined by ``ρ_Λ(z) = I(z) \\ ρ_{Λ,0}``.
"""
w(c::Union{FlatLCDM,OpenLCDM,ClosedLCDM},z) = -1
w(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM},z) = w0(c) + wa(c) * z / (1 + z)
w(z::Real,w0::Real,wa::Real) = w0 + wa * z / (1+z)

"""
    de_density_scale(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM},z)
    de_density_scale(z::Real,w0::Real,wa::Real)
Returns the redshift scaling of the dark energy density. See [`ρ_Λ`](@ref). 
"""
function de_density_scale(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM},z)
    zp1 = 1 + z
    return zp1^(3 * (1 + w0(c) + wa(c))) * exp(-3 * wa(c) * z / zp1)
end
de_density_scale(c::Union{FlatLCDM,OpenLCDM,ClosedLCDM},z) = 1
de_density_scale(z::Real,w0::Real,wa::Real) = (zp1 = 1+z; zp1^(3 * (1 + w0 + wa)) * exp(-3 * wa * z / zp1))

#############################################################################################
# Densities
#############################################################################################

"""
    ρ_c([u::UnitLike,], c::AbstractCosmology, z)
    ρ_c([u::UnitLike,], z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
    ρ_c(h::Number,E::Number)
The critical density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.

```math
\\begin{aligned}
\\rho_{0,c} &\\equiv \\frac{3 H_0^2}{8 \\pi G} = 1.878 \\times 10^{-29} \\ h^2 \\ \\text{g/cm}^3 \\newline
\\rho_{c}(z) &\\equiv \\frac{3 H(z)^2}{8 \\pi G} = 1.878 \\times 10^{-29} \\ \\left(E(z) \\times h\\right)^2 \\ \\text{g/cm}^3
\\end{aligned}
```
where [`E`](@ref Cosmology.E) is the Hubble factor in units of ``H_0``. See, e.g., Equation 1.5.28 on page 57 of Weinberg's "Cosmology" for more information.
"""
ρ_c(c::AbstractCosmology,z) = partype(c)(constants.RHO_C_Z0_CGS) * h(c)^2 * E(c,z)^2 * u.g/u.cm^3
ρ_c(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = constants.RHO_C_Z0_CGS * h^2 *
    E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2 * u.g/u.cm^3

"""
    ρ_m([u::UnitLike,], c::AbstractCosmology, z)
    ρ_m([u::UnitLike,], z,h,OmegaM)
The matter density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_m(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{m,0} \\ (1+z)^3 \\ \\text{g/cm}^3 
```
"""
ρ_m(c::AbstractCosmology,z) = partype(c)(constants.RHO_C_Z0_CGS) * h(c)^2 * Ω_m(c) * (1 + z)^3 * u.g/u.cm^3
ρ_m(z,h,Ω_m) = constants.RHO_C_Z0_CGS * h^2 * Ω_m * (1+z)^3 * u.g/u.cm^3

"""
    ρ_b([u::UnitLike,], c::AbstractCosmology, z)
    ρ_b([u::UnitLike,], z,h,OmegaB)
The baryon density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_b(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{b,0} \\ (1+z)^3 \\ \\text{g/cm}^3 
```
"""
ρ_b(c::AbstractCosmology,z) = partype(c)(constants.RHO_C_Z0_CGS) * h(c)^2 * Ω_b(c) * (1 + z)^3 * u.g/u.cm^3
ρ_b(z,h,Ω_b) = constants.RHO_C_Z0_CGS * h^2 * Ω_b * (1+z)^3 * u.g/u.cm^3

"""
    ρ_dm([u::UnitLike,], c::AbstractCosmology, z)
    ρ_dm([u::UnitLike,], z,h,OmegaDM)
The dark matter density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_{dm}(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{dm,0} \\ (1+z)^3 \\ \\text{g/cm}^3 
```
"""
ρ_dm(c::AbstractCosmology,z) = partype(c)(constants.RHO_C_Z0_CGS) * h(c)^2 * Ω_dm(c) * (1 + z)^3 * u.g/u.cm^3
ρ_dm(z,h,Ω_dm) = constants.RHO_C_Z0_CGS * h^2 * Ω_dm * (1+z)^3 * u.g/u.cm^3

"""
    ρ_Λ([u::UnitLike,], c::AbstractCosmology, z)
    ρ_Λ([u::UnitLike,], z,h,Ω_Λ,w0=-1,wa=0)
The dark energy density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_\\Lambda(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{\\Lambda,0} \\ \\text{g/cm}^3 
```
for a cosmological constant. For more general dark energy equations of state, see [`de_density_scale`](@ref Cosmology.de_density_scale).
"""
ρ_Λ(c::AbstractCosmology,z) = partype(c)(constants.RHO_C_Z0_CGS) * h(c)^2 * Ω_Λ(c) / de_density_scale(c,z) * u.g/u.cm^3
ρ_Λ(z,h,Ω_Λ,w0=-1,wa=0) = constants.RHO_C_Z0_CGS * h^2 * Ω_Λ / de_density_scale(z,w0,wa) * u.g/u.cm^3

"""
    ρ_γ([u::UnitLike,], c::AbstractCosmology, z)
    ρ_γ([u::UnitLike,], z,h,Ω_γ)
The photon matter density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_\\gamma(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{\\gamma,0} \\ (1+z)^4 \\ \\text{g/cm}^3 
```
"""
ρ_γ(c::AbstractCosmology,z) = partype(c)(constants.RHO_C_Z0_CGS) * h(c)^2 * Ω_γ(c) * (1 + z)^4 * u.g/u.cm^3
ρ_γ(z,h,Ω_γ) = constants.RHO_C_Z0_CGS * h^2 * Ω_γ * (1+z)^4 * u.g/u.cm^3

"""
    ρ_ν([u::UnitLike,], c::AbstractCosmology, z)
    ρ_ν([u::UnitLike,], z,h,Tcmb0,Neff,m_nu,N_nu=nothing)
The neutrino energy density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_\\nu(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{\\nu,0} \\ (1+z)^4 \\ \\text{g/cm}^3 
```
"""
ρ_ν(c::AbstractCosmology,z) = ρ_γ(c,z) * nu_relative_density(c, z)
ρ_ν(z,h,Tcmb0,Neff,m_nu,N_nu=nothing) = (Ω_γ = 4.481620089297254e-7 * Tcmb0^4 / h^2; ρ_γ(z,h,Ω_γ) * nu_relative_density(m_nu, Neff, u.ustrip(T_nu(Tcmb0, z)), N_nu) )

"""
    ρ_r([u::UnitLike,], c::AbstractCosmology, z)
    ρ_r([u::UnitLike,], z,h,Tcmb0,Neff,m_nu,N_nu=nothing)
The energy density of the universe in relativistic species at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_r(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{r,0} \\ (1+z)^4 \\ \\text{g/cm}^3 
```
"""
ρ_r(c::AbstractCosmology,z) = ρ_γ(c,z) * (1 + nu_relative_density(c, z))
ρ_r(z,h,Tcmb0,Neff,m_nu,N_nu=nothing) = (Ω_γ = 4.481620089297254e-7 * Tcmb0^4 / h^2; ρ_γ(z,h,Ω_γ) * (1.0 + nu_relative_density(m_nu, Neff, u.ustrip(T_nu(Tcmb0, z)), N_nu) ) )


#############################################################################################
# Omegas
#############################################################################################

""" 
    Ω_m(c::AbstractCosmology,z)
    Ω_m(c::AbstractCosmology)
    Ω_m(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Density of matter at redshift `z` in units of the critical density. When called without a redshift, returns the `z=0` value. """
Ω_m(c::AbstractCosmology) = c.Ω_m
Ω_m(c::AbstractCosmology,z) = Ω_m(c) * (1 + z)^3 / E(c,z)^2
Ω_m(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = OmegaM * (1+z)^3 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2
""" 
    Ω_b(c::AbstractCosmology,z)
    Ω_b(c::AbstractCosmology)
    Ω_b(z,h,OmegaM,OmegaB,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Density of baryons at redshift `z` in units of the critical density. When called without a redshift, returns the `z=0` value. """
Ω_b(c::AbstractCosmology) = c.Ω_b
Ω_b(c::AbstractCosmology,z) = Ω_b(c) * (1 + z)^3 / E(c,z)^2
Ω_b(z,h,OmegaM,OmegaB,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = OmegaB * (1+z)^3 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2
""" 
    Ω_dm(c::AbstractCosmology,z)
    Ω_dm(c::AbstractCosmology)
    Ω_dm(z,h,OmegaM,OmegaB,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Density of dark matter at redshift `z` in units of the critical density. When called without a redshift, returns the `z=0` value. """
Ω_dm(c::AbstractCosmology) = Ω_m(c) - Ω_b(c)
Ω_dm(c::AbstractCosmology,z) = Ω_dm(c) * (1 + z)^3 / E(c,z)^2
Ω_dm(z,h,OmegaM,OmegaB,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = (OmegaM - OmegaB) * (1+z)^3 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2
""" 
    Ω_k(c::AbstractCosmology,z)
    Ω_k(c::AbstractCosmology)
    Ω_k(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density of curvature at redshift `z` in units of the critical density. When called without a redshift, returns the `z=0` value.
```math
\\Omega_k = -\\frac{k}{H_0^2}
```
where `k` is the curvature parameter in the Friedmann equation. See Equations  1.5.19 and 1.5.40 on pages 56 and 60 in Weinberg's "Cosmology". """
Ω_k(c::AbstractCosmology) = c.Ω_k
Ω_k(c::AbstractCosmology,z) = Ω_k(c) * (1 + z)^2 / E(c,z)^2
Ω_k(c::AbstractFlatCosmology,z=0.0) = zero(z)
Ω_k(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = OmegaK * (1+z)^2 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2
""" 
    Ω_γ(c::AbstractCosmology,z)
    Ω_γ(c::AbstractCosmology)
    Ω_γ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density of photons at redshift `z` in units of the critical density. Calculated from [`T_cmb`](@ref). When called without a redshift, returns the `z=0` value. """
Ω_γ(c::AbstractCosmology) = partype(c)(4.481620089297254e-7) * u.ustrip(u.K,T_cmb(c))^4 / h(c)^2
Ω_γ(c::AbstractCosmology,z) = Ω_γ(c) * (1 + z)^4 / E(c,z)^2
Ω_γ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = (OmegaG = 4.481620089297254e-7 * Tcmb0^4 / h^2; OmegaG * (1+z)^4 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2 )
""" 
    Ω_ν(c::AbstractCosmology,z)
    Ω_ν(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density in neutrinos at redshift `z` in units of the critical density. When called without a redshift, returns the `z=0` value. """
Ω_ν(c::AbstractCosmology) = Ω_γ(c) * nu_relative_density(c)
Ω_ν(c::AbstractCosmology,z) = Ω_γ(c,z) * nu_relative_density(c,z)
Ω_ν(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = (nu_temp = u.ustrip(u.K,T_nu(Tcmb0, z)); Ω_γ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) * nu_relative_density(m_nu, Neff, nu_temp) )
""" 
    Ω_r(c::AbstractCosmology,z)
    Ω_r(c::AbstractCosmology)
    Ω_r(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density in radiation at redshift `z` in units of the critical density. Evaluated as `Ω_γ(c,z) + Ω_ν(c,z)`; sum of photons and neutrinos. When called without a redshift, returns the `z=0` value. """
Ω_r(c::AbstractCosmology) = Ω_γ(c) + Ω_ν(c)
Ω_r(c::AbstractCosmology,z) = Ω_γ(c,z) + Ω_ν(c,z)
Ω_r(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = Ω_γ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) + Ω_ν(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)
""" 
    Ω_Λ(c::AbstractCosmology,z)
    Ω_Λ(c::AbstractCosmology)
    Ω_Λ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density in dark energy at redshift `z` in units of the critical density. When called without a redshift, returns the `z=0` value. """
Ω_Λ(c::AbstractCosmology) = c.Ω_Λ
Ω_Λ(c::AbstractCosmology,z) = Ω_Λ(c) * de_density_scale(c,z) / E(c,z)^2
Ω_Λ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = OmegaL * de_density_scale(z,w0,wa) / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2

#############################################################################################
# Roots
#############################################################################################

"""
    z_at_value(c::AbstractCosmology, func::Function, fval; zmin=1e-8, zmax=1000.0,ztol=1e-8,maxfun=500)
    Find the redshift ``z`` at which ``func(z) = fval``.

    This finds the redshift at which one of the cosmology functions or
    methods (for example Planck13.distmod) is equal to a known value.

    .. warning::
      Make sure you understand the behavior of the function that you
      are trying to invert! Depending on the cosmology, there may not
      be a unique solution. For example, in the standard Lambda CDM
      cosmology, there are two redshifts which give an angular
      diameter distance of 1500 Mpc, z ~ 0.7 and z ~ 3.8. To force
      ``z_at_value`` to find the solution you are interested in, use the
      ``zmin`` and ``zmax`` keywords to limit the search range (see the
      example below).

    Parameters
    ----------
    func : function or method
       A function that takes a redshift as input.
    fval : astropy.Quantity instance
       The value of ``func(z)``.
    zmin : float, optional
       The lower search limit for ``z``.  Beware of divergences
       in some cosmological functions, such as distance moduli,
       at z=0 (default 1e-8).
    zmax : float, optional
       The upper search limit for ``z`` (default 1000).
    kws : keyword arguments to be passed to Roots.find_zero

    Returns
    -------
    z : float
      The redshift ``z`` satisfying ``zmin < z < zmax`` and ``func(z) =
      fval`` within ``ztol``.
"""
function z_at_value(c::AbstractCosmology, func::Function, fval; zmin=1e-8, zmax=1e5, kws...)
    fval_zmin = func(c,zmin)
    fval_zmax = func(c,zmax)
    if fval isa u.Quantity && !(fval_zmin isa u.Quantity)
        throw(DomainError(fval,"fval is a Unitful Quantity but the provided function does not provide units."))
        return
    elseif fval isa u.Quantity && fval_zmin isa u.Quantity
        fval = fval |> u.unit(fval_zmin) |> u.ustrip
        f = (z) -> u.ustrip(func(c,z)) - fval
    elseif !(fval isa u.Quantity) && fval_zmin isa u.Quantity
        f = (z) -> u.ustrip(func(c,z)) - fval
    else
        f = (z) -> func(c,z) - fval
    end

    result = find_zero(f,(zmin,zmax), Bisection(); kws...) # atol=ztol, maxevals=maxfun)
    return result
end

#########################################################################################
# Misc
#########################################################################################

"""
    sound_horizon(c::AbstractCosmology)
Return the sound horizon length (in Mpc), given by Equation 26 in Eisenstein & Hu 1998. 
"""
sound_horizon(c::AbstractCosmology) = (h2 = h(c)^2; 44.5 * log(9.83 / Ω_m(c) / h2) / sqrt(1.0 + 10.0 * (Ω_b(c) * h2)^0.75) * ua.Mpc)
# sound_horizon(c::AbstractCosmology) = 44.5 * log(9.83 / c.Ω_m / c.h^2) / sqrt(1.0 + 10.0 * (c.Ω_b * c.h^2)^0.75) * ua.Mpc

"""
    matter_radiation_equality(c::AbstractCosmology)
Return the redshift of matter-radiation equality. This is computed as Equation 2 in Eisenstein and Hu 1998. I previously had a different formula but couldn't figure out where it came from. 
"""
matter_radiation_equality(c::AbstractCosmology) = 2.5e4 * Ω_m(c) * h(c)^2 / (u.ustrip(u.K,T_cmb(c))/2.7)^4
# matter_radiation_equality(c::AbstractCosmology) = 3600 * (Ω_m(c) * h(c)^2 / 0.15) - 1

#########################################################################################
# Easily select a different unit
#############################################################################################

for f in (:hubble_dist0, :hubble_dist, :hubble_time0, :hubble_time, :comoving_radial_dist,
          :comoving_transverse_dist, :angular_diameter_dist, :luminosity_dist,
          :comoving_volume, :comoving_volume_element, :age, :lookback_time,
          :T_nu, :T_cmb, :ρ_c, :ρ_m, :ρ_b, :ρ_dm, :ρ_Λ, :ρ_γ, :ρ_ν, :ρ_r, :sound_horizon)
    @eval $f(uu::u.Unitlike, args...; kws...) = u.uconvert(uu, $f(args...; kws...))
end
