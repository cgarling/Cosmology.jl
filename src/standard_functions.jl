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

# Examples
```jldoctest
julia> T_cmb(Cosmology.Planck18) ≈ 2.7255 * Unitful.K
true

julia> T_cmb(Cosmology.Planck18, 1.0) ≈ 5.451 * Unitful.K
true

julia> T_cmb(Unitful.°C, Cosmology.Planck18, 1.0) ≈ -267.69899999999996 * Unitful.°C
true

julia> T_cmb(1.0, 2.7255) ≈ 5.451 * Unitful.K
true

julia> T_cmb(1.0, 2.7255 * Unitful.K) ≈ 5.451 * Unitful.K
true

julia> T_cmb(Unitful.°C, 1.0, 2.7255 * Unitful.K) ≈ -267.69899999999996 * Unitful.°C
true
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

# Examples
```jldoctest
julia> Cosmology.n_nu(Cosmology.Planck18)
3

julia> Cosmology.n_nu(3.046)
3
```
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

# Examples
```jldoctest
julia> T_nu(2.7255) ≈ 1.9453688391750839 * Unitful.K
true

julia> T_nu(Cosmology.Planck18) ≈ 1.9453688391750839 * Unitful.K
true

julia> T_nu(2.7255, 1) ≈ 3.8907376783501677 * Unitful.K
true

julia> T_nu(2.7255 * Unitful.K, 1) ≈ 3.8907376783501677 * Unitful.K
true

julia> T_nu(Cosmology.Planck18) ≈ 1.9453688391750839 * Unitful.K
true

julia> T_nu(Cosmology.Planck18, 1) ≈ 3.8907376783501677 * Unitful.K
true
```
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
 - `Neff`; effective number of neutrino species; see [`Cosmology.Neff`](@ref Cosmology.Neff).
 - `N_nu`; number of neutrino species; see [`Cosmology.n_nu`](@ref Cosmology.n_nu).
 - `nu_temp`; temperature of neutrino background in Kelvin; see [`T_nu`](@ref Cosmology.T_nu). This is the argument that carries the redshift dependence.
!!! note
    It is recommended that `length(m_nu) == N_nu`, unless `N_nu==0` in which case it doesn't matter. For example, if `N_nu==3` and you want one massive neutrino species with mass 0.06 eV, you should write `m_nu=(0.0,0.0,0.06)`. The current implementation is kind of stupid and can miscount the number of massless neutrinos if `length(m_nu) != N_nu`. """
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
Return the cosmological [`E`](@ref) factor times the square of the [`scale factor`](@ref Cosmology.scale_factor) `a`.

# Examples
```jldoctest
julia> Cosmology.a2E(Cosmology.Planck18,0.8)
0.7287593862161843
```
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
    scale_factor(z::Real)
    scale_factor(c::AbstractCosmology, z::Real)

Calculate the scale factor at redshift `z`. The scale factor is ``a=\\frac{1}{1+z}`` for cosmologies with FLRW metrics. The method that takes a cosmology is for compatibility with [`z_at_value`](@ref Cosmology.z_at_value). The redshift corresponding to a given scale factor can also be calculated with [`redshift`](@ref Cosmology.redshift). The derivative with respect to z is available as [`∇scale_factor`](@ref).

# Examples
```jldoctest
julia> scale_factor(1.0)
0.5

julia> scale_factor(Cosmology.Planck18, 1.0)
0.5
```
"""
scale_factor(z::Real) = 1 / (1 + z)
scale_factor(c::AbstractCosmology, z::Real) = scale_factor(z) #for compatibility with z_at_value
"""
    redshift(a::Real)

Calculate the redshift when the scale factor (see [`scale_factor`](@ref Cosmology.scale_factor)) is `a`. This is ``z = (1/a) - 1`` for cosmologies with FLRW metrics.
```jldoctest
julia> redshift(1.0)
0.0

julia> redshift(0.1)
9.0
```
"""
redshift(a::Real) = 1 / a - 1
""" 
    ∇scale_factor(z::Real)
    ∇scale_factor(c::AbstractCosmology, z::Real)
Calculate the derivative of the scale factor at redshift `z` with respect to `z` for cosmologies with FLRW metrics; ``\\frac{da}{dz} = -\\frac{1}{\\left(1+z\\right)^2}``.

# Examples
```jldoctest
julia> ∇scale_factor(1.0)
-0.25

julia> ∇scale_factor(Cosmology.Planck18, 1.0)
-0.25
```
"""
∇scale_factor(z::Real) = -1 / (1 + z)^2
∇scale_factor(c::AbstractCosmology, z) = ∇scale_factor(z)

"""
    E(c::AbstractCosmology, z::Real)
    E(z,h,OmegaM,OmegaK,OmegaL,OmegaG,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Return the Hubble parameter as a function of redshift, in units of H₀.
```math
E(z) \\equiv \\frac{H(z)}{H_0} = \\frac{H(z)}{h} \\ \\frac{1}{100 \\ \\text{km} \\ \\text{s}^{-1} \\ \\text{Mpc}^{-1}}
```
See also [`a2E`](@ref Cosmology.a2E), [`H`](@ref Cosmology.H).

# Examples
```jldoctest
julia> E(Cosmology.Planck18,1.0)
1.7828937335017068
```
"""
E(c::AbstractCosmology, z::Real) = (a = scale_factor(z); a2E(c, a) / a^2)
E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = (a=scale_factor(z); a2E(a,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) / a^2)
"""
    H(c::AbstractCosmology, z::Real)
Return the Hubble parameter as a function of redshift, in km / s / Mpc. From the Friedmann equation,

```math
H^2 = \\left( \\frac{\\dot{a}}{a} \\right) = \\frac{8 \\pi G \\rho}{3} - \\frac{kc^2}{a^2}
```

where `a` is the cosmological [`scale factor`](@ref Cosmology.scale_factor), ``\\dot{a}`` is the time derivative of `a`, `G` is Newton's gravitational constant, `c` is the speed of light, ``\\rho`` is a mass density, and `k` is the curvature parameter. `k` is typically -1, 0, or 1; `k=-1` corresponds to an open Universe (forever expanding), `k=0` corresponds to a "flat" Universe with [`ρ==ρ_c`](@ref Cosmology.ρ_c), and `k=1` corresponds to a closed Universe (eventual contraction). See also [`Ω_k`](@ref Cosmology.Ω_k).

# Examples
```jldoctest
julia> H(Cosmology.Planck18,1.0) ≈ 120.63059000872548 * Unitful.km / UnitfulAstro.Mpc / Unitful.s
true
```
"""
H(c::AbstractCosmology, z::Real) = 100 * h(c) * E(c, z) * u.km / u.s / ua.Mpc
"""
    hubble_dist0(c::AbstractCosmology)
Return the Hubble distance at present-day in Mpc, defined as the speed of light times the [`Hubble time`](@ref Cosmology.hubble_time0) at present-day.
```math
    D_0 = \\frac{c}{H_0}
```

# Examples
```jldoctest
julia> Cosmology.hubble_dist0(Cosmology.Planck18) ≈ 4430.866952409105 * UnitfulAstro.Mpc
true
```
"""
hubble_dist0(c::AbstractCosmology) = partype(c)(2997.92458) / h(c) * ua.Mpc # constant is speed of light in km/s divided by 100
"""
    hubble_dist(c::AbstractCosmology, z::Real)
Return the speed of light times the [`Hubble time`](@ref Cosmology.hubble_time) at redshift `z` in Mpc.
```math
    D(z) = \\frac{c}{H(z)}
```
# Examples
```jldoctest
julia> hubble_dist(Cosmology.Planck18,1.0) ≈ 2485.21090693758 * UnitfulAstro.Mpc
true
```
"""
hubble_dist(c::AbstractCosmology, z::Real) = hubble_dist0(c) / E(c, z)

"""
    hubble_time0(c::AbstractCosmology)
Return ``\\frac{1}{\\text{H}_0}`` in Gyr.

# Examples
```jldoctest
julia> Cosmology.hubble_time0(Cosmology.Planck18) ≈ 14.451555153425796 * Unitful.Gyr
true
```
"""
hubble_time0(c::AbstractCosmology) =  partype(c)(9.777922216807893) / h(c) * u.Gyr # 9.77814 # 9.777922216807893 / h(c) * u.Gyr
"""
    hubble_time(c::AbstractCosmology, z::Real)
Return ``\\frac{1}{\\text{H}\\left(z\\right)}`` in Gyr.

# Examples
```jldoctest
julia> hubble_time(Cosmology.Planck18, 1.0) ≈ 8.105673872689037 * Unitful.Gyr
true
```
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

# Examples
```jldoctest
julia> comoving_radial_dist(Cosmology.Planck18, 1.0) ≈ 3395.6344711515626 * UnitfulAstro.Mpc
true

julia> comoving_radial_dist(Cosmology.Planck18, 1.0, 2.0) ≈ 1912.5544127348157 * UnitfulAstro.Mpc
true

julia> comoving_radial_dist(UnitfulAstro.Gpc, Cosmology.Planck18, 1.0, 2.0) ≈ 1.9125544127348157 * UnitfulAstro.Gpc
true
```
"""
comoving_radial_dist(c::AbstractCosmology, z₁, z₂ = nothing; kws...) = hubble_dist0(c) * Z(c, z₁, z₂; kws...)

"""
    comoving_transverse_dist(c::AbstractCosmology, z₁, z₂ = nothing; kws...)

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

# Examples
```jldoctest
julia> comoving_transverse_dist(Cosmology.Planck18, 1.0) == comoving_radial_dist(Cosmology.Planck18, 1.0)
true

julia> comoving_transverse_dist(cosmology(OmegaK=0.1), 1.0) ≈ 3331.2531218753124 * UnitfulAstro.Mpc
true
```
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

```math
\\begin{aligned}
D_A(z_1,z_2) &= \\frac{D_T(z_1, z_2)}{1+z_2} = D_T(z_1, z_2) \\, a(z_2) \\newline
D_A(z) &= \\frac{D_T(z)}{1+z} = D_T(z) \\, a(z)
\\end{aligned}
```
where ``D_T`` is the [`comoving transverse distance`](@ref Cosmology.comoving_transverse_dist) and `a` is the [`scale factor`](@ref Cosmology.scale_factor). 

# Examples
```jldoctest
julia> angular_diameter_dist(Cosmology.Planck18, 1.0) ≈ 1697.8172355757813 * UnitfulAstro.Mpc
true

julia> angular_diameter_dist(Cosmology.Planck18, 1.0, 2.0) ≈ 637.5181375782719 * UnitfulAstro.Mpc
true

julia> angular_diameter_dist(UnitfulAstro.Gpc, Cosmology.Planck18, 1.0, 2.0) ≈ 0.6375181375782719 * UnitfulAstro.Gpc
true
```
"""
angular_diameter_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...) / (1 + z)
angular_diameter_dist(c::AbstractCosmology, z₁, z₂; kws...) =
    comoving_transverse_dist(c, z₁, z₂; kws...) / (1 + z₂)

"""
    luminosity_dist([u::Unitlike,] c::AbstractCosmology, z; kws...)

Bolometric luminosity distance in Mpc at redshift `z`. Will convert to compatible unit `u` if provided. `kws...` are integration options passed to `QuadGK.quadgk`.

```math
D_L(z) = D_T(z) \\times \\left(1 + z\\right) = \\frac{D_T(z)}{a\\left(z\\right)}
```
where ``D_T`` is the [`comoving transverse distance`](@ref Cosmology.comoving_transverse_dist) and `a` is the [`scale factor`](@ref Cosmology.scale_factor). 

# Examples
```jldoctest
julia> luminosity_dist(Cosmology.Planck18, 1.0) ≈ 6791.268942303125 * UnitfulAstro.Mpc
true

julia> luminosity_dist(UnitfulAstro.Gpc, Cosmology.Planck18, 1.0) ≈ 6.791268942303125 * UnitfulAstro.Gpc
true
```
"""
luminosity_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...) * (1 + z)

"""
    distmod(c::AbstractCosmology, z; kws...)

Distance modulus in magnitudes at redshift `z`. `kws...` are integration options passed to `QuadGK.quadgk`.

```math
\\mu(z) = 5 \\times \\text{log}_{10} \\left( D_L(z) \\right) + 25
```
where ``D_L(z)`` is the [`luminosity distance`](@ref Cosmology.luminosity_dist) in units of Mpc.

# Examples
```jldoctest
julia> distmod(Cosmology.Planck18,1.0)
44.159754646918806
```
"""
distmod(c::AbstractCosmology, z; kws...) = 5 * log10( u.ustrip(ua.Mpc, luminosity_dist(c, z; kws...))) + 25

#######################################################################################
# Volumes
#######################################################################################

"""
    comoving_volume([u::Unitlike,] c::AbstractCosmology, z; kws...)

Comoving volume in Gpc^3 at redshift `z`. Will convert to compatible unit `u` if provided. `kws...` are integration options passed to `QuadGK.quadgk`.

# Examples
```jldoctest
julia> comoving_volume(Cosmology.Planck18, 1.0) ≈ 164.00285577357855 * UnitfulAstro.Gpc^3
true

julia> comoving_volume(UnitfulAstro.Mpc^3, Cosmology.Planck18, 1.0) ≈ 1.6400285577357855e11 * UnitfulAstro.Mpc^3
true
```
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

Comoving volume element in Gpc per steradian out to redshift `z`. Will convert to compatible unit `u` if provided. `kws...` are integration options passed to `QuadGK.quadgk`.

# Examples
```jldoctest
julia> comoving_volume_element(Cosmology.Planck18, 1.0) ≈ 28.655310479576467 * UnitfulAstro.Gpc^3
true

julia> comoving_volume_element(UnitfulAstro.Mpc^3, Cosmology.Planck18, 1.0) ≈ 2.8655310479576466e10 * UnitfulAstro.Mpc^3
true
```
"""
comoving_volume_element(c::AbstractCosmology, z; kws...) =
    hubble_dist0(ua.Gpc, c) * angular_diameter_dist(ua.Gpc, c, z; kws...)^2 / a2E(c, scale_factor(z))

#############################################################################################
# Times
#############################################################################################

T(c::AbstractCosmology, a0, a1; kws...) = quadgk(x->x / a2E(c, x), a0, a1; kws...)[1]

"""
    age([u::Unitlike,] c::AbstractCosmology, z; kws...)

Return the age of the universe in Gyr at redshift `z`. Will convert to compatible unit `u` if provided. `kws...` are integration options passed to `QuadGK.quadgk`.

# Examples
```jldoctest
julia> age(Cosmology.Planck18, 0.0) ≈ 13.786885301987898 * Unitful.Gyr
true

julia> age(UnitfulAstro.Myr, Cosmology.Planck18, 0.0) ≈ 13786.885301987897 * Unitful.Myr
true
```
"""
age(c::AbstractCosmology, z; kws...) = hubble_time0(c) * T(c, 0, scale_factor(z); kws...)

"""
    lookback_time([u::Unitlike,] c::AbstractCosmology, z; kws...)

Difference between age at redshift 0 and age at redshift `z` in Gyr. Will convert to compatible unit `u` if provided. `kws...` are integration options passed to `QuadGK.quadgk`.

# Examples
```jldoctest
julia> lookback_time(Cosmology.Planck18, 1.0) ≈ 7.935542002084356 * Unitful.Gyr
true

julia> lookback_time(UnitfulAstro.Myr, Cosmology.Planck18, 1.0) ≈ 7935.542002084356 * Unitful.Myr
true
```
"""
lookback_time(c::AbstractCosmology, z; kws...) = hubble_time0(c) * T(c, scale_factor(z), 1; kws...)

#############################################################################################
# Dark energy
#############################################################################################

""" 
    w(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM,FlatLCDM,OpenLCDM,ClosedLCDM},z)
    w(z::Real,w0::Real,wa::Real)
Evaluates the redshift dependence of the dark energy equation of state,
```math
w \\equiv \\frac{p_\\Lambda}{\\rho_\\Lambda},
```
the ratio of the pressure to the energy density. The scaling factor, ``I(z)``, is defined by ``ρ_Λ(z) = I(z) \\ ρ_{Λ,0}``.
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
    ρ_c([u::UnitLike,] c::AbstractCosmology, z)
    ρ_c([u::UnitLike,] z, h, OmegaM, OmegaK, OmegaL, Tcmb0, m_nu, Neff, w0=-1, wa=0)
    ρ_c(h::Number,E::Number)
The critical density of the Universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.

```math
\\begin{aligned}
\\rho_{0,c} &\\equiv \\frac{3 H_0^2}{8 \\pi G} = 1.878 \\times 10^{-29} \\ h^2 \\ \\text{g/cm}^3 \\newline
\\rho_{c}(z) &\\equiv \\frac{3 H(z)^2}{8 \\pi G} = 1.878 \\times 10^{-29} \\ \\left(E(z) \\times h\\right)^2 \\ \\text{g/cm}^3
\\end{aligned}
```
where [`E`](@ref Cosmology.E) is the Hubble factor in units of ``H_0``. See, e.g., Equation 1.5.28 on page 57 of Weinberg's "Cosmology" for more information.

# Examples
```jldoctest
julia> ρ_c(Cosmology.Planck18, 0.0) ≈ 8.598814256619093e-30 * Unitful.g / Unitful.cm^3
true

julia> ρ_c(UnitfulAstro.Msun / UnitfulAstro.kpc^3, Cosmology.Planck18, 0.0) ≈ 127.05281539744222 * UnitfulAstro.Msun / UnitfulAstro.kpc^3
true
```

# Notes
If you want the critical density with the dimensionless hubble constant `h` factored out (i.e., in units of [density * h^2] you need to divide this result by `h^2`. 
"""
ρ_c(c::AbstractCosmology,z) = partype(c)(constants.RHO_C_Z0_CGS) * h(c)^2 * E(c,z)^2 * u.g/u.cm^3
ρ_c(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = constants.RHO_C_Z0_CGS * h^2 *
    E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2 * u.g/u.cm^3

"""
    ρ_m([u::UnitLike,] c::AbstractCosmology, z)
    ρ_m([u::UnitLike,] z, h, OmegaM)
The matter density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_m(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{m,0} \\ (1+z)^3 \\ \\text{g/cm}^3 
```

# Examples
```jldoctest
julia> ρ_m(Cosmology.Planck18, 0.0) ≈ 2.6627088227046682e-30 * Unitful.g / Unitful.cm^3
true

julia> ρ_m(UnitfulAstro.Msun / UnitfulAstro.kpc^3, Cosmology.Planck18, 0.0) ≈ 39.343174815971956 * UnitfulAstro.Msun / UnitfulAstro.kpc^3
true
```

# Notes
If you want the critical density with the dimensionless hubble constant `h` factored out (i.e., in units of [density * h^2] you need to divide this result by `h^2`. 
"""
ρ_m(c::AbstractCosmology,z) = partype(c)(constants.RHO_C_Z0_CGS) * h(c)^2 * Ω_m(c) * (1 + z)^3 * u.g/u.cm^3
ρ_m(z,h,Ω_m) = constants.RHO_C_Z0_CGS * h^2 * Ω_m * (1+z)^3 * u.g/u.cm^3

"""
    ρ_b([u::UnitLike,] c::AbstractCosmology, z)
    ρ_b([u::UnitLike,] z, h, OmegaB)
The baryon density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_b(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{b,0} \\ (1+z)^3 \\ \\text{g/cm}^3 
```

# Examples
```jldoctest
julia> ρ_b(Cosmology.Planck18, 0.0) ≈ 4.21083934146637e-31 * Unitful.g / Unitful.cm^3
true

julia> ρ_b(UnitfulAstro.Msun / UnitfulAstro.kpc^3, Cosmology.Planck18, 0.0) ≈ 6.221776370012746 * UnitfulAstro.Msun / UnitfulAstro.kpc^3
true
```

# Notes
If you want the critical density with the dimensionless hubble constant `h` factored out (i.e., in units of [density * h^2] you need to divide this result by `h^2`. 
"""
ρ_b(c::AbstractCosmology,z) = partype(c)(constants.RHO_C_Z0_CGS) * h(c)^2 * Ω_b(c) * (1 + z)^3 * u.g/u.cm^3
ρ_b(z,h,Ω_b) = constants.RHO_C_Z0_CGS * h^2 * Ω_b * (1+z)^3 * u.g/u.cm^3

"""
    ρ_dm([u::UnitLike,] c::AbstractCosmology, z)
    ρ_dm([u::UnitLike,] z, h, OmegaDM)
The dark matter density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_{dm}(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{dm,0} \\ (1+z)^3 \\ \\text{g/cm}^3 
```

# Examples
```jldoctest
julia> ρ_dm(Cosmology.Planck18, 0.0) ≈ 2.2416248885580312e-30 * Unitful.g / Unitful.cm^3
true

julia> ρ_dm(UnitfulAstro.Msun / UnitfulAstro.kpc^3, Cosmology.Planck18, 0.0) ≈ 33.12139844595921 * UnitfulAstro.Msun / UnitfulAstro.kpc^3
true
```

# Notes
If you want the critical density with the dimensionless hubble constant `h` factored out (i.e., in units of [density * h^2] you need to divide this result by `h^2`. 
"""
ρ_dm(c::AbstractCosmology,z) = partype(c)(constants.RHO_C_Z0_CGS) * h(c)^2 * Ω_dm(c) * (1 + z)^3 * u.g/u.cm^3
ρ_dm(z,h,Ω_dm) = constants.RHO_C_Z0_CGS * h^2 * Ω_dm * (1+z)^3 * u.g/u.cm^3

"""
    ρ_Λ([u::UnitLike,] c::AbstractCosmology, z)
    ρ_Λ([u::UnitLike,] z, h, Ω_Λ, w0=-1, wa=0)
The dark energy density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided. For a general dark energy equation of state ``w(z)``,
```math
\\rho_\\Lambda(z) = \\rho_{\\Lambda,0} \\ \\exp \\left[  3 \\int_0^z  \\frac{1+w \\left( z^\\prime\\right)}{1+z^\\prime} \\ dz^\\prime \\ \\right]
```
For constant ``w``, this reduces to
```math
\\rho_\\Lambda(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{\\Lambda,0} \\ \\left(1+z\\right)^{3\\left(1+w\\right)} \\ \\text{g/cm}^3 
```
and for a cosmological constant ``w=-1``,
```math
\\rho_\\Lambda(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{\\Lambda,0} \\ \\text{g/cm}^3 
```

See also [`de_density_scale`](@ref Cosmology.de_density_scale).

# Examples
```jldoctest
julia> ρ_Λ(Cosmology.Planck18, 0.0) ≈ 5.923261432735804e-30 * Unitful.g / Unitful.cm^3
true

julia> ρ_Λ(UnitfulAstro.Msun / UnitfulAstro.kpc^3, Cosmology.Planck18, 0.0) ≈ 87.51986249556084 * UnitfulAstro.Msun / UnitfulAstro.kpc^3
true

julia> ρ_Λ(Cosmology.Planck18, 0.0) == ρ_Λ(Cosmology.Planck18, 1.0)
true
```

# Notes
If you want the critical density with the dimensionless hubble constant `h` factored out (i.e., in units of [density * h^2] you need to divide this result by `h^2`. 
"""
ρ_Λ(c::AbstractCosmology,z) = partype(c)(constants.RHO_C_Z0_CGS) * h(c)^2 * Ω_Λ(c) / de_density_scale(c,z) * u.g/u.cm^3
ρ_Λ(z,h,Ω_Λ,w0=-1,wa=0) = constants.RHO_C_Z0_CGS * h^2 * Ω_Λ / de_density_scale(z,w0,wa) * u.g/u.cm^3

"""
    ρ_γ([u::UnitLike,] c::AbstractCosmology, z)
    ρ_γ([u::UnitLike,] z, h, Ω_γ)
The photon matter density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_\\gamma(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{\\gamma,0} \\ (1+z)^4 \\ \\text{g/cm}^3 
```

# Examples
```jldoctest
julia> ρ_γ(Cosmology.Planck18, 0.0) ≈ 4.645092477570597e-34 * Unitful.g / Unitful.cm^3
true

julia> ρ_γ(UnitfulAstro.Msun / UnitfulAstro.kpc^3, Cosmology.Planck18, 0.0) ≈ 0.006863412319931541 * UnitfulAstro.Msun / UnitfulAstro.kpc^3
true
```

# Notes
If you want the critical density with the dimensionless hubble constant `h` factored out (i.e., in units of [density * h^2] you need to divide this result by `h^2`. 
"""
ρ_γ(c::AbstractCosmology,z) = partype(c)(constants.RHO_C_Z0_CGS) * h(c)^2 * Ω_γ(c) * (1 + z)^4 * u.g/u.cm^3
ρ_γ(z,h,Ω_γ) = constants.RHO_C_Z0_CGS * h^2 * Ω_γ * (1+z)^4 * u.g/u.cm^3

"""
    ρ_ν([u::UnitLike,] c::AbstractCosmology, z)
    ρ_ν([u::UnitLike,] z, h, Tcmb0, Neff, m_nu, N_nu=Cosmology.n_nu(Neff))
The neutrino energy density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_\\nu(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{\\nu,0} \\ (1+z)^4 \\ \\text{g/cm}^3 
```

# Examples
```jldoctest
julia> ρ_ν(Cosmology.Planck18, 0.0) ≈ 1.2379491930863195e-32 * Unitful.g / Unitful.cm^3
true

julia> ρ_ν(UnitfulAstro.Msun / UnitfulAstro.kpc^3, Cosmology.Planck18, 0.0) ≈ 0.1829146735894845 * UnitfulAstro.Msun / UnitfulAstro.kpc^3
true
```

# Notes
If you want the critical density with the dimensionless hubble constant `h` factored out (i.e., in units of [density * h^2] you need to divide this result by `h^2`. 
"""
ρ_ν(c::AbstractCosmology,z) = ρ_γ(c,z) * nu_relative_density(c, z)
ρ_ν(z,h,Tcmb0,Neff,m_nu,N_nu=n_nu(Neff)) = (Ω_γ = 4.481620089297254e-7 * Tcmb0^4 / h^2; ρ_γ(z,h,Ω_γ) * nu_relative_density(m_nu, Neff, u.ustrip(T_nu(Tcmb0, z)), N_nu) )

"""
    ρ_r([u::UnitLike,] c::AbstractCosmology, z)
    ρ_r([u::UnitLike,] z, h, Tcmb0, Neff, m_nu, N_nu=Cosmology.n_nu(Neff))
The energy density of the universe in relativistic species at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
```math
\\rho_r(z) = 1.878 \\times 10^{-29} \\ h^2 \\ \\Omega_{r,0} \\ (1+z)^4 \\ \\text{g/cm}^3 
```

# Examples
```jldoctest
julia> ρ_r(Cosmology.Planck18, 0.0) ≈ 1.2844001178620257e-32 * Unitful.g / Unitful.cm^3
true

julia> ρ_r(UnitfulAstro.Msun / UnitfulAstro.kpc^3, Cosmology.Planck18, 0.0) ≈ 0.18977808590941606 * UnitfulAstro.Msun / UnitfulAstro.kpc^3
true

julia> ρ_r(Cosmology.Planck18, 0.0) ≈ ρ_γ(Cosmology.Planck18, 0.0) + ρ_ν(Cosmology.Planck18, 0.0)
true
```

# Notes
If you want the critical density with the dimensionless hubble constant `h` factored out (i.e., in units of [density * h^2] you need to divide this result by `h^2`. 
"""
ρ_r(c::AbstractCosmology,z) = ρ_γ(c,z) * (1 + nu_relative_density(c, z))
ρ_r(z,h,Tcmb0,Neff,m_nu,N_nu=n_nu(Neff)) = (Ω_γ = 4.481620089297254e-7 * Tcmb0^4 / h^2; ρ_γ(z,h,Ω_γ) * (1.0 + nu_relative_density(m_nu, Neff, u.ustrip(T_nu(Tcmb0, z)), N_nu) ) )

#############################################################################################
# Quantities Involving Densities
#############################################################################################
"""
    lagrangianR([u::UnitLike,] M::Union{Real,u.Mass}, c::AbstractCosmology, z::Real=0.0)
The Lagrangian radius (in Mpc) of a sphere of mass `M` (in solar masses) at redshift `z`; i.e. the radius of a sphere that would enclose the mass `M` given the mean density of the universe at redshift `z`. 

```math
R(z) = \\left( \\frac{3 M}{4π \\ ρ_m(z)} \\right)^{1/3}
```

# Examples
```jldoctest
julia> lagrangianR(10^12, Cosmology.Planck18, 0.0) == lagrangianR(10^12, Cosmology.Planck18) == lagrangianR(10^12 * UnitfulAstro.Msun, Cosmology.Planck18) ≈ 1.8239544820629736 * UnitfulAstro.Mpc
true

julia> lagrangianR(Unitful.m, 10^12, Cosmology.Planck18, 0.0) ≈ 5.628135454962416e22 * Unitful.m
true
```

# Notes
If you want to provide `M` in units of `Msun / h` and get out `R` in units of `Mpc / h`, then you need to multiply the result by `h^(2/3)`; i.e. `R[Mpc/h] = lagrangianR(M[Msun/h], c, z) * h^(2/3)`.
"""
lagrangianR(M::Real, c::AbstractCosmology, z::Real=0.0) = cbrt(3 * M / (4π * u.ustrip(ua.Msun/ua.Mpc^3, ρ_m(c,z)))) * ua.Mpc
lagrangianR(M::u.Mass, c::AbstractCosmology, z::Real=0.0) = lagrangianR(u.ustrip(ua.Msun,M), c, z)

"""
    lagrangianM([u::UnitLike,] R::Union{Number,u.Length}, c::AbstractCosmology, z::Number=0.0)
The Lagrangian mass of a sphere of radius `R` in Mpc at redshift `z`; i.e. the mass enclosed by a sphere of radius `R` at redshift `z` given the mean density of the universe at redshift `z`.

```math
M(z) = \\frac{4π}{3} R^3 ρ_m(z)
```

# Examples
```jldoctest
julia> lagrangianM(8.0, Cosmology.Planck18, 0.0) == lagrangianM(8.0, Cosmology.Planck18) == lagrangianM(8.0 * UnitfulAstro.Mpc, Cosmology.Planck18) ≈ 8.437775631070308e13 * UnitfulAstro.Msun
true

julia> lagrangianM(Unitful.kg, 8.0, Cosmology.Planck18, 0.0) ≈ 1.6777756351555676e44 * Unitful.kg
true
```

# Notes
If you want to provide `R` in units of `Mpc / h` and get out `M` in units of `Msun / h`, then you need to divide the result by `h^2`; i.e. `M[Msun/h] = lagrangianR(R[Mpc/h], c, z) / h^2`.
"""
lagrangianM(R::Number, c::AbstractCosmology, z::Number=0.0) = 4π * R^3 * u.ustrip(ua.Msun/ua.Mpc^3, ρ_m(c,z)) / 3 * ua.Msun
lagrangianM(R::u.Length, c::AbstractCosmology, z::Number=0.0) = lagrangianM(u.ustrip(ua.Mpc,R), c, z)

#############################################################################################
# Omegas
#############################################################################################

""" 
    Ω_m(c::AbstractCosmology,z)
    Ω_m(c::AbstractCosmology)
    Ω_m(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Density of matter at redshift `z` in units of the critical density. When called without a redshift, returns the `z=0` value.

# Examples
```jldoctest
julia> Ω_m(Cosmology.Planck18)
0.30966

julia> Ω_m(Cosmology.Planck18, 1.0)
0.7793349973337195
```
"""
Ω_m(c::AbstractCosmology) = c.Ω_m
Ω_m(c::AbstractCosmology,z) = Ω_m(c) * (1 + z)^3 / E(c,z)^2
Ω_m(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = OmegaM * (1+z)^3 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2
""" 
    Ω_b(c::AbstractCosmology,z)
    Ω_b(c::AbstractCosmology)
    Ω_b(z,h,OmegaM,OmegaB,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Density of baryons at redshift `z` in units of the critical density. When called without a redshift, returns the `z=0` value.

# Examples
```jldoctest
julia> Ω_b(Cosmology.Planck18)
0.04897

julia> Ω_b(Cosmology.Planck18, 1.0)
0.1232449616335085
```
"""
Ω_b(c::AbstractCosmology) = c.Ω_b
Ω_b(c::AbstractCosmology,z) = Ω_b(c) * (1 + z)^3 / E(c,z)^2
Ω_b(z,h,OmegaM,OmegaB,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = OmegaB * (1+z)^3 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2
""" 
    Ω_dm(c::AbstractCosmology,z)
    Ω_dm(c::AbstractCosmology)
    Ω_dm(z,h,OmegaM,OmegaB,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Density of dark matter at redshift `z` in units of the critical density. When called without a redshift, returns the `z=0` value.

# Examples
```jldoctest
julia> Ω_dm(Cosmology.Planck18)
0.26069

julia> Ω_dm(Cosmology.Planck18, 1.0)
0.6560900357002108
```
"""
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
where `k` is the curvature parameter in the Friedmann equation. See Equations  1.5.19 and 1.5.40 on pages 56 and 60 in Weinberg's "Cosmology".

# Examples
```jldoctest
julia> Ω_k(Cosmology.Planck18) == Ω_k(Cosmology.Planck18, 1.0) == 0.0
true

julia> Ω_k( cosmology(OmegaK = 0.1 ), 1.0)
0.11498515039500401

julia> Ω_k( cosmology(OmegaK = -0.1 ), 1.0)
-0.13895112427920248
```
"""
Ω_k(c::AbstractCosmology) = c.Ω_k
Ω_k(c::AbstractCosmology,z) = Ω_k(c) * (1 + z)^2 / E(c,z)^2
Ω_k(c::AbstractFlatCosmology,z=0.0) = zero(z)
Ω_k(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = OmegaK * (1+z)^2 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2
""" 
    Ω_γ(c::AbstractCosmology,z)
    Ω_γ(c::AbstractCosmology)
    Ω_γ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density of photons at redshift `z` in units of the critical density. Calculated from [`T_cmb`](@ref). When called without a redshift, returns the `z=0` value.

# Examples
```jldoctest
julia> Ω_γ(Cosmology.Planck18) ≈ 5.4020151371393475e-5
true

julia> Ω_γ(Cosmology.Planck18, 10000.0) ≈ 0.44150560915009124
true
```
"""
Ω_γ(c::AbstractCosmology) = partype(c)(4.481620089297254e-7) * u.ustrip(u.K,T_cmb(c))^4 / h(c)^2
Ω_γ(c::AbstractCosmology,z) = Ω_γ(c) * (1 + z)^4 / E(c,z)^2
Ω_γ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = (OmegaG = 4.481620089297254e-7 * Tcmb0^4 / h^2; OmegaG * (1+z)^4 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2 )
""" 
    Ω_ν(c::AbstractCosmology,z)
    Ω_ν(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density in neutrinos at redshift `z` in units of the critical density. When called without a redshift, returns the `z=0` value.

# Examples
```jldoctest
julia> Ω_ν(Cosmology.Planck18) ≈ 0.0014396743040860382
true

julia> Ω_ν(Cosmology.Planck18, 10000.0) ≈ 0.30543520244776484
true
```
"""
Ω_ν(c::AbstractCosmology) = Ω_γ(c) * nu_relative_density(c)
Ω_ν(c::AbstractCosmology,z) = Ω_γ(c,z) * nu_relative_density(c,z)
Ω_ν(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = (nu_temp = u.ustrip(u.K,T_nu(Tcmb0, z)); Ω_γ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) * nu_relative_density(m_nu, Neff, nu_temp) )
""" 
    Ω_r(c::AbstractCosmology,z)
    Ω_r(c::AbstractCosmology)
    Ω_r(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density in radiation at redshift `z` in units of the critical density. Evaluated as `Ω_γ(c,z) + Ω_ν(c,z)`; sum of photons and neutrinos. When called without a redshift, returns the `z=0` value.

# Examples
```jldoctest
julia> Ω_r(Cosmology.Planck18) ≈ 0.0014936944554574316
true

julia> Ω_r(Cosmology.Planck18, 10000.0) ≈ 0.7469408115978561
true

julia> Ω_r(Cosmology.Planck18, 10000.0) == Ω_γ(Cosmology.Planck18, 10000.0) + Ω_ν(Cosmology.Planck18, 10000.0)
true
```
"""
Ω_r(c::AbstractCosmology) = Ω_γ(c) + Ω_ν(c)
Ω_r(c::AbstractCosmology,z) = Ω_γ(c,z) + Ω_ν(c,z)
Ω_r(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = Ω_γ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) + Ω_ν(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)
""" 
    Ω_Λ(c::AbstractCosmology,z)
    Ω_Λ(c::AbstractCosmology)
    Ω_Λ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density in dark energy at redshift `z` in units of the critical density. When called without a redshift, returns the `z=0` value.

# Examples
```jldoctest
julia> Ω_Λ(Cosmology.Planck18)
0.6888463055445425

julia> Ω_Λ(Cosmology.Planck18, 1.0)
0.21670623978512665
```
"""
Ω_Λ(c::AbstractCosmology) = c.Ω_Λ
Ω_Λ(c::AbstractCosmology,z) = Ω_Λ(c) * de_density_scale(c,z) / E(c,z)^2
Ω_Λ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = OmegaL * de_density_scale(z,w0,wa) / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2

"""
    δc(c::AbstractCosmology, z::Real)
    δc(z::Real, h, OmegaM, OmegaK, OmegaL, Tcmb0, m_nu, Neff, w0=-1, wa=0)
The linear overdensity threshold for halo collapse as calculated from the spherical top-hat collapse model. The canonical value for an Einstein-de-Sitter (EdS) cosmology is ``3/5(3π/2)^{2/3} ≈ 1.686`` but there is a redshift dependence for non-EdS cosmologies. We adopt the fitting function from [Nakamura & Suto 1997](https://ui.adsabs.harvard.edu/abs/1997PThPh..97...49N/abstract) (Equation C-28 in the ArXiv version), also given by Equation 56 in [Chisari et al. 2019](https://ui.adsabs.harvard.edu/abs/2019ApJS..242....2C/abstract).

# Examples
```jldoctest
julia> δc(Cosmology.Planck18, 0.0) ≈ 1.675910191226453
true

julia> δc(Cosmology.Planck18, 10.0) ≈ 1.686388790734125
true

julia> δc(Cosmology.Planck18,10^10,0.0,1.0,1.5) ≈ 1.712312778883257
true
```
"""
δc(c::AbstractCosmology, z::Real) = 3 * (12π)^(2/3) / 20 * (1 + 0.012299 * log10(Ω_m(c,z)))
δc(z::Real,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = 3 * (12π)^(2/3) / 20 * (1 + 0.012299 * log10(Ω_m(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)))
""" 
    δc(c::AbstractCosmology, M::Union{Real,AbstractArray}, z::Real, mx::Real, gx::Real=1.5)
Return the linear overdensity threshold for collapse in warm dark matter cosmologies as formulated by Equations 7--10 in [Benson 2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.428.1774B/abstract). `M` is the halo mass, `z` is the redshift of evaluation, `mx` is the mass of the WDM particle in keV, and `gx` is the effective degrees of freedom with 1.5 being the expected value for a fermionic spin-1/2 particle.

# Examples
```jldoctest
julia> δc(Cosmology.Planck18, 10^10, 0.0, 1.0, 1.5) ≈ 1.712312778883257
true
```
"""
function δc(c::AbstractCosmology, M::Union{Real,AbstractArray}, z::Real, mx::Real, gx::Real=1.5)
    lcdm_val = δc(c,z)
    OmegaM = Ω_m(c)
    littleh = h(c)
    zeq = 3600 * (OmegaM * littleh^2 / 0.15) - 1
    # jeans_mass = 3.06e8 * ((1+zeq)/3000)^1.5 * sqrt(OmegaM * littleh^2 / 0.15) * (1.5/gx) * (1/mx)^4
    jeans_mass = 3.06e8 * sqrt(((1+zeq)/3000)^3) * sqrt(OmegaM * littleh^2 / 0.15) * (1.5/gx) * (1/mx)^4
    
    x = @. log(M/jeans_mass)
    hh = @. (1+exp((x+2.4)/0.1))^-1
    return @. lcdm_val * ((hh * 0.04 / exp(2.3*x)) + ((1-hh) * exp(0.31687/exp(0.809*x))))
end

#############################################################################################
# Roots
#############################################################################################

"""
    z_at_value(c::AbstractCosmology, func::Function, fval; zmin=1e-8, zmax=1000.0, kws...)

Find the redshift `z` at which `func(c,z) == fval` for cosmology instance `c`. This uses a numerical root finder and searches between `zmin` and `zmax`. Additional `kws...` are passed through to [`Roots.find_zero`](https://juliamath.github.io/Roots.jl/stable/reference/#Roots.find_zero).

# Examples
```julia
julia> z_at_value(Cosmology.Planck18, scale_factor, 0.8) ≈ 0.25
true

julia> z_at_value(Cosmology.Planck18, distmod, 50.0) ≈ 9.5f0
true
```
!!! warning
    Not all cosmological methods defined in this module are monotonic with redshift, such that there may not be a unique solution over the default `zmin->zmax` (e.g., [`angular_diameter_dist`](@ref Cosmology.angular_diameter_dist)). In this case, you need to make sure that `zmin` and `zmax` properly bracket the solution you are looking for.
```julia
julia> z_at_value(Cosmology.Planck18, angular_diameter_dist, 1250.0; zmin=1e-5, zmax=2.5)
0.46668775101654764

julia> z_at_value(Cosmology.Planck18, angular_diameter_dist, 1250.0; zmin=2.5, zmax=10.0)
5.595635655898187
```
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

    result = find_zero(f, (zmin,zmax), Bisection(); kws...) # atol=ztol, maxevals=maxfun)
    return result
end

#########################################################################################
# Misc
#########################################################################################

"""
    sound_horizon(c::AbstractCosmology)
Return the sound horizon length (in Mpc), given by Equation 26 in Eisenstein & Hu 1998.

```math
s = \\frac{44.5 \\, \\text{log} \\left( 9.83 / \\Omega_0 / h^2 \\right)}{ \\sqrt{ 1 + 10 \\, \\left( \\Omega_b \\, h^2 \\right)^{3/4} } } \\ \\text{Mpc}
```

# Examples
```jldoctest
julia> sound_horizon(Cosmology.Planck18) ≈ 150.10339082768203 * UnitfulAstro.Mpc
true
```
"""
sound_horizon(c::AbstractCosmology) = (h2 = h(c)^2; 44.5 * log(9.83 / Ω_m(c) / h2) / sqrt(1.0 + 10.0 * (Ω_b(c) * h2)^0.75) * ua.Mpc)

"""
    matter_radiation_equality(c::AbstractCosmology)
Return the redshift of matter-radiation equality. This is computed as Equation 2 in Eisenstein and Hu 1998. I previously had a different formula but couldn't figure out where it came from.

# Examples
```jldoctest
julia> matter_radiation_equality(Cosmology.Planck18) ≈ 3413.1817608491015
true
```
"""
matter_radiation_equality(c::AbstractCosmology) = 2.5e4 * Ω_m(c) * h(c)^2 / (u.ustrip(u.K,T_cmb(c))/2.7)^4
# cant figure out where this came from
# matter_radiation_equality(c::AbstractCosmology) = 3600 * (Ω_m(c) * h(c)^2 / 0.15) - 1

#########################################################################################
# Easily select a different unit
#############################################################################################

for f in (:hubble_dist0, :hubble_dist, :hubble_time0, :hubble_time, :comoving_radial_dist,
          :comoving_transverse_dist, :angular_diameter_dist, :luminosity_dist,
          :comoving_volume, :comoving_volume_element, :age, :lookback_time,
          :T_nu, :T_cmb, :ρ_c, :ρ_m, :ρ_b, :ρ_dm, :ρ_Λ, :ρ_γ, :ρ_ν, :ρ_r, :lagrangianR, :lagrangianM, :sound_horizon)
    @eval $f(uu::u.Unitlike, args...; kws...) = u.uconvert(uu, $f(args...; kws...))
end
