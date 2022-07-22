# import QuadGK: quadgk
# import Unitful as u
# import UnitfulAstro as ua
# import Roots: find_zero, Bisection

# include("types_base.jl")

#############################################################################################
# neutrinos
"""
    n_nu(c::AbstractCosmology)
    n_nu(Neff::Real)
Returns the number of discrete neutrinos in the cosmology.
"""
n_nu(c::AbstractCosmology) = Int(floor(c.Neff)) #number of neutrinos
n_nu(Neff::Real) = Int(floor(Neff))
"""
    T_nu([u::Unitlike,], c::AbstractCosmology, z::Real)
    T_nu(Tcmb0::Real, z::Real)
    T_nu(Tcmb0::u.Quantity, z::Real)
The neutrino temperature of the cosmology at redshift z in Kelvin. Will convert to compatible unit `u` if provided.
"""
T_nu(c::AbstractCosmology, z::Real) = 0.7137658555036082 * c.Tcmb0 * (1.0 + z) * u.K #neutrino temperature
T_nu(Tcmb0::Real, z::Real) = 0.7137658555036082 * Tcmb0 * (1.0 + z) * u.K #neutrino temperature
T_nu(Tcmb0::u.Quantity, z::Real) = 0.7137658555036082 * (Tcmb0 |> u.K) * (1.0 + z) #neutrino temperature

""" 
    nu_relative_density(m_nu::Real, Neff, N_nu, nu_temp)
    nu_relative_density(m_nu::AbstractArray, Neff, N_nu, nu_temp)
    nu_relative_density(c::AbstractCosmology, z)
Neutrino density function relative to the energy density in photons. """
@inline function nu_relative_density(m_nu::Real, Neff, nu_temp, N_nu=nothing)
    N_nu === nothing && (N_nu = n_nu(Neff))
    prefac = 0.22710731766023898 #7/8 * (4/11)^(4/3)
    m_nu==0 && return prefac * Neff
    p = 1.83
    invp = 0.54644808743  # 1.0 / p
    k = 0.3173
    nu_y = m_nu / (8.617333262145179e-5 * nu_temp )
    rel_mass = (1.0 + (k * nu_y)^p)^invp
    return prefac * m_nu * Neff / N_nu
end
@inline function nu_relative_density(m_nu::AbstractArray, Neff, nu_temp, N_nu=nothing)
    N_nu === nothing && (N_nu = n_nu(Neff))
    prefac = 0.22710731766023898 #7/8 * (4/11)^(4/3)
    p = 1.83
    invp = 0.54644808743  # 1.0 / p
    k = 0.3173
    #rel_mass_per=zeros(eltype(m_nu), length(m_nu)) #allows use of staticarrays for m_nu
    rel_mass_per_sum=0.0
    massless=0
    for i in eachindex(m_nu)
        if m_nu[i]==0
            massless+=1
        else
            nu_y = m_nu[i] / (8.617333262145179e-5 * nu_temp )
            # rel_mass_per[i] = (1.0 + (k * nu_y)^p)^invp
            rel_mass_per_sum += (1.0 + (k * nu_y)^p)^invp
        end
    end
    # return prefac * ( sum(rel_mass_per) + massless ) * Neff / N_nu
    return prefac * ( rel_mass_per_sum + massless ) * Neff / N_nu
end

@inline function nu_relative_density(c::AbstractCosmology, z)
    (c.Neff == 0 || c.Neff === nothing) && return 0

    # this is more concise but slower due to the broadcast operations, which are slower than coded loops
    # prefac = 0.22710731766023898 #7/8 * (4/11)^(4/3)
    # if 0 in (c.m_nu .== 0) #if at least one neutrino has mass
    #     p = 1.83
    #     invp = 0.54644808743  # 1.0 / p
    #     k = 0.3173
    #     nu_y = c.m_nu ./ (8.617333262145179e-5 * (T_nu(c,z) |> u.ustrip))  #the constant is k_B in ev / K
    #     rel_mass_per = @. (1 + (k * nu_y[nu_y!=0])^p)^invp
    #     rel_mass = sum(rel_mass_per) + length(c.m_nu[c.m_nu.==0])
    #     return prefac * rel_mass * (c.Neff / n_nu(c) )
        
    # else #if all neutrinos are massless
    #     return prefac * c.Neff
    # end

    if isscalar(c.m_nu)
        # return _nu_relative_density_scalar(c,z)
        return nu_relative_density(c.m_nu, c.Neff, T_nu(c,z)|>u.ustrip)
    else
        # return _nu_relative_density_array(c,z)
        return nu_relative_density(c.m_nu, c.Neff, T_nu(c,z)|>u.ustrip)

    end
end

##############################################################################
# a2E(c::FlatLCDM, a) = sqrt(c.Ω_r + c.Ω_m * a + c.Ω_Λ * a^4)
"""
    a2E(c::AbstractCosmology,a)
    a2E(a,OmegaM,OmegaK,OmegaL,OmegaG,Tcmb0,m_nu,Neff,w0=-1,wa=0)
The cosmological `E` factor times the scale factor `a` squared. See also the `E` method. 
"""
function a2E(c::FlatLCDM, a)
    z = 1/a-1
    Or = c.Ω_γ * (1+nu_relative_density(c,z))
    return sqrt(Or + c.Ω_m * a + c.Ω_Λ * a^4)
end

function a2E(c::Union{ClosedLCDM,OpenLCDM}, a)
    a2 = a * a
    z = 1/a-1
    Or = c.Ω_γ * (1+nu_relative_density(c,z))
    sqrt(Or + c.Ω_m * a + (c.Ω_k + c.Ω_Λ * a2) * a2)
end

function a2E(c::Union{FlatWCDM,ClosedWCDM,OpenWCDM}, a)
    z = 1/a-1
    Or = c.Ω_γ * (1+nu_relative_density(c,z))
    ade = exp((1 - 3 * (c.w0 + c.wa)) * log(a) + 3 * c.wa * (a - 1))
    sqrt(Or + (c.Ω_m + c.Ω_k * a) * a + c.Ω_Λ * ade)
end

function a2E(a,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
    z = 1/a-1
    OmegaG = 4.481620089297254e-7 * Tcmb0^4 / h^2
    nu_temp = T_nu(Tcmb0, z) |> u.ustrip
    Or = OmegaG * (1+nu_relative_density(m_nu, Neff, nu_temp))
    ade = exp((1 - 3 * (w0 + wa)) * log(a) + 3 * wa * (a - 1))
    sqrt(Or + (OmegaM + OmegaK * a) * a + OmegaL * ade)
end

###############################################################################
# hubble rate
""" 
    scale_factor(z::Real) or scale_factor(c::AbstractCosmology, z::Real)

Scale factor at redshift ``z``. The scale factor is defined as :math:`a = 1 / (1 + z)`.
Derivative with respect to z is available as ∇scale_factor(z)
"""
scale_factor(z::Real) = 1 / (1 + z)
scale_factor(c::AbstractCosmology, z::Real) = scale_factor(z) #for compatibility with z_at_value
""" 
    ∇scale_factor(z::Real) or ∇scale_factor(c::AbstractCosmology, z::Real)
Derivative of the scale factor at redshift ``z``. """
∇scale_factor(z::Real) = -1 / (1 + z)^2
∇scale_factor(c::AbstractCosmology,z) = ∇scale_factor(z)

"""
    E(c::AbstractCosmology, z::Real)
    E(z,h,OmegaM,OmegaK,OmegaL,OmegaG,Tcmb0,m_nu,Neff,w0=-1,wa=0)
The cosmological `E` factor used in the density evolutions and other quantities. See also a2E. 
"""
E(c::AbstractCosmology, z::Real) = (a = scale_factor(z); a2E(c, a) / a^2)
E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = (a=scale_factor(z); a2E(a,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) / a^2)
H(c::AbstractCosmology, z::Real) = 100 * c.h * E(c, z) * u.km / u.s / ua.Mpc

hubble_dist0(c::AbstractCosmology) = 2997.92458 / c.h * ua.Mpc
hubble_dist(c::AbstractCosmology, z::Real) = hubble_dist0(c) / E(c, z)

hubble_time0(c::AbstractCosmology) =  9.777922216807893 / c.h * u.Gyr # 9.77814
hubble_time(c::AbstractCosmology, z::Real) = hubble_time0(c) / E(c, z)

###############################################################################
# distances

Z(c::AbstractCosmology, z::Real, ::Nothing; kws...) =
    quadgk(a->1 / a2E(c, a), scale_factor(z), 1; kws...)[1]
Z(c::AbstractCosmology, z₁::Real, z₂::Real; kws...) =
    quadgk(a->1 / a2E(c, a), scale_factor(z₂), scale_factor(z₁); kws...)[1]


"""
    comoving_radial_dist([u::Unitlike,] c::AbstractCosmology, [z₁,] z₂; kws...)

Comoving radial distance in Mpc at redshift `z₂` as seen by an observer at `z₁`.  Redshift `z₁` defaults to 0 if omitted.  Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
comoving_radial_dist(c::AbstractCosmology, z₁, z₂ = nothing; kws...) = hubble_dist0(c) * Z(c, z₁, z₂; kws...)

"""
    comoving_transverse_dist(c::T<:AbstractCosmology, z₁, z₂ = nothing; kws...)

Comoving transverse distance at a given redshift, or between two redshifts.

This value is the transverse comoving distance at redshift ``z`` corresponding to an angular separation of 1 radian. This is the same as the comoving distance if omega_k is zero (as in the current concordance lambda CDM model). Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
comoving_transverse_dist(c::AbstractFlatCosmology, z₁, z₂ = nothing; kws...) =
    comoving_radial_dist(c, z₁, z₂; kws...)
function comoving_transverse_dist(c::AbstractOpenCosmology, z₁, z₂ = nothing; kws...)
    sqrtok = sqrt(c.Ω_k)
    hubble_dist0(c) * sinh(sqrtok * Z(c, z₁, z₂; kws...)) / sqrtok
end
function comoving_transverse_dist(c::AbstractClosedCosmology, z₁, z₂ = nothing; kws...)
    sqrtok = sqrt(abs(c.Ω_k))
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
# volumes

"""
    comoving_volume([u::Unitlike,] c::AbstractCosmology, z; kws...)

Comoving volume in cubic Gpc out to redshift `z`. Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
comoving_volume(c::AbstractFlatCosmology, z; kws...) =
    (4pi / 3) * (comoving_radial_dist(ua.Gpc, c, z; kws...))^3
function comoving_volume(c::AbstractOpenCosmology, z; kws...)
    DH = hubble_dist0(ua.Gpc, c)
    x = comoving_transverse_dist(ua.Gpc, c, z; kws...) / DH
    sqrtok = sqrt(c.Ω_k)
    2π * (DH)^3 * (x * sqrt(1 + c.Ω_k * x^2) - asinh(sqrtok * x) / sqrtok) / c.Ω_k
end
function comoving_volume(c::AbstractClosedCosmology, z; kws...)
    DH = hubble_dist0(ua.Gpc, c)
    x = comoving_transverse_dist(ua.Gpc, c, z; kws...) / DH
    sqrtok = sqrt(abs(c.Ω_k))
    2pi * (DH)^3 * (x * sqrt(1 + c.Ω_k * x^2) - asin(sqrtok * x) / sqrtok) / c.Ω_k
end

"""
    comoving_volume_element([u::Unitlike,] c::AbstractCosmology, z; kws...)

Comoving volume element in Gpc out to redshift `z`. Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
comoving_volume_element(c::AbstractCosmology, z; kws...) =
    hubble_dist0(ua.Gpc, c) * angular_diameter_dist(ua.Gpc, c, z; kws...)^2 / a2E(c, scale_factor(z))

#############################################################################################
# times

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
# dark energy
""" 
    w(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM,FlatLCDM,OpenLCDM,ClosedLCDM},z)
    w(z::Real,w0::Real,wa::Real)
Evaluates the redshift dependence of the dark energy density. The scaling factor, I, is defined by :math:`\\rho(z) = \\rho_0 I`.
"""
w(c::Union{FlatLCDM,OpenLCDM,ClosedLCDM},z) = -1.0
w(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM},z) = c.w0 + c.wa * z / (1.0+z)
w(z::Real,w0::Real,wa::Real) = w0 + wa * z / (1+z)

"""
    de_density_scale(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM},z)
    de_density_scale(z::Real,w0::Real,wa::Real)
"""
function de_density_scale(c::Union{FlatWCDM,OpenWCDM,ClosedWCDM},z)
    zp1 = 1.0 + z
    return zp1^(3 * (1 + c.w0 + c.wa)) * exp(-3 * c.wa * z / zp1)
end
de_density_scale(c::Union{FlatLCDM,OpenLCDM,ClosedLCDM},z) = 1.0
de_density_scale(z::Real,w0::Real,wa::Real) = (zp1 = 1+z; zp1^(3 * (1 + w0 + wa)) * exp(-3 * wa * z / zp1))

#############################################################################################
# densities
"""
    ρ_c([u::UnitLike,], c::AbstractCosmology, z)
    ρ_c([u::UnitLike,], z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
The critical density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
"""
ρ_c(c::AbstractCosmology,z) = constants.RHO_C_Z0_CGS * c.h^2 * E(c,z)^2 * u.g/u.cm^3
ρ_c(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = constants.RHO_C_Z0_CGS * h^2 *
    E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2 * u.g/u.cm^3

"""
    ρ_m([u::UnitLike,], c::AbstractCosmology, z)
    ρ_m([u::UnitLike,], z,h,OmegaM)
The matter density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
"""
ρ_m(c::AbstractCosmology,z) = constants.RHO_C_Z0_CGS * c.h^2 * c.Ω_m * (1.0 + z)^3 * u.g/u.cm^3
ρ_m(z,h,Ω_m) = constants.RHO_C_Z0_CGS * h^2 * Ω_m * (1+z)^3 * u.g/u.cm^3

"""
    ρ_b([u::UnitLike,], c::AbstractCosmology, z)
    ρ_b([u::UnitLike,], z,h,OmegaB)
The baryon density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
"""
ρ_b(c::AbstractCosmology,z) = constants.RHO_C_Z0_CGS * c.h^2 * c.Ω_b * (1.0 + z)^3 * u.g/u.cm^3
ρ_b(z,h,Ω_b) = constants.RHO_C_Z0_CGS * h^2 * Ω_b * (1+z)^3 * u.g/u.cm^3

"""
    ρ_dm([u::UnitLike,], c::AbstractCosmology, z)
    ρ_dm([u::UnitLike,], z,h,OmegaDM)
The dark matter density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
"""
ρ_dm(c::AbstractCosmology,z) = constants.RHO_C_Z0_CGS * c.h^2 * c.Ω_dm * (1.0 + z)^3 * u.g/u.cm^3
ρ_dm(z,h,Ω_dm) = constants.RHO_C_Z0_CGS * h^2 * Ω_dm * (1+z)^3 * u.g/u.cm^3

"""
    ρ_Λ([u::UnitLike,], c::AbstractCosmology, z)
    ρ_Λ([u::UnitLike,], z,h,Ω_Λ,w0=-1,wa=0)
The dark energy density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
"""
ρ_Λ(c::AbstractCosmology,z) = constants.RHO_C_Z0_CGS * c.h^2 * c.Ω_Λ / de_density_scale(c,z) * u.g/u.cm^3
ρ_Λ(z,h,Ω_Λ,w0=-1,wa=0) = constants.RHO_C_Z0_CGS * h^2 * Ω_Λ / de_density_scale(z,w0,wa) * u.g/u.cm^3

"""
    ρ_γ([u::UnitLike,], c::AbstractCosmology, z)
    ρ_γ([u::UnitLike,], z,h,Ω_γ)
The photon matter density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
"""
ρ_γ(c::AbstractCosmology,z) = constants.RHO_C_Z0_CGS * c.h^2 * c.Ω_γ * (1.0 + z)^4 * u.g/u.cm^3
ρ_γ(z,h,Ω_γ) = constants.RHO_C_Z0_CGS * h^2 * Ω_γ * (1+z)^4 * u.g/u.cm^3

"""
    ρ_ν([u::UnitLike,], c::AbstractCosmology, z)
    ρ_ν([u::UnitLike,], z,h,Tcmb0,Neff,m_nu,N_nu=nothing)
The neutrino energy density of the universe at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
"""
ρ_ν(c::AbstractCosmology,z) = ρ_γ(c,z) * nu_relative_density(c, z)
ρ_ν(z,h,Tcmb0,Neff,m_nu,N_nu=nothing) = (Ω_γ = 4.481620089297254e-7 * Tcmb0^4 / h^2; ρ_γ(z,h,Ω_γ) * nu_relative_density(m_nu, Neff, u.ustrip(T_nu(Tcmb0, z)), N_nu) )

"""
    ρ_r([u::UnitLike,], c::AbstractCosmology, z)
    ρ_r([u::UnitLike,], z,h,Tcmb0,Neff,m_nu,N_nu=nothing)
The energy density of the universe in relativistic species at redshift z, in g / cm^3. Will convert to compatible unit `u` if provided.
"""
ρ_r(c::AbstractCosmology,z) = ρ_γ(c,z) * (1.0 + nu_relative_density(c::AbstractCosmology, z))
ρ_r(z,h,Tcmb0,Neff,m_nu,N_nu=nothing) = (Ω_γ = 4.481620089297254e-7 * Tcmb0^4 / h^2; ρ_γ(z,h,Ω_γ) * (1.0 + nu_relative_density(m_nu, Neff, u.ustrip(T_nu(Tcmb0, z)), N_nu) ) )


#############################################################################################

""" 
    Ω_m(c::AbstractCosmology,z)
    Ω_m(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Density of matter in units of the critical density. """
Ω_m(c::AbstractCosmology,z) = c.Ω_m * (1 + z)^3 / E(c,z)^2
Ω_m(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = OmegaM * (1+z)^3 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2
""" 
    Ω_dm(c::AbstractCosmology,z)
    Ω_dm(z,h,OmegaM,OmegaB,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Density of dark matter in units of the critical density. """
Ω_dm(c::AbstractCosmology,z) = (c.Ω_m - c.Ω_b) * (1 + z)^3 / E(c,z)^2
Ω_dm(z,h,OmegaM,OmegaB,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = (OmegaM - OmegaB) * (1+z)^3 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2
""" 
    Ω_b(c::AbstractCosmology,z)
    Ω_b(z,h,OmegaM,OmegaB,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Density of baryons in units of the critical density. """
Ω_b(c::AbstractCosmology,z) = c.Ω_b * (1 + z)^3 / E(c,z)^2
Ω_b(z,h,OmegaM,OmegaB,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = OmegaB * (1+z)^3 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2
""" 
    Ω_k(c::AbstractCosmology,z)
    Ω_k(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density of curvature in units of the critical density. """
Ω_k(c::Union{AbstractOpenCosmology,AbstractClosedCosmology},z) = c.Ω_k * (1 + z)^2 / E(c,z)^2
Ω_k(c::AbstractFlatCosmology,z) = 0.0
Ω_k(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = OmegaK * (1+z)^2 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2
""" 
    Ω_γ(c::AbstractCosmology,z)
    Ω_γ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density of photons in units of the critical density. """
Ω_γ(c::AbstractCosmology,z) = c.Ω_γ * (1 + z)^4 / E(c,z)^2
Ω_γ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = (OmegaG = 4.481620089297254e-7 * Tcmb0^4 / h^2; OmegaG * (1+z)^4 / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2 )
""" 
    Ω_ν(c::AbstractCosmology,z)
    Ω_ν(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density in neutrinos in units of the critical density. """
Ω_ν(c::AbstractCosmology,z) = Ω_γ(c,z) * nu_relative_density(c,z)
Ω_ν(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = (nu_temp = u.ustrip(T_nu(Tcmb0, z)); Ω_γ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) * nu_relative_density(m_nu, Neff, nu_temp) )
""" 
    Ω_r(c::AbstractCosmology,z)
    Ω_r(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density in radiation (photons+neutrinos) in units of the critical density. """
Ω_r(c::AbstractCosmology,z) = Ω_γ(c,z) + Ω_ν(c,z)
Ω_r(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = Ω_γ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) + Ω_ν(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)
""" 
    Ω_Λ(c::AbstractCosmology,z)
    Ω_Λ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Energy density in dark energy in units of the critical density. """
Ω_Λ(c::AbstractCosmology,z) = c.Ω_Λ * de_density_scale(c,z) / E(c,z)^2
Ω_Λ(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = OmegaL * de_density_scale(z,w0,wa) / E(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)^2
#############################################################################################
# temperature
"""
    T_cmb([u::Unitlike,], c::AbstractCosmology, z)
    T_cmb([u::Unitlike,], z,Tcmb0)
The temperature of the CMB at redshift z, in Kelvin. Will convert to compatible unit `u` if provided.
"""
T_cmb(c::AbstractCosmology, z) = c.Tcmb0 * (1 + z) * u.K
T_cmb(z,Tcmb0) = (!isa(Tcmb0,u.Quantity) && (Tcmb0 *= u.K); Tcmb0 * (1 + z))

#############################################################################################
# roots
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
# misc
sound_horizon(c::AbstractCosmology) = 44.5 * log(9.83 / c.Ω_m / c.h^2) / sqrt(1.0 + 10.0 * (c.Ω_b * c.h^2)^0.75) * ua.Mpc

#########################################################################################
# Easily select a different unit
for f in (:hubble_dist0, :hubble_dist, :hubble_time0, :hubble_time, :comoving_radial_dist,
          :comoving_transverse_dist, :angular_diameter_dist, :luminosity_dist,
          :comoving_volume, :comoving_volume_element, :age, :lookback_time,
          :T_nu, :T_cmb,:ρ_c, :ρ_m, :ρ_b, :ρ_dm, :ρ_Λ, :ρ_γ, :ρ_ν, :ρ_r, :sound_horizon)
    @eval $f(uu::u.Unitlike, args...; kws...) = u.uconvert(uu, $f(args...; kws...))
end
