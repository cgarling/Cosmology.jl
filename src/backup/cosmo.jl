module cosmo

using QuadGK, Unitful
import Unitful: km, s
using UnitfulAstro: Mpc, Gpc, Gyr

export cosmology,
       age,
       angular_diameter_dist,
       comoving_radial_dist,
       comoving_transverse_dist,
       comoving_volume,
       comoving_volume_element,
       distmod,
       H,
       hubble_dist,
       hubble_time,
       luminosity_dist,
       lookback_time,
       scale_factor


include("types_base.jl")

a2E(c::FlatLCDM, a) = sqrt(c.Ω_r + c.Ω_m * a + c.Ω_Λ * a^4)

function a2E(c::Union{ClosedLCDM,OpenLCDM}, a)
    a2 = a * a
    sqrt(c.Ω_r + c.Ω_m * a + (c.Ω_k + c.Ω_Λ * a2) * a2)
end


function a2E(c::Union{FlatWCDM,ClosedWCDM,OpenWCDM}, a)
    ade = exp((1 - 3 * (c.w0 + c.wa)) * log(a) + 3 * c.wa * (a - 1))
    sqrt(c.Ω_r + (c.Ω_m + c.Ω_k * a) * a + c.Ω_Λ * ade)
end


# hubble rate

scale_factor(z) = 1 / (1 + z)
E(c::AbstractCosmology, z) = (a = scale_factor(z); a2E(c, a) / a^2)
H(c::AbstractCosmology, z) = 100 * c.h * E(c, z) * km / s / Mpc

hubble_dist0(c::AbstractCosmology) = 2997.92458 / c.h * Mpc
hubble_dist(c::AbstractCosmology, z) = hubble_dist0(c) / E(c, z)

hubble_time0(c::AbstractCosmology) = 9.77814 / c.h * Gyr
hubble_time(c::AbstractCosmology, z) = hubble_time0(c) / E(c, z)

# distances

Z(c::AbstractCosmology, z::Real, ::Nothing; kws...) =
    QuadGK.quadgk(a->1 / a2E(c, a), scale_factor(z), 1; kws...)[1]
Z(c::AbstractCosmology, z₁::Real, z₂::Real; kws...) =
    QuadGK.quadgk(a->1 / a2E(c, a), scale_factor(z₂), scale_factor(z₁); kws...)[1]

comoving_radial_dist(c::AbstractCosmology, z₁, z₂ = nothing; kws...) = hubble_dist0(c) * Z(c, z₁, z₂; kws...)

"""
    comoving_radial_dist([u::Unitlike,] c::AbstractCosmology, [z₁,] z₂; kws...)

Comoving radial distance in Mpc at redshift `z₂` as seen by an observer at `z₁`.  Redshift `z₁` defaults to 0 if omitted.  Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
comoving_radial_dist


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

angular_diameter_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...) / (1 + z)
angular_diameter_dist(c::AbstractCosmology, z₁, z₂; kws...) =
    comoving_transverse_dist(c, z₁, z₂; kws...) / (1 + z₂)

"""
    angular_diameter_dist([u::Unitlike,] c::AbstractCosmology, [z₁,] z₂; kws...)

Ratio of the proper transverse size in Mpc of an object at redshift `z₂` to its angular size in radians, as seen by an observer at `z₁`.  Redshift `z₁` defaults to 0 if omitted.  Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
angular_diameter_dist

luminosity_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...) * (1 + z)

"""
    luminosity_dist([u::Unitlike,] c::AbstractCosmology, z; kws...)

Bolometric luminosity distance in Mpc at redshift `z`. Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
luminosity_dist

"""
    distmod(c::AbstractCosmology, z; kws...)

Distance modulus in magnitudes at redshift `z`. kws are integration options passed to quadgk.
"""
distmod(c::AbstractCosmology, z; kws...) =
    5 * log10(luminosity_dist(c, z; kws...) / Mpc) + 25

# volumes

"""
    comoving_volume([u::Unitlike,] c::AbstractCosmology, z; kws...)

Comoving volume in cubic Gpc out to redshift `z`. Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
comoving_volume(c::AbstractFlatCosmology, z; kws...) =
    (4pi / 3) * (comoving_radial_dist(Gpc, c, z; kws...))^3
function comoving_volume(c::AbstractOpenCosmology, z; kws...)
    DH = hubble_dist0(Gpc, c)
    x = comoving_transverse_dist(Gpc, c, z; kws...) / DH
    sqrtok = sqrt(c.Ω_k)
    2pi * (DH)^3 * (x * sqrt(1 + c.Ω_k * x^2) - asinh(sqrtok * x) / sqrtok) / c.Ω_k
end
function comoving_volume(c::AbstractClosedCosmology, z; kws...)
    DH = hubble_dist0(Gpc, c)
    x = comoving_transverse_dist(Gpc, c, z; kws...) / DH
    sqrtok = sqrt(abs(c.Ω_k))
    2pi * (DH)^3 * (x * sqrt(1 + c.Ω_k * x^2) - asin(sqrtok * x) / sqrtok) / c.Ω_k
end

"""
    comoving_volume_element([u::Unitlike,] c::AbstractCosmology, z; kws...)

Comoving volume element in Gpc out to redshift `z`. Will convert to compatible unit `u` if provided. kws are integration options passed to quadgk.
"""
comoving_volume_element(c::AbstractCosmology, z; kws...) =
    hubble_dist0(Gpc, c) * angular_diameter_dist(Gpc, c, z; kws...)^2 / a2E(c, scale_factor(z))

# times

T(c::AbstractCosmology, a0, a1; kws...) = QuadGK.quadgk(x->x / a2E(c, x), a0, a1; kws...)[1]
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

# Easily select a different unit
for f in (:hubble_dist0, :hubble_dist, :hubble_time0, :hubble_time, :comoving_radial_dist,
          :comoving_transverse_dist, :angular_diameter_dist, :luminosity_dist,
          :comoving_volume, :comoving_volume_element, :age, :lookback_time)
    @eval $f(u::Unitful.Unitlike, args...; kws...) = uconvert(u, $f(args...; kws...))
end

end # module
