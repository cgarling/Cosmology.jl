
module peaks

export lagrangianR, lagrangianM, δc, νpeak2, mass_from_νpeak, r_from_νpeak, nonlinear_mass, nonlinear_scale

import UnitfulAstro as ua
import Unitful as u
import Dierckx: Spline1D, derivative
import ..Cosmology: AbstractCosmology, Ω_m, ρ_m, σ2, PkFilter, DefaultFilter
# include("constants.jl")
# import ..transfer: σ2
include("constants.jl")
import .constants: DELTA_COLLAPSE

"""
    lagrangianR(c::AbstractCosmology,M::Union{Real,u.Quantity})
The Lagrangian radius (in Mpc) of a halo of mass `M` (in Msun) at z=0; that is the radius of a sphere that would enclose the halo's mass at the mean density of the universe at z=0. """
lagrangianR(c::AbstractCosmology,M::Real) = cbrt(3 * M / (4π * (ρ_m(c,0.0) |> ua.Msun/ua.Mpc^3 |> u.ustrip) ) )
lagrangianR(c::AbstractCosmology,M::u.Quantity) = lagrangianR(c,M|>ua.Msun|>u.ustrip)

"""
    lagrangianM(c::AbstractCosmology,R::Union{Real,u.Quantity})
The Lagrangian mass of a halo of radius `R`; that is the mass enclosed by a sphere of radius `R` at z=0 with at the mean density of the universe.
"""
lagrangianM(c::AbstractCosmology,R::Real) = 4π * R^3 * (ρ_m(c,0.0) |> ua.Msun/ua.Mpc^3 |> u.ustrip) / 3
lagrangianM(c::AbstractCosmology,R::u.Quantity) = lagrangianR(c,R|>ua.Mpc|>u.ustrip)

"""
    δc(z::Real,c::AbstractCosmology)
    δc(z::Real,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
The linear overdensity threshold for halo collapse. This uses the fitting function from Nakamura & Suto 1997, or Equation 56 in Chisari 2019. 
"""
δc(c::AbstractCosmology,z::Real) = 3 * (12π)^(2/3) / 20 * (1 + 0.012299 * log10(Ω_m(c,z)))
δc(c::AbstractCosmology,z::AbstractArray) = @. 3 * (12π)^(2/3) / 20 * (1 + 0.012299 * log10(Ω_m(c,z)))
δc(z::Real,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = 3 * (12π)^(2/3) / 20 * (1 + 0.012299 * log10(Ω_m(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)))
""" 
    δc(c::AbstractCosmology,M::Union{Real,AbstractArray},z::Real,mx::Real,gx::Real=1.5)
Equations 7--10 in Benson 2013 that give the revised critical density for collapse in warm dark matter cosmologies. `M` is the halo mass, `z` is the redshift of evaluation, `mx` is the mass of the WDM particle in keV, and `gx` is the effective degrees of freedom with 1.5 being the expected value for a fermionic spin-1/2 particle.
"""
function δc(c::AbstractCosmology,M::Union{Real,AbstractArray},z::Real,mx::Real,gx::Real=1.5)
    lcdm_val = δc(c,z)
    zeq = 3600 * (c.Ω_m * c.h^2 / 0.15) - 1
    jeans_mass = 3.06e8 * ((1+zeq)/3000)^1.5 * sqrt(c.Ω_m * c.h^2 / 0.15) * (gx/1.5)^-1 * (1/mx)^4
    
    x = @. log(M/jeans_mass)
    h = @. (1+exp((x+2.4)/0.1))^-1
    return @. lcdm_val * ((h * 0.04 / exp(2.3*x)) + ((1-h) * exp(0.31687/exp(0.809*x))))
end
"""
    νpeak2(R,Pk_func,δc=1.68647; s2=nothing, growth_factor=1.0, k_min=9.99e-3, k_max=1.99e4, j=0, filt::Type{<:PkFilter}=DefaultFilter, kws...)
    νpeak2(R::Real,c::AbstractCosmology,z; s2=nothing, j=0, filt::Type{<:PkFilter}=DefaultFilter, kws...)
    νpeak2(R::Union{Real,AbstractArray},k::AbstractVector,Pk::AbstractVector,deltac=DELTA_COLLAPSE; s2=nothing, growth_factor=1.0, j=0, filt::Type{<:PkFilter}=DefaultFilter)
Peak height squared, ν^2 = δc^2 / σ2(R), where `R` is in Mpc. Options are as in σ2, except for arg δc and kws s2=nothing. δc can be changed (e.g., evolved with z) if desired. If s2 is not nothing, this is taken to be the value of σ2 at radius R so it doesn't need to be recomputed. """
function νpeak2(R::Real,Pk_func::Union{Function,Spline1D},deltac=DELTA_COLLAPSE; s2=nothing, growth_factor=1.0, k_min=9.99e-3, k_max=1.99e4, j=0, filt::Type{<:PkFilter}=DefaultFilter, kws...)
    s2===nothing && (s2 = σ2(R, Pk_func; growth_factor, k_min, k_max, j, filt, kws...) )
    return deltac^2 / s2
end
νpeak2(R::AbstractArray,Pk_func::Union{Function,Spline1D},deltac=DELTA_COLLAPSE; kws...) = [νpeak2(i, Pk_func, deltac; kws...) for i in R]
function νpeak2(R::Union{Real,AbstractArray},c::AbstractCosmology,z=0.0; s2=nothing, j=0, filt=nothing, kws...)
    s2===nothing && (s2 = σ2(R,c,z; j=j, filt=filt, kws...) )
    return δc(c,z).^2 ./ s2
end
function νpeak2(R::Union{Real,AbstractArray},k::AbstractVector,Pk::AbstractVector,deltac=DELTA_COLLAPSE; s2=nothing, growth_factor=1.0, j=0, filt::Type{<:PkFilter}=DefaultFilter)
    s2===nothing && (s2 = σ2(R, k, Pk; growth_factor, j, filt) )
    return deltac^2 ./ s2
end
"""
    mass_from_νpeak(c::AbstractCosmology, ν::Union{Real,AbstractArray}, z::Union{Real,AbstractArray})
Calculates present-day halo mass from peak height at redshift `z`. Peak height is defined as ν(z) = δc(z) / σ(z). This is defined as in νpeak2, so the provided `ν` should be properly scaled to the redshift `z`. Explicitly, 
```julia
    R = lagrangianR(c,M)
    ν = sqrt(νpeak2(R,c,z)) 
    mass_from_νpeak(c,ν,z) ≈ M
```
Currently uses the inverse interpolator saved in the AbstractCosmology object, which uses a `tophat` filter. 
"""
function mass_from_νpeak(c::AbstractCosmology, ν::Union{Real,AbstractArray}, z::Union{Real,AbstractArray})
    sigma = δc(c,z) ./ ν ./ c.growth_function(z)
    R = c.σ2_interp_inv(sigma.^2)
    M = lagrangianM.(c,R)
    return M
end
"""
    r_from_νpeak(c::AbstractCosmology,ν::Union{Real,AbstractArray}, z::Union{Real,AbstractArray})
Calculates the present-day radius (in Mpc/h) corresponding to peak height `ν` at redshift `z`. See also `mass_from_νpeak`. """
function r_from_νpeak(c::AbstractCosmology,ν::Union{Real,AbstractArray}, z::Union{Real,AbstractArray})
    sigma = δc(c,z) ./ ν ./ c.growth_function(z)
    R = c.σ2_interp_inv(sigma.^2)
    return R
end

"""
    nonlinear_mass(c::AbstractCosmology, z::Union{Real,AbstractArray})
The non-linear mass, defined as the mass for which the variance is equal to the collapse threshold,
```julia
    σ(M,z) = δc(z) 
    ν = δc(z) / σ(M,z)
    ν = 1
```
"""
nonlinear_mass(c::AbstractCosmology, z::Union{Real,AbstractArray}) = mass_from_νpeak(c, 1.0, z)
""" 
    nonlinear_scale(c::AbstractCosmology,z::Union{Real,AbstractArray})
Calculates the scale (in Mpc/h) at which the variance is equal to the collapse threshold. See also `nonlinear_mass`."""
nonlinear_scale(c::AbstractCosmology,z::Union{Real,AbstractArray}) = r_from_νpeak(c,1.0,z)


end # module
