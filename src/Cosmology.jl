module Cosmology

import Unitful as u
import UnitfulAstro as ua
import UnitfulEquivalences as ue
import Roots: find_zero, Bisection
import QuadGK: quadgk

export cosmology,
    age,
    angular_diameter_dist,
    comoving_radial_dist,
    comoving_transverse_dist,
    comoving_volume,
    comoving_volume_element,
    distmod,
    H,
    E,
    hubble_dist,
    hubble_time,
    luminosity_dist,
    lookback_time,
    scale_factor, ∇scale_factor,
    ρ_c, ρ_m, ρ_b, ρ_dm, ρ_Λ, ρ_γ, ρ_ν, ρ_r,
    # n_nu,
    T_nu, T_cmb,
    # nu_relative_density,
    # de_density_scale,
    z_at_value,
    Ω_m, Ω_dm, Ω_b, Ω_k, Ω_γ, Ω_ν, Ω_r, Ω_Λ, 
    #m_nu, Neff, h, w,
    sound_horizon, matter_radiation_equality

include("utils.jl")
include("constants.jl")
import ..constants

#### type definitions #########################################################################################
""" `AbstractCosmology` is the base type for all cosmologies. """
abstract type AbstractCosmology end
Base.Broadcast.broadcastable(m::AbstractCosmology) = Ref(m)

""" `AbstractClosedCosmology` is the base type for all closed cosmologies (Ω_k<0). """
abstract type AbstractClosedCosmology <: AbstractCosmology end
""" `AbstractFlatCosmology` is the base type for all flat cosmologies (Ω_k=0). """
abstract type AbstractFlatCosmology <: AbstractCosmology end
""" `AbstractOpenCosmology` is the base type for all open cosmologies (Ω_k<0). """
abstract type AbstractOpenCosmology <: AbstractCosmology end

#############################################################################
"""
    FlatLCDM(h::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν)
Type for flat (Ω_k=0) ΛCDM cosmologies (w0=-1, wa=0). """
struct FlatLCDM{T <: Real, N} <: AbstractFlatCosmology
    h::T
    Ω_Λ::T
    Ω_m::T
    Ω_b::T
    Tcmb0::T # in kelvin
    Neff::T
    m_nu::NTuple{N,T} # If x is a vector of m_nu, ntuple(i->x[i],Val(length(x)))
end                   # Also ntuple(i->x[i],length(x)) if length(x) not known at compile-time. 
function FlatLCDM(h::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν)
    # Promote the scalars
    h, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff = promote(h, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff)
    T = typeof(h)   # Get the type determined by the scalar promotion
    m_ν = m_nu(m_ν) # Convert m_ν to an NTuple
    m_ν = convert(NTuple{length(m_ν),T},m_ν) # Convert the eltype of the NTuple to the correct thing
    return FlatLCDM(h, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, m_ν) 
end
#############################################################################
"""
    ClosedLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν)
Type for closed (Ω_k<0) ΛCDM cosmologies (w0=-1, wa=0). """
struct ClosedLCDM{T <: Real, N} <: AbstractClosedCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_b::T
    Tcmb0::T
    Neff::T
    m_nu::NTuple{N,T}
end
function ClosedLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν)
    # Promote the scalars
    h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff = promote(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff)
    T = typeof(h)   # Get the type determined by the scalar promotion
    m_ν = m_nu(m_ν) # Convert m_ν to an NTuple
    m_ν = convert(NTuple{length(m_ν),T},m_ν) # Convert the eltype of the NTuple to the correct thing
    return ClosedLCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, m_ν) 
end
#############################################################################
"""
    OpenLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν)
Type for open (Ω_k>0) ΛCDM cosmologies (w0=-1, wa=0). """
struct OpenLCDM{T <: Real, N} <: AbstractOpenCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_b::T
    Tcmb0::T
    Neff::T
    m_nu::NTuple{N,T}
end
function OpenLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν)
    # Promote the scalars
    h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff = promote(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff)
    T = typeof(h)   # Get the type determined by the scalar promotion
    m_ν = m_nu(m_ν) # Convert m_ν to an NTuple
    m_ν = convert(NTuple{length(m_ν),T},m_ν) # Convert the eltype of the NTuple to the correct thing
    return OpenLCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, m_ν) 
end
#############################################################################
# the same cosmologies but with w0-wa-state dark energy
for c in ("Flat", "Open", "Closed")
    name = Symbol("$(c)WCDM")
    @eval begin
        """
                FlatWCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν, w0::Real, wa::Real)
            $($c) cosmology with the w0wa dark energy equation of state,
            ```math
            w(z) = w_0 + w_a \\, \\frac{z}{1+z} = w_0 + w_a \\, (1-a)
            ```
            (Equations 6 and 7 in [Linder 2003](https://ui.adsabs.harvard.edu/abs/2003PhRvL..90i1301L/abstract)). """
        struct $(name){TTT <: Real, N} <: $(Symbol("Abstract$(c)Cosmology"))
            h::TTT
            Ω_k::TTT
            Ω_Λ::TTT
            Ω_m::TTT
            Ω_b::TTT
            Tcmb0::TTT
            Neff::TTT
            m_nu::NTuple{N,TTT}
            w0::TTT
            wa::TTT
        end
    end
end

function FlatWCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν, w0::Real, wa::Real)
    @assert iszero(Ω_k)
    # Promote the scalars
    h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, w0, wa = promote(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, w0, wa)
    T = typeof(h)   # Get the type determined by the scalar promotion
    m_ν = m_nu(m_ν) # Convert m_ν to an NTuple
    m_ν = convert(NTuple{length(m_ν),T},m_ν) # Convert the eltype of the NTuple to the correct thing
    return FlatWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, m_ν, w0, wa)     
end

function OpenWCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν, w0::Real, wa::Real)
    @assert Ω_k > 0
    # Promote the scalars
    h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, w0, wa = promote(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, w0, wa)
    T = typeof(h)   # Get the type determined by the scalar promotion
    m_ν = m_nu(m_ν) # Convert m_ν to an NTuple
    m_ν = convert(NTuple{length(m_ν),T},m_ν) # Convert the eltype of the NTuple to the correct thing
    return OpenWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, m_ν, w0, wa)     
end

function ClosedWCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν, w0::Real, wa::Real)
    @assert Ω_k < 0
    # Promote the scalars
    h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, w0, wa = promote(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, w0, wa)
    T = typeof(h)   # Get the type determined by the scalar promotion
    m_ν = m_nu(m_ν) # Convert m_ν to an NTuple
    m_ν = convert(NTuple{length(m_ν),T},m_ν) # Convert the eltype of the NTuple to the correct thing
    return OpenWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, m_ν, w0, wa)     
end

"""
    WCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν, w0::Real, wa::Real)
Constructs a cosmology with the w0wa dark energy equation of state,
```math
w(z) = w_0 + w_a \\, \\frac{z}{1+z} = w_0 + w_a \\, (1-a)
```
and returns one of [`FlatWCDM`](@ref Cosmology.FlatWCDM), [`OpenWCDM`](@ref Cosmology.OpenWCDM), [`ClosedWCDM`](@ref Cosmology.ClosedWCDM) depending the value of `Ω_k`. This function is type-unstable by design. 
"""
function WCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν, w0::Real, wa::Real)
    # Promote the scalars
    h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, w0, wa = promote(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, w0, wa)
    T = typeof(h)   # Get the type determined by the scalar promotion
    m_ν = m_nu(m_ν) # Convert m_ν to an NTuple
    m_ν = convert(NTuple{length(m_ν),T},m_ν) # Convert the eltype of the NTuple to the correct thing
    return WCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, m_ν, w0, wa) 
end
function WCDM(h::T, Ω_k::T, Ω_Λ::T, Ω_m::T, Ω_b::T, Tcmb0::T, Neff::T, m_nu::NTuple{N,T},
              w0::T, wa::T) where {T<:Real,N}
    if Ω_k < 0
        ClosedWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, m_nu, w0, wa)
    elseif Ω_k > 0
        OpenWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, m_nu, w0, wa)
    else
        FlatWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_b, Tcmb0, Neff, m_nu, w0, wa)
    end
end

#############################################################################
include("standard_functions.jl")
#############################################################################

"""
    m_nu(x::AbstractArray)
    m_nu(x::AbstractArray{<:Unitful.Energy})
    m_nu(x::AbstractArray{<:Unitful.Mass})
    m_nu(x::Tuple)
    m_nu(x::NTuple)
    m_nu(x::NTuple{N,T}) where {N,T<:u.Energy}
    m_nu(x::NTuple{N,T}) where {N,T<:u.Mass}
    m_nu(x::Number)
    m_nu(x::Nothing)
Converts various forms of `m_nu` to `NTuple`s so they can be stored in the cosmology structs. Currently mixed units in a `Tuple` is not supported. If providing units in the argument, they should all be the same. 

# Examples
```jldoctest; setup = :(import Cosmology: m_nu; import Unitful)
julia> x=[0.0,0.0,0.06];

julia> m_nu(x)
(0.0, 0.0, 0.06)

julia> x=[0.0,0.0,0.06]*Unitful.eV;

julia> m_nu(x)
(0.0, 0.0, 0.06)

julia> x=(0.0,0.0,0.06);

julia> m_nu(x) === x
true

julia> x=(0.0,0.0,0.06) .* Unitful.eV;

julia> m_nu(x)
(0.0, 0.0, 0.06)

julia> m_nu(0.06)
(0.06,)

julia> m_nu(nothing)
(0,)

```
"""
m_nu(x::AbstractArray) = ntuple(i->x[i],length(x))
m_nu(x::AbstractArray{<:u.Energy}) = ntuple(i->u.ustrip(u.eV,x[i]),length(x))
m_nu(x::AbstractArray{<:u.Mass}) = ntuple(i->u.ustrip(u.eV,x[i],ue.MassEnergy()),length(x))
m_nu(x::Tuple) = m_nu(promote(x...)) # this doesn't work well for mixed units in m_nu. 
m_nu(x::NTuple) = x
m_nu(x::NTuple{N,T}) where {N,T<:u.Energy} = ntuple(i->u.ustrip(u.eV,x[i]),Val(N))
m_nu(x::NTuple{N,T}) where {N,T<:u.Mass} = ntuple(i->u.ustrip(u.eV,x[i],ue.MassEnergy()),Val(N))
m_nu(x::Number) = (x,)
m_nu(x::Nothing) = (0,)

#############################################################################
# Main cosmology function
#############################################################################

"""
    cosmology(;h = 0.6766, OmegaK = 0, OmegaM = 0.30966, OmegaB = 0.04897, OmegaG = nothing, Tcmb0 = 2.7255, w0 = -1, wa = 0, N_eff = 3.046, m_ν=(0.0,0.0,0.06))

Constructs the proper `AbstractCosmology` type depending on the passed parameters. To specify the photon density of the cosmology, either `OmegaG` or `Tcmb0` can be provided; if `OmegaG` is provided, it takes precedence over Tcmb0. To turn neutrinos off, you can set `N_eff=nothing`, `m_ν=nothing`; then `OmegaG==OmegaR`, all radiation is in the form of photons.  

# Parameters
 - `h` - Dimensionless Hubble constant
 - `OmegaK` - Curvature density (Ω_k)
 - `OmegaM` - Matter density (Ω_m)
 - `OmegaB` - Baryon density (Ω_b)
 - `OmegaG` - Photon density (Ω_γ)  
 - `Tcmb0` - CMB temperature in Kelvin; used to compute Ω_γ if not provided
 - `w0` - CPL dark energy equation of state; `w = w0 + wa(1-a)`
 - `wa` - CPL dark energy equation of state; `w = w0 + wa(1-a)`
 - `N_eff` - Effective number of massless neutrino species; used to compute Ω_ν
 - `m_ν` - Neutrino masses; should have length equal to `Int(floor(N_eff))` but this is not explicitly checked. If `N_eff==0` then `m_ν` isn't used for anything.

# Examples

!!! info
    This function is not type-stable; it will return different concrete subtypes of `AbstractCosmology` depending on the parameters you pass. This incurs a performance penalty. If you want to, you can access constructors of the concrete types directly; this can make creation of new instances much faster. For example, `c=cosmology()` takes ~300 ns, while `FlatLCDM(c.h,c.Ω_Λ,c.Ω_m,c.Ω_b,c.Tcmb0,c.Neff,c.m_nu)` takes ~1 ns. If you REALLY need this speed, you should use the basic constructors.
!!! note
    Inclusion of massive neutrinos is expensive. For example, for the default massive neutrino parameters `c=cosmology()`, the evaluation of `E(c, 0.8)` takes 114.613 ns, while `E( cosmology(m_ν=(0.0,),N_eff=3.046), 0.8)` takes 6.986 ns and `E( cosmology(m_ν=(0.0,),N_eff=0), 0.8)` takes 6.095 ns. This makes a significant difference in methods that involve integrations (e.g., [`comoving_radial_dist`](@ref)). If speed is a concern, consider if you can neglect neutrinos for your calculation. """
function cosmology(;h::Number = 0.6766,         # Scalar; Hubble constant at z = 0 / 100 [km/sec/Mpc]
                   OmegaK::Number = 0,          # Energy density of curvature in units of the critical density at z=0. 
                   OmegaM::Number = 0.30966,    # Energy density of matter in units of the critical density at z=0. 
                   OmegaB::Number = 0.04897,    # Energy density of matter in units of the critical density at z=0. 
                   OmegaG::Union{Number,Nothing} = nothing,    # Energy density of photons in units of the critical density at z=0. 
                   Tcmb0::Union{Number,Nothing} = 2.7255,      # Temperature of CMB at z=0 in Kelvin; can have units
                   w0 = -1,             
                   wa = 0,
                   N_eff::Union{Real,Nothing} = 3.046,
                   m_ν=(0.0,0.0,0.06))

    # Functions for dispatch and type checking
    # Resolve Ω_γ from OmegaG and Tcmb0
    Ω_γ0(h,OmegaG::Nothing,Tcmb0::Nothing) = zero(h)
    Ω_γ0(h,OmegaG::Number,Tcmb0::Nothing) = OmegaG
    Ω_γ0(h,OmegaG::Nothing,Tcmb0::Number) = 4.481620089297254e-7 * Tcmb0^4 / h^2
    function Ω_γ0(h,OmegaG::Number,Tcmb0::Number)
        if !isapprox(OmegaG,4.481620089297254e-7 * Tcmb0^4 / h^2)
            @warn "Both Tcmb0 and OmegaG are provided but they do not agree; OmegaG will take precedence."
        end
        return OmegaG
    end
    # Resolve Tcmb0 from OmegaG and Tcmb0; OmegaG takes precedence and mutates Tcmb0
    T_cmb0(h,Tcmb0::Nothing,OmegaG::Number) = sqrt(sqrt(OmegaG * h^2 * 2.2313359456508644e6))
    T_cmb0(h,Tcmb0::Number,OmegaG::Number) = T_cmb0(h, nothing, Ω_γ0(h,OmegaG,Tcmb0))
    T_cmb0(h,Tcmb0::u.Temperature,OmegaG::Number) = T_cmb0(h,u.ustrip(u.K,Tcmb0),OmegaG)
    T_cmb0(h,Tcmb0::Number,OmegaG::Nothing) = Tcmb0
    T_cmb0(h,Tcmb0::u.Temperature,OmegaG::Nothing) = u.ustrip(u.K,Tcmb0)

    ######################################################################################
    # Domain checks; not completely done. 
    OmegaM < 0 && throw(DomainError(OmegaM,"OmegaM must be greater than or equal to 0."))
    OmegaB < 0 && throw(DomainError(OmegaB,"OmegaB must be greater than or equal to 0."))
    ((Tcmb0 !== nothing) && (Tcmb0 < zero(Tcmb0))) && throw(DomainError(Tcmb0,"Tcmb0 must be greater than or equal to 0."))
    OmegaB < 0 && throw(DomainError(OmegaB,"OmegaB must be greater than or equal to 0."))
    OmegaB > OmegaM && throw(DomainError(OmegaB,"OmegaB cannot be greater than OmegaM."))
    N_eff < 0 && throw(DomainError(N_eff,"Effective number of neutrino species cannot be negative."))
  
    # Do some type conversions; see above definitions
    m_ν = m_nu(m_ν)
    (!iszero(N_eff) && (length(m_ν) != n_nu(N_eff))) && @warn "Length of `m_ν` should typically be equal to `Int(floor(N_eff))` but is not. Did you make a mistake?" 
    Tcmb0 = T_cmb0(h,Tcmb0,OmegaG)
    OmegaG = Ω_γ0(h,OmegaG,Tcmb0)

    iszero(OmegaG) ? (OmegaR=zero(OmegaM)) : (OmegaR=OmegaG * (1 + nu_relative_density(m_ν, N_eff, u.ustrip(u.K,T_nu(Tcmb0)), n_nu(N_eff))))
    OmegaL = 1 - OmegaK - OmegaM - OmegaR # calculate the dark energy density
    ## Type conversions ###############################################################################################
    # I think I am moving these the to the other type constructors? Not 100% sure if I need to do a type conversion here or not
    ### Initializing structs; this is not type-stable but that's okay
    if !(w0 == -1 && wa == 0)
        return WCDM(h, OmegaK, OmegaL, OmegaM, OmegaB, Tcmb0, N_eff, m_ν, w0, wa)
    end
    if OmegaK < 0
        return ClosedLCDM(h, OmegaK, OmegaL, OmegaM, OmegaB, Tcmb0, N_eff, m_ν)
    elseif OmegaK > 0
        return OpenLCDM(h, OmegaK, OmegaL, OmegaM, OmegaB, Tcmb0, N_eff, m_ν)
    else
        return FlatLCDM(h, OmegaL, OmegaM, OmegaB, Tcmb0, N_eff, m_ν)
    end
end

#############################################################################
include("default_cosmologies.jl")
#############################################################################

end # module
