module halo

import Unitful as u
import UnitfulAstro as ua
import QuadGK: quadgk
import Roots: find_zero, Bisection
import Dierckx: Spline1D
import FiniteDifferences: central_fdm # FiniteDiff is faster, but less flexible.
c_fdm = central_fdm(3,1)
# import LoopVectorization: @turbo
using ..cosmo: AbstractCosmology, Ω_m, ρ_c, ρ_m, hubble_time
include("utils.jl")
include("constants.jl")
import .constants: G, KPC_KM

export NFW, ρ, ∇ρ, dlogρ_dlogr, ρmean, ∇ρmean, enclosed_mass, Σ, projected_mass, Φ, ∇Φ, Δvir, M_to_R, R_to_M, dynamical_time, ∇∇Φ, ∇²Φ, circular_velocity, escape_velocity, V_max, RΔ, MΔ, evolve_so, convert_mdef, convert_z


abstract type AbstractHaloProfile end
Base.Broadcast.broadcastable(m::AbstractHaloProfile) = Ref(m) # maybe look into the broadcasting rules
# https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting

# struct NFW{T <: Union{Real,u.Quantity,Array{<:Real},Array{<:u.Quantity}},
#            S <: Union{Real,u.Quantity,Array{<:Real},Array{<:u.Quantity}}} <: AbstractHaloProfile
#     rs::T
#     rhos::S
# end
struct NFW{T <: Union{Real,AbstractArray}, S <: Union{Real,AbstractArray}, Q <: Union{Real,AbstractArray}, V <: Union{Real,AbstractArray}} <: AbstractHaloProfile
    rs::T
    rhos::S
    conc::V
    mdef::String
    z::Q
end

"""
    NFW(; c::Union{AbstractCosmology,Nothing}=nothing, rs::Union{Nothing,Real,u.Quantity,Array{<:Real},Array{<:u.Quantity}}=nothing, rhos::Union{Nothing,Real,u.Quantity,Array{<:Real},Array{<:u.Quantity}}=nothing, conc::Union{Nothing,Real,Array{<:Real}}=nothing,M::Union{Nothing,Real,Array{<:Real},u.Quantity,Array{<:u.Quantity}}=nothing,mdef::Union{AbstractString,Nothing}=nothing,z::Union{Real,Array{<:Real}})
Constructor for the NFW struct, which is used as input for all the methods below. Can either specify rs and rhos directly (assumed kpc and Msun, unless Unitful units as provided; must be same length), or M (halo mass), conc (halo concentration), c (AbstractCosmology object), and z (redshift). If using the mass interface, there are options for lengths of M, conc, and z. Each can be scalar; M and conc can be of same length and z can be scalar; M and conc can be scalar and z can have length; or M, conc, and z can all have same length. Available functions for the NFW struct are
```julia
RΔ(c::NFW,cosmo::AbstractCosmology,z,mdef)
MΔ(c::NFW,cosmo::AbstractCosmology,z,mdef)
evolve_so(c::NFW, cosmo::AbstractCosmology, z, mdef)
convert_mdef(c::NFW, cosmo::AbstractCosmology, mdef)
convert_z(c::NFW, cosmo::AbstractCosmology, z)
ρ(r::Union{Real,u.Quantity},c::NFW) 
∇ρ(r::Union{Real,u.Quantity},c::NFW) 
dlogρ_dlogr(r::Union{Real,u.Quantity},c::NFW)
ρmean(r::Union{Real,u.Quantity},c::NFW)
∇ρmean(r::Union{Real,u.Quantity},c::NFW)
enclosed_mass(r::Union{Real,u.Quantity},c::NFW)
Σ(r::Union{Real,u.Quantity},c::NFW)
projected_mass(r::Union{Real,u.Quantity},c::NFW)
Φ(r::Union{Real,u.Quantity},c::NFW)
∇Φ(r::Union{Real,u.Quantity},c::NFW)
∇∇Φ(r::Union{Real,u.Quantity},c::NFW)
∇²Φ(r::Union{Real,u.Quantity},c::AbstractHaloProfile)
circular_velocity(r::Union{Real,u.Quantity},c::AbstractHaloProfile)
escape_velocity(r::Union{Real,u.Quantity},c::AbstractHaloProfile)
V_max(c::NFW)
```
"""
function NFW(;c::Union{AbstractCosmology,Nothing}=nothing, rs::Union{Nothing,Real,u.Quantity,Array{<:Real},Array{<:u.Quantity}}=nothing, rhos::Union{Nothing,Real,u.Quantity,Array{<:Real},Array{<:u.Quantity}}=nothing, conc::Union{Nothing,Real,Array{<:Real}}=nothing,M::Union{Nothing,Real,Array{<:Real},u.Quantity,Array{<:u.Quantity}}=nothing,mdef::Union{AbstractString,Nothing}=nothing,z::Union{Real,Array{<:Real}}=0.0)
    # convert M, concentration, and mdef to rs and rhos.
    if M!==nothing && conc!==nothing && mdef !== nothing && c!==nothing
        if eltype(M) <: u.Quantity # deal with units
            M = M .|> ua.Msun .|> u.ustrip
        end
        if isscalar(M) && isscalar(conc) && isscalar(z)
            rs,rhos=NFW_fundamental_parameters(c, M, conc, z, mdef)
        elseif isscalar(z) && !isscalar(M) && !isscalar(conc)
            len=length(M)
            len != length(conc) && throw(DomainError(conc,"M and conc must have same length."))
            rho_thresh = density_threshold(c,z,mdef)
            p_type=Base.promote_eltype(M,conc,z,rho_thresh,[1.0])
            rs=zeros(p_type,len) 
            rhos=zeros(p_type,len)
            @inbounds for i in 1:len
                rs[i], rhos[i] = NFW_fundamental_parameters(M[i], conc[i], rho_thresh)
            end
        elseif !isscalar(z) && isscalar(M) && isscalar(conc)
            p_type=Base.promote_eltype(M,conc,z,[1.0])
            len=length(z)
            rs=zeros(p_type,len) 
            rhos=zeros(p_type,len) 
            @inbounds for i in 1:len
                rs[i], rhos[i] = NFW_fundamental_parameters(c, M, conc, z[i], mdef)
            end
        elseif !isscalar(z) && !isscalar(M) && !isscalar(conc)
            p_type=Base.promote_eltype(M,conc,z,[1.0])
            len=length(M)
            len != length(conc) || len != length(z) && throw(DomainError(conc,"if z is an Array, M, z, and conc must have same length."))
            rs=zeros(p_type,len) 
            rhos=zeros(p_type,len) 
            @inbounds for i in 1:len
                rs[i], rhos[i] = NFW_fundamental_parameters(c, M[i], conc[i], z[i], mdef)
            end  
        end
    
    elseif rs!==nothing && rhos!==nothing && mdef !== nothing && c!==nothing && conc!==nothing # if rs and rhos are provided
        length(rs)!=length(rhos) && throw(DomainError(rs,"rs and rhos must have same length"))
        # convert units to make following functions simpler; the NFW struct saves rs in kpc and rhos in Msun/kpc^3.
        eltype(rs) <: u.Quantity && (rs = (@. rs |> ua.kpc |> u.ustrip))
        eltype(rhos) <: u.Quantity && (rhos = (@. rhos |> (ua.Msun/ua.kpc^3) |> u.ustrip))
        tt = Base.promote_eltype(rs,rhos)
        rs = tt.(rs)
        rhos = tt.(rhos)
    else
        error("Incorrect input values.")
    end
    # rs,rhos=promote_eltype(rs,rhos)
    return NFW(rs,rhos,conc,mdef,z)
end

"""
    NFW_fundamental_parameters(c::AbstractCosmology, M::Union{Real,AbstractVector}, conc::Union{Real,AbstractVector}, z::Real, mdef::String)
    NFW_fundamental_parameters(M::Real, conc::Real, rho::Real)
    NFW_fundamental_parameters(M::AbstractVector, conc::AbstractVector, z::Real, mdef::String, h, OmegaM, OmegaK, OmegaL, Tcmb0, m_nu, Neff, w0=-1, wa=0)
The fundamental NFW parameters, rhos and rs, calculated from (mass, concentration, z, and mdef) or (mass, concentration, and the density threshold (rho)), which is used for multiple masses but uniform z, so density threshold doesn't get recomputed. Mass should be in Msun, or have a Unitful unit, conc is the concentration, defined as R(mdef)/R_s, z is redshift, and mdef is the mass definition (see parse_mass_definition). """
function NFW_fundamental_parameters(c::AbstractCosmology, M::Real, conc::Real, z::Real, mdef::String)
    M isa u.Quantity && (M = M |> ua.Msun |> u.ustrip)
    rs = M_to_R(c, M, z, mdef) / conc
    rhos = M / rs^3 / 4.0 / π / NFWmu(conc) 
    return rs, rhos
end
function NFW_fundamental_parameters(c::AbstractCosmology, M::AbstractVector, conc::AbstractVector, z::Real, mdef::String)
    @assert length(M) == length(conc)
    eltype(M) <: u.Quantity && (M .= M .|> ua.Msun .|> u.ustrip)
    rs,rhos=similar(M),similar(M)
    for i in eachindex(M,conc)
        rs[i] = M_to_R(c, M[i], z, mdef) / conc[i]
        rhos[i] = M[i] / rs[i]^3 / 4.0 / π / NFWmu(conc[i])
    end
    return rs, rhos
end
NFW_fundamental_parameters(M::Real, conc::Real, rho::Real) = (M isa u.Quantity && (M = M |> ua.Msun |> u.ustrip); rs=M_to_R(M,rho)/conc; rhos=M / rs^3 / 4.0 / π / NFWmu(conc); return rs, rhos)
function NFW_fundamental_parameters(M::Real, conc::Real, z::Real, mdef::String, h, OmegaM, OmegaK, OmegaL, Tcmb0, m_nu, Neff, w0=-1, wa=0)
    M isa u.Quantity && (M = M |> ua.Msun |> u.ustrip)
    rs = M_to_R(M,z,mdef,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) / conc
    rhos = M / rs^3 / 4.0 / π / NFWmu(conc) 
    return rs, rhos
end
function NFW_fundamental_parameters(M::AbstractVector, conc::AbstractVector, z::Real, mdef::String, h, OmegaM, OmegaK, OmegaL, Tcmb0, m_nu, Neff, w0=-1, wa=0)
    @assert length(M) == length(conc)
    M isa u.Quantity && (M = M |> ua.Msun |> u.ustrip)
    rs,rhos = similar(M), similar(M)
    for i in eachindex(M,conc)
        rs[i] = M_to_R(M[i],z,mdef,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) / conc[i]
        rhos[i] = M[i] / rs[i]^3 / 4.0 / π / NFWmu(conc[i])
    end
    return rs, rhos
end

"""
    Δvir(c::AbstractCosmology,z)
    Δvir(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
The virial overdensity in units of the critical density. This function uses the fitting formula of Bryan & Norman 1998 to determine the virial overdensity. While the universe is dominated by matter, this overdensity is about 178. Once dark energy starts to matter, it decreases. `z` is the redshift. """
Δvir(c::AbstractCosmology,z) = (x = Ω_m(c,z) - 1.0; 18 * π^2 + 82.0 * x - 39.0 * x^2)
Δvir(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0) = (x = Ω_m(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) - 1.0; 18 * π^2 + 82.0 * x - 39.0 * x^2)

"""
    parse_mass_definition(mdef)
Parses the type and overdensity of a given spherical overdensity mass definition. mdef is a string, either 'c', 'm', or "vir" for critical density, mean density, or virial density. Returns a String, either "crit", "mean", or "vir", and the overdensity factor, which is an Int64. """
@views function parse_mass_definition(mdef::String)
    if mdef[end] === 'c'
        return "crit", parse(Int64,mdef[1:end-1])
    elseif mdef[end] === 'm'
        return "mean", parse(Int64,mdef[1:end-1])
    elseif mdef === "vir"
        return "vir", 0
    else
        Base.throw(Base.DomainError(mdef,"Invalid mass definition"))
    end
end

"""
    density_threshold(c::AbstractCosmology,z,mdef)
    density_threshold(z,mdef::String,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
The threshold density for a given spherical overdensity mass definition. `z` is the redshift, `mdef` is the mass definition (see parse_mass_definition). Returns the threshold density in Msun/kpc^3. """
function density_threshold(c::AbstractCosmology,z,mdef)
    mdef_type, mdef_delta = parse_mass_definition(mdef)
    if mdef_type === "crit"
	return mdef_delta * ρ_c(ua.Msun/ua.kpc^3,c,z) |> u.ustrip
    elseif mdef_type === "mean"
	return mdef_delta * ρ_m(ua.Msun/ua.kpc^3,c,z) |> u.ustrip
    elseif mdef_type === "vir"
	return Δvir(c,z) * ρ_c(ua.Msun/ua.kpc^3,c,z) |> u.ustrip
    else
        Base.throw(Base.DomainError(mdef,"Invalid mass definition"))
    end
end
function density_threshold(z,mdef::String,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
    mdef_type, mdef_delta = parse_mass_definition(mdef)
    if mdef_type === "crit"
	return mdef_delta * ρ_c(ua.Msun/ua.kpc^3,z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,
                                Neff,w0,wa) |> u.ustrip
    elseif mdef_type === "mean"
	return mdef_delta * ρ_m(ua.Msun/ua.kpc^3,z,h,Ω_m) |> u.ustrip
    elseif mdef_type === "vir"
	return Δvir(z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) *
            ρ_c(ua.Msun/ua.kpc^3,z,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa) |> u.ustrip
    else
        Base.throw(Base.DomainError(mdef,"Invalid mass definition"))
    end
end
"""
    M_to_R(c::AbstractCosmology, M, z, mdef)
    M_to_R(M::Real, rho::Real)
    M_to_R(M,z,mdef::String,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Convert spherical overdensity mass to radius. `M` is the halo mass in Msun (can be Unitful Quantity), `z` is the redshift, `mdef` is a String specifying the mass definition (see parse_mass_definition). Alternatively, you can pass only the mass and the density threshold value. Returns the halo radius for the mass definition in kpc.
"""
function M_to_R(c::AbstractCosmology, M, z, mdef)
    M isa u.Quantity && (M = M |> ua.Msun |> u.ustrip)
    rho = density_threshold(c, z, mdef)
    R = cbrt(M * 3.0 / 4.0 / π / rho)
    return R
end
M_to_R(M::Real, rho::Real) = (M isa u.Quantity && (M = M |> ua.Msun |> u.ustrip); cbrt(M * 3.0 / 4.0 / π / rho))
function M_to_R(M,z,mdef::String,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
    M isa u.Quantity && (M = M |> ua.Msun |> u.ustrip)
    rho = density_threshold(z,mdef,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)
    R = cbrt(M * 3.0 / 4.0 / π / rho)
    return R
end
"""
    R_to_M(c::AbstractCosmology, R, z, mdef)
    R_to_M(R, rho)
    R_to_M(R,z,mdef::String,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
Convert spherical overdensity radius to mass. R is the halo radius in kpc (can be Unitful Quantity), z is the redshift (Number), mdef is a string specifying the mass definition (see parse_mass_definition). Alternatively, you can pass only the radius and the density threshold value. Returns the halo mass in Msun. """
function R_to_M(c::AbstractCosmology, R, z, mdef)
    R isa u.Quantity && (R = R |> ua.kpc |> u.ustrip)
    rho = density_threshold(c, z, mdef)
    M = 4.0 / 3.0 * π * rho * R^3
    return M
end
R_to_M(R, rho) = (R isa u.Quantity && (R = R |> ua.kpc |> u.ustrip); 4.0 / 3.0 * π * rho * R^3)
function R_to_M(R,z,mdef::String,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0=-1,wa=0)
    R isa u.Quantity && (R = R |> ua.kpc |> u.ustrip)
    rho = density_threshold(z,mdef,h,OmegaM,OmegaK,OmegaL,Tcmb0,m_nu,Neff,w0,wa)
    M = 4.0 / 3.0 * π * rho * R^3
    return M
end
NFWmu(r,rs) = (x=r/rs; log(1+x) - x / (1+x))
NFWmu(x) = log(1+x) - x / (1+x)


#########################################################################################################
##### functions on AbstractHaloProfile 's


#### radius and mass conversions
const NFW_RΔ_xarr = logspace(1e4,1e-4,length=1000) 
const NFW_RΔ_interp = Spline1D((@. NFWmu(NFW_RΔ_xarr) * 3 / NFW_RΔ_xarr^3),
                            NFW_RΔ_xarr; k=3, bc="error", s=0.0)

"""
    RΔ(c::NFW,cosmo::AbstractCosmology,z,mdef)
Find the r where the enclosed density has a particular value. Takes an NFW object, an AbstractCosmology object, a redshift, and a mass definition (see parse_mass_definition). See also NFWRΔ.
"""
function RΔ(c::NFW,cosmo::AbstractCosmology,z::Union{Real,AbstractArray},mdef::String)
    rho_thresh = density_threshold.(cosmo,z,mdef) #density_threshold(c::AbstractCosmology,z,mdef)
    # x_arr = logspace(1e4,1e-4,length=1000)
    # y_arr = @. NFWmu(x_arr) * 3 / x_arr^3
    # x_interp = Spline1D(y_arr,x_arr; k=3, bc="error", s=0.0)
    # y = @. rho_thresh / c.rhos
    # return x_interp(y) .* c.rs
    y = rho_thresh ./ c.rhos
    return NFW_RΔ_interp(y) .* c.rs
end

"""
    NFWRΔ(rs::Union{Real,AbstractArray},rhos::Union{Real,AbstractArray},cosmo::AbstractCosmology,z::Union{Real,AbstractArray},mdef::String)
Convenience function for concentration to provide NFW structural parameters and compute the r where the enclosed density has a particular value. See also RΔ.
"""
function NFWRΔ(rs::Union{Real,AbstractArray},rhos::Union{Real,AbstractArray},cosmo::AbstractCosmology,z::Union{Real,AbstractArray},mdef::String)
    rho_thresh = density_threshold.(cosmo,z,mdef) #density_threshold(c::AbstractCosmology,z,mdef)
    # x_arr = logspace(1e4,1e-4,length=1000)
    # y_arr = @. NFWmu(x_arr) * 3 / x_arr^3
    # x_interp = Spline1D(y_arr,x_arr; k=3, bc="error", s=0.0)
    # y = @. rho_thresh / rhos
    # return x_interp(y) .* rs
    y = rho_thresh ./ rhos
    return NFW_RΔ_interp(y) .* rs
end

"""
    MΔ(c::AbstractHaloProfile,cosmo::AbstractCosmology,z,mdef)
Find the spherical overdensity mass for a given mass definition. Takes an NFW object, an AbstractCosmology object, a redshift, and a mass definition (see parse_mass_definition). Calls RΔ, which must have a method for the given AbstractHaloProfile.
"""
MΔ(c::AbstractHaloProfile,cosmo::AbstractCosmology,z,mdef) = R_to_M.(cosmo,RΔ(c,cosmo,z,mdef),z,mdef)
NFWMΔ(rs::Union{Real,AbstractArray},rhos::Union{Real,AbstractArray},cosmo::AbstractCosmology,z,mdef) = R_to_M.(cosmo,NFWRΔ(rs,rhos,cosmo,z,mdef),z,mdef)

"""
    evolveSO(c::AbstractHaloProfile, cosmo::AbstractCosmology, z, mdef)
Evolve the spherical overdensity radius for a fixed density profile. This function computes the evolution of spherical overdensity mass and radius due to a changing  reference density, redshift, or both. The base AbstractHaloProfile should contain the mass, concentration, redshift, and mass definition used to compute the structural parameters contained in the array. New structural parameters will be computed at the new redshift and mass definition provided.
Returns a tuple of M_new, R_new, and c_new. See convert_mdef or convert_z to return updated structs.  """
function evolve_so(c::AbstractHaloProfile, cosmo::AbstractCosmology, z, mdef)
    rho_thresh = density_threshold(cosmo, z, mdef)
    R_new = RΔ(c,cosmo,z,mdef)
    c_new = R_new ./ c.rs
    M_new = R_to_M(R_new, rho_thresh)
    return M_new, R_new, c_new
end

"""
    convert_mdef(c::NFW, cosmo::AbstractCosmology, mdef)
Convert an AbstractHaloProfile, given a cosmology object, and convert it to the same profile with a different mass definition. """
function convert_mdef(c::NFW, cosmo::AbstractCosmology, mdef)
    rho_thresh = density_threshold(cosmo, c.z, mdef)
    R_new = RΔ(c,cosmo,c.z,mdef)
    c_new = R_new ./ c.rs
    M_new = R_to_M.(R_new, Ref(rho_thresh))
    return NFW(c=cosmo,M=M_new, conc=c_new, z=c.z, mdef=mdef)
end

""" Convenience function used in concentration. """
function convert_mdef_NFW(rs::Union{Real,AbstractArray},rhos::Union{Real,AbstractArray},cosmo::AbstractCosmology,z::Real, mdef::String)
    rho_thresh = density_threshold.(cosmo, z, mdef)
    R_new = NFWRΔ.(rs,rhos,cosmo,z,mdef)
    c_new = R_new ./ rs
    M_new = R_to_M.(R_new, rho_thresh)
    return M_new, R_new, c_new
end

"""
    convert_z(c::NFW, cosmo::AbstractCosmology, z)
Convert an AbstractHaloProfile, given a cosmology object, and convert it to the same profile at a different redshift. """
function convert_z(c::NFW, cosmo::AbstractCosmology, z)
    rho_thresh = density_threshold(cosmo, z, c.mdef)
    R_new = RΔ(c,cosmo,z,c.mdef)
    c_new = R_new ./ c.rs
    M_new = R_to_M.(R_new, Ref(rho_thresh))
    return NFW(c=cosmo,M=M_new, conc=c_new, z=z, mdef=c.mdef)
end

################################################################################################################
#### functions based on densities
"""
    ρ(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
Returns the density of the halo at radius r in Msun/kpc^3. 
"""
function ρ(r::Union{Real,u.Quantity},c::NFW)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    if isscalar(c.rs) && isscalar(c.rhos)
        x = r / c.rs
        return c.rhos / x / (1.0 + x)^2
    else
        result=zeros( Base.promote_eltype(c.rs,c.rhos,[1.0]), length(c.rs) )
        for i in eachindex(c.rs,c.rhos)
            x = r/c.rs[i]
            result[i] = c.rhos[i] / x / (1 + x)^2
        end
        return result
    end
end
ρ(r::AbstractArray,c::AbstractHaloProfile) = [ρ(i,c) for i in r]

"""
    ∇ρ(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
Returns the gradient density of the halo at radius r. 
"""
function ∇ρ(r::Union{Real,u.Quantity},c::NFW)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    if isscalar(c.rs) && isscalar(c.rhos)
        return -c.rs^3 * c.rhos * (3 * r + c.rs) / (r^2 * (r + c.rs)^3)
    else
        result=zeros( Base.promote_eltype(c.rs,c.rhos,[1.0]), length(c.rs) )
        for i in eachindex(c.rs,c.rhos)
            result[i] = -c.rs[i]^3 * c.rhos[i] * (3 * r + c.rs[i]) /
                (r^2 * (r + c.rs[i])^3)
        end
        return result
    end
end
∇ρ(r::Union{Real,u.Quantity}, c::AbstractHaloProfile) = c_fdm(x->ρ(x,c), r)
∇ρ(r::AbstractArray,c::AbstractHaloProfile) = [∇ρ(i,c) for i in r]

"""
    dlogρ_dlogr(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
Returns the dimensionless logarithmic derivative of the halo density at radius r. Specifically, this is dln(ρ)/dln(r) = dρ/dr / (r * ρ(r))
"""
function dlogρ_dlogr(r::Union{Real,u.Quantity},c::NFW)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    if isscalar(c.rs) && isscalar(c.rhos)
        x = r / c.rs
        return -(1 + 2 * x / (1 + x) )
    else
        result=zeros( Base.promote_eltype(c.rs,c.rhos,[1.0]), length(c.rs) )
        for i in eachindex(c.rs,c.rhos)
            x = r/c.rs[i]
            result[i] = -(1 + 2 * x / (1 + x) )
        end
        return result
    end
end
dlogρ_dlogr(r::Union{Real,u.Quantity},c::AbstractHaloProfile) = (r isa u.Quantity && (r=r|>ua.kpc|>u.ustrip); ∇ρ(r,c) * r / ρ(r,c) )
dlogρ_dlogr(r::AbstractArray,c::AbstractHaloProfile) = [dlogρ_dlogr(i,c) for i in r]

"""
    ρmean(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
Returns the mean density inside r. Equivalent to the enclosed mass divided by the volume. """
function ρmean(r::Union{Real,u.Quantity},c::NFW)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    if isscalar(c.rs) && isscalar(c.rhos)
        return 3 * c.rs^3 * c.rhos * (-r + ( r + c.rs ) * log( ( r + c.rs ) /
            c.rs ) ) / (r^3 * ( r + c.rs ) )
    else
        result=zeros( Base.promote_eltype(c.rs,c.rhos,[1.0]), length(c.rs) )
        for i in eachindex(c.rs,c.rhos)
            result[i] = 3 * c.rs[i]^3 * c.rhos[i] * (-r + ( r + c.rs[i] ) * log( ( r + c.rs[i] ) /
                c.rs[i] ) ) / (r^3 * ( r + c.rs[i] ) )
        end
        return result
    end
end
ρmean(r::Union{Real,u.Quantity},c::AbstractHaloProfile; kws...) = enclosed_mass(r,c; kws...) / (4.0 / 3.0 * π * r^3)
ρmean(r::AbstractArray,c::AbstractHaloProfile) = [ρmean(i,c) for i in r]

"""
    ∇ρmean(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
The gradient of the average density inside R."""
function ∇ρmean(r::Union{Real,u.Quantity},c::NFW)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    if isscalar(c.rs) && isscalar(c.rhos)
        return -3 * c.rs^3 * c.rhos * (-r*(4*r+3*c.rs) + 3*log((r+c.rs)/c.rs)*(r+c.rs)^2) /
		(r^4 * (r+c.rs)^2)
    else
        result=zeros( Base.promote_eltype(c.rs,c.rhos,[1.0]), length(c.rs) )
        for i in eachindex(c.rs,c.rhos)
            result[i] = -3 * c.rs[i]^3 * c.rhos[i] * (-r*(4*r+3*c.rs[i]) + 3*log((r+c.rs[i])/c.rs[i])*
                (r+c.rs[i])^2) / (r^4 * (r+c.rs[i])^2)
        end
        return result
    end
end
∇ρmean(r::Union{Real,u.Quantity},c::AbstractHaloProfile) = c_fdm(x->ρmean(x,c),r)
∇ρmean(r::AbstractArray,c::AbstractHaloProfile) = [∇ρmean(i,c) for i in r]

"""
    enclosed_mass(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
The total mass enclosed within r in Msun."""
function enclosed_mass(r::Union{Real,u.Quantity},c::NFW)
    if isscalar(c.rs) && isscalar(c.rhos)
        return 4 * π * c.rs^3 * c.rhos * NFWmu(r|>u.unit(c.rs),c.rs)
    else
        return [4 * π * c.rs[i]^3 * c.rhos[i] * NFWmu(r,c.rs[i]) for i in eachindex(c.rs,c.rhos)]
    end
end
enclosed_mass(r::Union{Real,u.Quantity},c::AbstractHaloProfile; kws...) = 4 * π * quadgk(x-> ρ(x,c) * x^2, 0*u.unit(r), r; kws...)[1]
enclosed_mass(r::AbstractArray,c::AbstractHaloProfile) = [enclosed_mass(i,c) for i in r]

"""
    Σ(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
The projected surface density at radius r in Msun / kpc^2. Uses the analytic formula of Lokas & Mamon 2001 for NFW. """
function Σ(r::Union{Real,u.Quantity},c::NFW)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    if isscalar(c.rs) && isscalar(c.rhos)
        x = r / c.rs
        if x ≈ 1
            return c.rhos * c.rs * 2 / 3
        elseif x<1
            x2m1 = x^2 - 1
            return c.rs * c.rhos * 2 / x2m1 * (1 - 2 / sqrt(-x2m1) * atanh( sqrt( (1 - x) / (x + 1) ) ) )
        else#if x>1
            x2m1 = x^2 - 1
            return c.rs * c.rhos * 2 / x2m1 * (1 - 2 / sqrt(x2m1) * atan( sqrt( (x - 1) / (x + 1) ) ) )
        end
    else
        result=zeros( Base.promote_eltype(c.rs,c.rhos,[1.0]), length(c.rs) )
        for i in eachindex(c.rs,c.rhos)
            x = r / c.rs[i]
            if x ≈ 1
                result[i] = c.rhos[i] * c.rs[i] * 2 / 3
            elseif x<1
                x2m1 = x^2 - 1
                result[i] = c.rs[i] * c.rhos[i] * 2 / x2m1 * (1 - 2 / sqrt(-x2m1) * atanh( sqrt( (1 - x) / (x + 1) ) ) )
            else#if x>1
                x2m1 = x^2 - 1
                result[i] = c.rs[i] * c.rhos[i] * 2 / x2m1 * (1 - 2 / sqrt(x2m1) * atan( sqrt( (x - 1) / (x + 1) ) ) )
            end
        end
        return result
    end
end
Σ(r::Union{Real,u.Quantity}, c::AbstractHaloProfile; kws...) = (r isa u.Quantity && (r=r|>ua.kpc|>u.ustrip); 2 * quadgk(x -> exp(x)^2 * ρ(exp(x),c) / sqrt(exp(x)^2 - r^2), log(r), log(1e6); maxevals=500, kws...)[1]) #poorly optimized, should overwrite
Σ(r::AbstractArray,c::AbstractHaloProfile) = [Σ(i,c) for i in r]

"""
    projected_mass(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
The total projected mass interior to radius r (kpc or Unitful) in Msun. Uses the analytical form of Lokas & Mamon 2001 (Equation 43) for NFW.
"""
function projected_mass(r::Union{Real,u.Quantity},c::NFW)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    if isscalar(c.rs) && isscalar(c.rhos)
        x = r / c.rs
        if x ≈ 1
            return 4 * π * c.rhos * c.rs^3 * ( x + log( x / 2 ) )
            # return 4 * π * c.rhos * c.rs^3 * (1 - log(2) )
        elseif x<1
            return 4 * π * c.rhos * c.rs^3 * ( acosh( 1 / x ) / sqrt( 1 - x^2 ) + log( x / 2 ) )
        else#if x>1
            return 4 * π * c.rhos * c.rs^3 * ( acos( 1 / x ) / sqrt( x^2 - 1 ) + log( x / 2 ) )
        end
    else
        result=zeros( Base.promote_eltype(c.rs,c.rhos,[1.0]), length(c.rs) )
        for i in eachindex(c.rs,c.rhos)
            x = r / c.rs[i]
            if x ≈ 1
                result[i] = 4 * π * c.rhos[i]* c.rs[i]^3 * ( x + log( x / 2 ) )
            elseif x<1
                result[i] = 4 * π * c.rhos[i] * c.rs[i]^3 * ( acosh( 1 / x ) / sqrt( 1 - x^2 ) + log( x / 2 ) )
            else#if x>1
                result[i] = 4 * π * c.rhos[i] * c.rs[i]^3 * ( acos( 1 / x ) / sqrt( x^2 - 1 ) + log( x / 2 ) )
            end
        end
        return result
    end
end
projected_mass(r::Union{Real,u.Quantity}, c::AbstractHaloProfile; kws...) = 2π * quadgk(x -> exp(x)^2 * Σ(exp(x),c), 0, log(r); maxevals=1000, kws...)[1] #poorly optimized, should overwrite
projected_mass(r::AbstractArray,c::AbstractHaloProfile) = [projected_mass(i,c) for i in r]

# """
#     Φ(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
# The gravitational potential. If no units are provided, assumes r and rs in kpc and rhos in solMass/kpc^3, and returns (km/s)^2."""
# function Φ(r::Union{Real,u.Quantity},c::NFW)
#     if isscalar(c.rs) && isscalar(c.rhos)
#         typeof(r) <: u.Quantity && (r = r|> ua.kpc |> u.ustrip)
#         typeof(c.rs) <: u.Quantity ? rs = c.rs |> ua.kpc |> u.ustrip : rs=c.rs
#         typeof(c.rhos) <: u.Quantity ? rhos = c.rhos |> (ua.Msun/ua.kpc^3) |> u.ustrip : rhos=c.rhos
#         return 4 * π * (4.30091727003628e-6) * rs^3 * rhos * ( log( rs ) - log( r + rs ) ) / r * (u.km/u.s)^2
#         # constant is G in km^2 kpc s^-2 M⊙^-1
#     else
#         tt = promote_type(typeof(r|>u.ustrip), typeof(c.rs[1]|>u.ustrip), typeof(c.rhos[1]|>u.ustrip), Float64)
#         result = zeros(typeof(tt(1.0)*(u.km/u.s)^2), length(c.rhos))
#         typeof(r) <: u.Quantity && (r = r|> ua.kpc |> u.ustrip)
#         eltype(c.rs) <: u.Quantity ? rs_conv=true : rs_conv=false
#         eltype(c.rhos) <: u.Quantity ? rhos_conv=true : rhos_conv=false
#         for i in eachindex(c.rs,c.rhos)
#             rs_conv ? (rs=c.rs[i] |> ua.kpc |> u.ustrip) : rs=c.rs[i]
#             rhos_conv ? (rhos=c.rhos[i] |> (ua.Msun/ua.kpc^3) |> u.ustrip) : rhos=c.rhos[i]
#             result[i] = 4 * π * (4.30091727003628e-6) * rs^3 * rhos *
#                 ( log( rs ) - log( r + rs ) ) / r * (u.km/u.s)^2
#         end
#         return result
#     end
# end
"""
    Φ(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
The gravitational potential. If no units are provided, assumes r and rs in kpc and rhos in solMass/kpc^3, and returns (km/s)^2."""
function Φ(r::Union{Real,u.Quantity},c::NFW)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    if isscalar(c.rs) && isscalar(c.rhos)
        return 4 * π * G * c.rs^3 * c.rhos * ( log( c.rs ) - log( r + c.rs ) ) / r
    else
        result=zeros( Base.promote_eltype(c.rs,c.rhos,[1.0]), length(c.rs) )
        for i in eachindex(c.rs,c.rhos)
            result[i] = 4 * π * G * c.rs[i]^3 * c.rhos[i] * ( log( c.rs[i] ) - log( r + c.rs[i] ) ) / r
        end
        return result
    end
end
Φ(r::AbstractArray,c::AbstractHaloProfile) = [Φ(i,c) for i in r]

"""
    ∇Φ(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
The gradient of the gravitational potential. If no units are provided, assumes r and rs in kpc and rhos in solMass/kpc^3, and returns km / s^2."""
function ∇Φ(r::Union{Real,u.Quantity},c::NFW)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    if isscalar(c.rs) && isscalar(c.rhos)
        return 4. * π * G * c.rs^3 * c.rhos / KPC_KM * ( log( r + c.rs ) - log( c.rs ) - r / (r+c.rs) ) / r^2
    else
        result=zeros( Base.promote_eltype(c.rs,c.rhos,[1.0]), length(c.rs) )
        for i in eachindex(c.rs,c.rhos)
            result[i] = 4. * π * G * c.rs[i]^3 * c.rhos[i] / KPC_KM * ( log( r + c.rs[i] ) - log( c.rs[i] ) - r / (r+c.rs[i]) ) / r^2
        end
        return result
    end
end
∇Φ(r::AbstractArray,c::AbstractHaloProfile) = [∇Φ(i,c) for i in r]

"""
    ∇∇Φ(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
The second gradient of the gravitational potential. If no units are provided, assumes r and rs in kpc and rhos in solMass/kpc^3, and returns s^-2."""
function ∇∇Φ(r::Union{Real,u.Quantity},c::NFW)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    if isscalar(c.rs) && isscalar(c.rhos)
        return ( 4 * π * G * c.rs^3 * c.rhos / KPC_KM^2 / (r * (r + c.rs)^2) ) - 
		( 8 * G * π * c.rs^3 * c.rhos / KPC_KM^2 * ( log(r + c.rs) - log(c.rs) 
		- r/(r+c.rs) ) / r^3 )
    else
        result=zeros( Base.promote_eltype(c.rs,c.rhos,[1.0]), length(c.rs) )
        for i in eachindex(c.rs,c.rhos)
            result[i] = ( 4 * π * G * c.rs[i]^3 * c.rhos[i] / KPC_KM^2 / (r * (r + c.rs[i])^2) ) - 
		( 8 * G * π * c.rs[i]^3 * c.rhos[i] / KPC_KM^2 * ( log(r + c.rs[i]) - log(c.rs[i]) 
		- r/(r+c.rs[i]) ) / r^3 )
        end
        return result
    end
end
∇∇Φ(r::AbstractArray,c::AbstractHaloProfile) = [∇∇Φ(i,c) for i in r]

"""
    ∇²Φ(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
The Laplacian (divergence of the gradient) of the gravitational potential. This is expressly equal to 4 * π * G * ρ. If no units are provided, assumes r and rs in kpc and rhos in solMass/kpc^3, and returns s^-2."""
function ∇²Φ(r::Union{Real,u.Quantity},c::AbstractHaloProfile)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    dens=ρ(r,c)
    return @. 4 * π * G / KPC_KM^2 * dens
end
∇²Φ(r::AbstractArray,c::AbstractHaloProfile) = [∇²Φ(i,c) for i in r]

"""
    circular_velocity(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
The circular velocity of the halo at radius r. If no units are provided, assumes r and rs in kpc and rhos in solMass/kpc^3, and returns km/s."""
function circular_velocity(r::Union{Real,u.Quantity},c::AbstractHaloProfile)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    M_enc = enclosed_mass(r,c)
    return @. sqrt( G * M_enc / r)
end
circular_velocity(r::AbstractArray,c::AbstractHaloProfile) = [circular_velocity(i,c) for i in r]

"""
    escape_velocity(r::Union{Real,u.Quantity,AbstractArray{<:Real},AbstractArray{<:u.Quantity}},c::AbstractHaloProfile)
The escape velocity of the halo at radius r. If no units are provided, assumes r and rs in kpc and rhos in solMass/kpc^3, and returns km/s."""
function escape_velocity(r::Union{Real,u.Quantity},c::AbstractHaloProfile)
    eltype(r) <: u.Quantity && (r = r .|> ua.kpc .|> u.ustrip)
    pot = Φ(r,c)
    return @. sqrt( 2 * -pot )
end
escape_velocity(r::AbstractArray,c::AbstractHaloProfile) = [escape_velocity(i,c) for i in r]

"""
    V_max(c::NFW)
Returns (rmax, vmax); the radius in kpc of maximum circular velocity, and the corresponding velocity in km/s. """
V_max(c::NFW) = (rmax = 2.16258 .* c.rs; vmax = circular_velocity.(rmax,Ref(c)); return rmax, vmax )

########################################################################################
########################################################################################
########################################################################################
#### extra


"""
The dynamical time of a halo.
	
The dynamical time can be estimated in multiple ways, but is almost always based on the ratio
of a distance to circular velocity. The relevant distance can be defined in different ways as
indicated with the ``definition`` parameter. The dynamical time is more succinctly expressed 
as a multiple of the Hubble time which depends on the overdensity threshold and redshift.
	
Parameters
-----------------------------------------------------------------------------------------------
c::AbstractCosmology
z: array_like
	Redshift; can be a number or a numpy array.
mdef: str
	The mass definition. 
definition: str
	An identifier for a definition of the dynamical time. Valid definitions are ``crossing``
	(the crossing time), ``peri`` (the time to reach the halo center, half the crossing time)
	and ``orbit`` (the time to orbit around the halo, crossing time times π).
		
Returns
-----------------------------------------------------------------------------------------------
t_dyn: array_like
	Dynamical time in Gyr; has the same dimensions as ``z``. """
function dynamical_time(c::AbstractCosmology, z, mdef, definition = "crossing")
    t_cross = 2^1.5 * hubble_time(c, z) * (density_threshold(c, z, mdef) / ρ_c(ua.Msun/ua.kpc^3,c,z))^(-0.5)
    if definition == "crossing"
	t_dyn = t_cross
    elseif definition == "peri"
	t_dyn = t_cross / 2.0
    elseif definition == "orbit"
	t_dyn = t_cross * π
    else
	throw(DomainError(definition,"Invalid definition."))
    end
    return t_dyn
end



end
