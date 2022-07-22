module cosmo

import Unitful as u
import UnitfulAstro as ua
import UnitfulEquivalences as ue
import Dierckx: Spline1D, derivative
import OrdinaryDiffEq: ODEProblem, solve, Tsit5
import SpecialFunctions: gamma
import Roots: find_zero, Bisection
import QuadGK: quadgk
import Trapz: trapz
import DoubleExponentialFormulas: quaddeo, quadde
using PyCall
# import Trapz: trapz

__precompile__(true)

const camb = PyNULL()
function __init__()
    copy!(camb, pyimport_conda("camb","camb","conda-forge"))
end

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
    n_nu,
    T_nu, T_cmb,
    # nu_relative_density,
    # de_density_scale,
    z_at_value,
    Ω_m, Ω_dm, Ω_b, Ω_k, Ω_γ, Ω_ν, Ω_r, Ω_Λ, # Ω_sum,
    w,
    sound_horizon,
    setup_growth, growth_integral,
    transfer_function, Sugiyama1995,Eisenstein1998,Eisenstein1998_ZB,Eisenstein1999,CAMB,TopHat,Gaussian,SharpK,
    power_spectrum,
    σ2, filter_function_k, filter_function_r, dw_dlnkr, dw_dkr,
    dlnσ2_dlnr, dlnσ_dlnr, dlnr_dlnm, dlnσ2_dlnm, dlnσ_dlnm,
    ξmm, # dξmm_dr,
    δc, νpeak2,
    concentration,
    massfunc_dndm
    # NFW, ρ, ∇ρ, ρmean, ∇ρmean, enclosedMass, Φ, ∇Φ, Δvir # halo is now a separate module so it should be imported separately. 


include("utils.jl")
include("constants.jl")
import ..constants
# include("types_base.jl")
# include("standard_functions.jl")
# include("default_cosmologies.jl")
# these two are included in types_base.jl
# include("transfer.jl")
# include("growth.jl")

#### type definitions #########################################################################################
abstract type AbstractCosmology end
Base.Broadcast.broadcastable(m::AbstractCosmology) = Ref(m)
# we could do more efficient, automatic broadcasting with StructArrays, but already hand-coded the iterables which is perhaps regrettable. Could do a refactoring if necessary. 
abstract type AbstractClosedCosmology <: AbstractCosmology end
abstract type AbstractFlatCosmology <: AbstractCosmology end
abstract type AbstractOpenCosmology <: AbstractCosmology end

include("transfer.jl") # include here because it requires AbstractCosmology, but provides the filter types needed to define the concrete cosmologies

#############################################################################
struct FlatLCDM{T <: Real, S <: Union{T,Array{<:T}}} <: AbstractFlatCosmology
# struct FlatLCDM{T <: Real} <: AbstractFlatCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_dm::T
    Ω_b::T
    Ω_r::T
    Ω_γ::T
    Ω_n::T
    Tcmb0::T # in kelvin
    Neff::T
    m_nu::S # Union{T,Array{<:T}} # in eV
    growth_function::Function
    growth_function_inv::Function
    # vel_div::Spline1D
    ns::T
    sigma8::T
    z_eq::T
    k::Vector{<:T}
    Pk::Vector{<:T}
    Pk_interp::Spline1D
    Pk_normalize::T
    σ2_interp::Function
    σ2_interp_inv::Function
    filt::Type{<:PkFilter}
end
# FlatLCDM(h::Real, Ω_Λ::Real, Ω_m::Real, Ω_dm::Real, Ω_b::Real, Ω_r::Real, Tcmb0::Real, Neff::Real,
#          m_nu::Union{Real, Array{<:Real}}) =
    # FlatLCDM(float(h), float(Ω_Λ), float(Ω_m), float(Ω_dm), float(Ω_b),
             #                  float(Ω_r), float(Tcmb0), float(Neff), float(m_nu))
        # try
        #     FlatLCDM(promote(float(h), float(Ω_Λ), float(Ω_m), float(Ω_dm), float(Ω_b),
        #                      float(Ω_r), float(Tcmb0), float(Neff), float(m_nu))...)

        # catch x
        #     FlatLCDM(float(h), float(Ω_Λ), float(Ω_m), float(Ω_dm), float(Ω_b),
        #              float(Ω_r), float(Tcmb0), float(Neff), float(m_nu))
             # end


#############################################################################
struct ClosedLCDM{T <: Real, S <: Union{T,Array{<:T}}} <: AbstractClosedCosmology
# struct ClosedLCDM{T <: Real} <: AbstractClosedCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_dm::T
    Ω_b::T
    Ω_r::T
    Ω_γ::T
    Ω_n::T
    Tcmb0::T
    Neff::T
    m_nu::S # Union{T,Array{<:T}} # in eV
    growth_function::Function
    growth_function_inv::Function
    # vel_div::Spline1D
    ns::T
    sigma8::T
    z_eq::T
    k::Vector{<:T}
    Pk::Vector{<:T}
    Pk_interp::Spline1D
    Pk_normalize::T
    σ2_interp::Function
    σ2_interp_inv::Function
    filt::Type{<:PkFilter}
end

#############################################################################
struct OpenLCDM{T <: Real, S <: Union{T,Array{<:T}}} <: AbstractOpenCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_dm::T
    Ω_b::T
    Ω_r::T
    Ω_γ::T
    Ω_n::T
    Tcmb0::T
    Neff::T
    m_nu::S # Union{T,Array{<:T}} # in eV
    growth_function::Function
    growth_function_inv::Function
    # vel_div::Spline1D
    ns::T
    sigma8::T
    z_eq::T
    k::Vector{<:T}
    Pk::Vector{<:T}
    Pk_interp::Spline1D
    Pk_normalize::T
    σ2_interp::Function
    σ2_interp_inv::Function
    filt::Type{<:PkFilter}
end

#############################################################################
# the same cosmologies but with w0-wa-state dark energy
for c in ("Flat", "Open", "Closed")
    name = Symbol("$(c)WCDM")
    @eval begin
        struct $(name){TTT <: Real, SSS <: Union{Real,Array{TTT}}} <: $(Symbol("Abstract$(c)Cosmology"))
        # struct $(name){TTT <: Real} <: $(Symbol("Abstract$(c)Cosmology"))
            h::TTT
            Ω_k::TTT
            Ω_Λ::TTT
            Ω_m::TTT
            Ω_dm::TTT
            Ω_b::TTT
            Ω_r::TTT
            Ω_γ::TTT
            Ω_n::TTT
            Tcmb0::TTT
            Neff::TTT
            m_nu::SSS # Union{TTT,Array{<:TTT}} # in eV
            growth_function::Function
            growth_function_inv::Function
            # vel_div::Spline1D
            ns::TTT
            sigma8::TTT
            z_eq::TTT
            k::Vector{<:TTT}
            Pk::Vector{<:TTT}
            Pk_interp::Spline1D
            Pk_normalize::TTT
            σ2_interp::Function
            σ2_interp_inv::Function
            filt::Type{<:PkFilter}
            w0::TTT
            wa::TTT
        end
    end
end
#############################################################################

function WCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_dm::Real, Ω_b::Real, Ω_r::Real, Ω_γ::Real, Ω_n::Real,
              Tcmb0::Real, Neff::Real, m_nu::Union{Real, Array{<:Real}}, growth_function::Function,
              growth_function_inv::Function, ns::Real, sigma8::Real, z_eq::Real,
              k::Union{Real, Vector{<:Real}}, Pk::Union{Real, Vector{<:Real}},
              Pk_interp::Spline1D, Pk_normalize::Real, σ2_interp::Function, σ2_interp_inv::Function, filt::Type{<:PkFilter}, w0::Real, wa::Real)
    if Ω_k < 0
        ClosedWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_dm, Ω_b, Ω_r, Ω_γ, Ω_n, Tcmb0, Neff, m_nu, growth_function, growth_function_inv,
                   ns, sigma8, z_eq, k, Pk, Pk_interp, Pk_normalize, σ2_loginterp, σ2_loginterp_inv, filt, w0, wa)
    elseif Ω_k > 0
        OpenWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_dm, Ω_b, Ω_r, Ω_γ, Ω_n, Tcmb0, Neff, m_nu, growth_function, growth_function_inv,
                 ns, sigma8, z_eq, k, Pk, Pk_interp, Pk_normalize, σ2_loginterp, σ2_loginterp_inv, filt, w0, wa)
    else
        FlatWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_dm, Ω_b, Ω_r, Ω_γ, Ω_n, Tcmb0, Neff, m_nu, growth_function, growth_function_inv,
                 ns, sigma8, z_eq, k, Pk, Pk_interp, Pk_normalize, σ2_loginterp, σ2_loginterp_inv, filt, w0, wa)
    end
end

#############################################################################

include("growth.jl")
include("standard_functions.jl")

#############################################################################

"""
    cosmology(;h = 0.6766, OmegaK = 0, OmegaM = 0.30966, OmegaR = nothing, Tcmb0 = 2.7255, w0 = -1, wa = 0, Neff = 3.046, m_nu=[0., 0., 0.06] * u.eV, ns=0.9665, sigma8=0.8102, tau=0.0561, z_reion=7.82, a_min=1e-5, transfer_model=DefaultTransferModel, k_min=1e-20, k_max=1e20, dlogk=0.01,filt::Type{<:PkFilter}=DefaultFilter)



# Parameters
* `h` - Dimensionless Hubble constant
* `OmegaK` - Curvature density (Ω_k)
* `OmegaM` - Matter density (Ω_m)
* `OmegaR` - Radiation density (Ω_r) - includes photons Ω_γ and massless neutrinos
* `Tcmb0` - CMB temperature in Kelvin; used to compute Ω_γ
* `OmegaG` - Photon density (Ω_γ)
* `w0` - CPL dark energy equation of state; `w = w0 + wa(1-a)`
* `wa` - CPL dark energy equation of state; `w = w0 + wa(1-a)`
* `Neff` - Effective number of massless neutrino species; used to compute Ω_ν
* `m_nu` - Neutrino masses 
* `ns` - The tilt of the primordial power spectrum.
* `sigma8` - The normalization of the power spectrum, i.e. the variance when the field is filtered with a top hat filter of radius 8 Mpc/h. 
* `tau`
* `z_reion` - The redshift of reionization. 
* `a_min` - Minimum scale factor at which to solve the ODE for the growth factor.
* `transfer_model` - string specifying which of the transfer function models to use.
* `k_min` - Minimum wavenumber k in 1/Mpc at which to calculate the power spectrum
* `k_max` - Minimum wavenumber k in 1/Mpc at which to calculate the power spectrum
* `dlnk` - spacing of k in natural log k

# Examples
```jldoctest
julia> c = cosmology()
Cosmology.FlatLCDM{Float64}(0.69, 0.7099122024007928, 0.29, 8.77975992071536e-5)

julia> c = cosmology(OmegaK=0.1)
Cosmology.OpenLCDM{Float64}(0.69, 0.1, 0.6099122024007929, 0.29, 8.77975992071536e-5)

julia> c = cosmology(w0=-0.9, OmegaK=-0.1)
Cosmology.ClosedWCDM{Float64}(0.69, -0.1, 0.8099122024007929, 0.29, 8.77975992071536e-5, -0.9, 0.0)
```
"""
function cosmology(;h = 0.6766,         # Scalar; Hubble constant at z = 0 / 100 [km/sec/Mpc]
                   OmegaK = 0,          # Energy density of curvature in units of the critical density at z=0. 
                   OmegaM = 0.30966,    # Energy density of matter in units of the critical density at z=0. 
                   OmegaB = 0.04897,    # Energy density of matter in units of the critical density at z=0. 
                   OmegaR = nothing,    # Energy density of radiation in units of the critical density at z=0. 
                   Tcmb0 = 2.7255*u.K,  # Temperature of CMB at z=0 in Kelvin; can have units
                   OmegaG = nothing,    # Energy density in radiation; if nothing, computed from neutrinos and Tcmb. 
                   w0 = -1,             
                   wa = 0,
                   Neff = 3.046,
                   m_nu=[0, 0, 0.06] * u.eV,
                   ns=0.9665,
                   sigma8=0.8102,
                   tau=0.0561,
                   z_reion=7.82,
                   a_min=1e-5,            # minimum scale factor for the growth factor calculation.
                   # transfer_model="eisenstein98",    # model to compute the transfer function
                   transfer_model::Type{<:TransferModel}=DefaultTransferModel, # model to compute the transfer function, defined in transfer.py
                   k_min=1e-20,           # minimum wavenumber [h/Mpc] for calculation of power spectrum
                   k_max=1e20,            # maximum wavenumber [h/Mpc] for calculation of power spectrum
                   dlogk=0.01,            # spacing in log10(k)
                   filt::Type{<:PkFilter}=DefaultFilter,    # filter to use for computing the σ2 interp, defined in transfer.jl
                   camb_nonlinear_model=nothing)            # string specifying the halofit version to use; mead2020 is usually good. 
                   
    # value checks; not completely done. 
    if OmegaM < 0
        throw(DomainError(OmegaM,"OmegaM must be greater than or equal to 0."))
    end
    if OmegaB < 0
        throw(DomainError(OmegaB,"OmegaB must be greater than or equal to 0."))
    end
    if Tcmb0 |> u.ustrip < 0
        throw(DomainError(Tcmb0,"Tcmb0 must be greater than or equal to 0."))
    end
    if OmegaB < 0
        throw(DomainError(OmegaB,"OmegaB must be greater than or equal to 0."))
    end
    if OmegaB > OmegaM
        throw(DomainError(OmemaB,"OmegaB cannot be greater than OmegaM."))
    end
    if Neff < 0
        throw(DomainError(Neff,"Effective number of neutrino species cannot be negative."))
    end

    # compute Ω_dm from OmegaM and OmegaB
    OmegaDM = OmegaM - OmegaB
    
    if Tcmb0 isa u.Quantity #check if Tcmb has a unit
        Tcmb0=Tcmb0 |> u.K |> u.ustrip
    end
    # else
    #     Tcmb0*=u.K
    # end

    # for simplicity, pass m_nu on in eV but without units.
    if m_nu !=0 && m_nu !== nothing
        try #if mass units, convert with UnitfulEquivalences
            m_nu = u.ustrip.(u.uconvert.(u.eV, m_nu, ue.MassEnergy()))
        catch xxx #else just convert energy units
            m_nu = m_nu .|> u.eV .|> u.ustrip
        end
    end

    if OmegaG === nothing && (Tcmb0 != 0 && Tcmb0 !==nothing) #compute Omega due to photons (cmb)
        # this is the full equation, but the units and constants are a bit slow. 
        # OmegaG = (4 * σ / c_0^ 3 |> u.g / u.cm^3 / u.K^4) * Tcmb0^4 / 
        #     (3.0 / (8.0 * π * G) * (h * 100 * u.km / u.s / ua.Mpc)^2 |> u.g/u.cm^3)
        OmegaG = 4.481620089297254e-7 * Tcmb0^4 / h^2
    elseif OmegaG !== nothing && (Tcmb0==0 || Tcmb0===nothing)
        Tcmb0 = (OmegaG / 4.481620089297254e-7 * h^2)^(1/4)
    end
    
    if Tcmb0 !== nothing && Neff !== nothing
        if m_nu===nothing; m_nu=0.; end
        Tnu0 = 0.7137658555036082 * Tcmb0
        n_nu = Int(floor(Neff)) # number of neutrino species
        prefac = 0.22710731766023898 #7/8 * (4/11)^(4/3)
        if 0 in (m_nu .== 0) #if at least one neutrino has mass
            p = 1.83
            invp = 0.54644808743  # 1.0 / p
            k = 0.3173
            nu_y = m_nu ./ (8.617333262145179e-5 * Tnu0)  #the constant is k_B in ev / K
            rel_mass_per = @. (1 + (k * nu_y[nu_y!=0])^p)^invp
            rel_mass = sum(rel_mass_per) + length(m_nu[m_nu.==0])
            nu_rel_density0 = prefac * rel_mass * (Neff / n_nu )
        else #if all neutrinos are massless
            nu_rel_density0 = prefac * Neff
        end
        OmegaN = OmegaG * nu_rel_density0
    end
    if OmegaR===nothing
        OmegaR = OmegaG + OmegaN
    end
    
    # if OmegaR === nothing
    #     OmegaG = 4.48131e-7 * Tcmb0^4 / h^2
    #     OmegaN = Neff * OmegaG * (7 / 8) * (4 / 11)^(4 / 3)
    #     OmegaR = OmegaG + OmegaN
    # end

    OmegaL = 1 - OmegaK - OmegaM - OmegaR # calculate the dark energy density
    z_eq = 3600 * (OmegaM * h^2 / 0.15) - 1 #  calculate redshift of matter-radiation equality
    ## type conversions ###############################################################################################
    # for each component that can be a vector, test if it is and if so copy out its first element to test promotion
    if length(m_nu)>1
        m_nu_test = deepcopy(m_nu[1])
    else
        m_nu_test = deepcopy(m_nu)
    end
    # doing type conversions; do all the scalars here
    h, OmegaK, OmegaM, OmegaB, OmegaR, OmegaL, OmegaDM, Tcmb0, OmegaG, OmegaN, w0, wa, Neff, ns, sigma8, tau, z_reion, m_nu_test, k_min, k_max, dlogk, z_eq =
        promote(h, OmegaK, OmegaM, OmegaB, OmegaR, OmegaL, OmegaDM, Tcmb0|>u.ustrip, OmegaG, OmegaN, w0, wa, Neff, ns, sigma8, tau, z_reion, m_nu_test, k_min, k_max, dlogk, z_eq)
    # now do the (possible) vectors, since promote doesn't really play nice with them.
    m_nu=typeof(h).(m_nu)

    ### initialize the precomputed functions ##########################################################################
    # growth_function, vel_div, growth_function_inv = setup_growth(OmegaM,OmegaL,OmegaK,OmegaR,a_min)
    growth_function_tmp, growth_function_inv_tmp = setup_growth(OmegaM,OmegaL,OmegaK,OmegaG,Neff,m_nu,Tcmb0,w0,wa,a_min)
    # the above are interpolated in scale factor, convert them to accept redshifts. 
    growth_function(z) = growth_function_tmp((1 ./ (1 .+ z)))
    growth_function_inv(g) = 1 ./ growth_function_inv_tmp(g) .- 1
    # vel_div=Spline1D(1:10,rand(10)) # placeholder till API settles

    if transfer_model===CAMB
        k_max>1e3 && println("When using CAMB for power spectra, the runtime scales strongly with k_max. Your k_max of $k_max may take a long time to run or crash.")
        camb_params=camb.CAMBparams(H0=h*100,ombh2=OmegaB*h^2,omch2=OmegaDM*h^2, omk=OmegaK, nnu=Neff, standard_neutrino_neff=Neff, num_massive_neutrinos=length(m_nu[m_nu.!=0]), m_nu=sum(m_nu),TCMB=Tcmb0,InitPower=camb.initialpower.InitialPowerLaw(ns=ns),WantTransfer=true, Transfer=camb.model.TransferParams(high_precision=true, accurate_massive_neutrinos=true,kmin=k_min,kmax=k_max,PK_num_redshifts=1,PK_redshifts=[0.],npoints=Int((log10(k_max)-log10(k_min))÷dlogk)),WantCMB=false,WantScalars=false,WantDerivedParameters=false,Want_cl_2D_array=false,Want_CMB_lensing=false,DoLensing=false)
        camb_nonlinear_model != nothing && (camb_params.NonLinear="NonLinear_pk"; camb_params.NonLinearModel.set_params(halofit_version=camb_nonlinear_model))
        camb_result = camb.get_results(camb_params)
        camb_sigma8 = camb_result.get_sigma8()[1]
        k, _, unnormalized_pk = camb_result.get_matter_power_spectrum(minkh=k_min,maxkh=k_max,npoints=Int((log10(k_max)-log10(k_min))÷dlogk))
        unnormalized_pk=vec(unnormalized_pk)
        normalized_pk = @. unnormalized_pk * sigma8^2 / camb_sigma8^2 # renormalize to correct sigma8
        Pk_normalize=camb_result.get_matter_transfer_data().transfer_data[7,:][1] * sigma8^2 / camb_sigma8^2 / (2π)
        # not sure this normalization is correct, but it's somewhat close and not likely to be used
        # Pk_normalize=promote_type(typeof(Pk_normalize),typeof(h))(Pk_normalize)
    else
        k = logspace(k_min, k_max, step=dlogk)
        unnormalized_T = transfer_function(k, h, OmegaM, OmegaB, OmegaN, Neff, Tcmb0, transfer_model)
        # unnormalized_T /= maximum(unnormalized_T)
        unnormalized_pk = power_spectrum(k,unnormalized_T, ns)
        sigma2_8_norm = σ2(8.0, Spline1D(k,unnormalized_pk; k=4, bc="error", s=0.0); growth_factor=1.0,
                       k_min=minimum(k), k_max=maximum(k), j=0, filt = TopHat)
        Pk_normalize = sigma8^2 / sigma2_8_norm
        normalized_pk = unnormalized_pk * Pk_normalize
    end
    pk_interp = Spline1D(k,normalized_pk; k=4, bc="zero", s=0.0)
    sigma2_r = logspace(1/k_max,1/k_min,length=100)
    s2_arr = σ2(sigma2_r, pk_interp; growth_factor=1.0, k_min=k_min, k_max=k_max, j=0, filt=filt, maxevals=1000)
    sigma2_interp_tmp = Spline1D(log10.(sigma2_r), log10.(s2_arr); k=4, bc="zero",s=0.0) # this is the slowest part so far but it is not TERRIBLY slow and having a sigma2 interpolator is very useful for other functions.
    sigma2_interp(x) = exp10.(sigma2_interp_tmp(log10.(x)))
    sigma2_interp_inv_tmp = Spline1D(log10.(reverse(s2_arr)), log10.(reverse(sigma2_r)); k=4, bc="error",s=0.0) # inverse of above
    sigma2_interp_inv(x) = exp10.(sigma2_interp_inv_tmp(log10.(x)))
    # generally I would prefer to use bc="error" for the interpolation, but quaddeo requires the ability to evaluate from 0 to infinity, unfortunately. 
    # pk_interp(k) = power_spectrum(k,
    #                               transfer_function(k, h, OmegaM, OmegaB, OmegaN, Neff, Tcmb0, transfer_model),
    #                               ns, 1.0) * Pk_normalize
                   
    ### initializing structs ##########################################################################################
    if !(w0 == -1 && wa == 0)
        return WCDM(h, OmegaK, OmegaL, OmegaM, OmegaDM, OmegaB, OmegaR, OmegaG, OmegaN, Tcmb0, Neff, m_nu,
                    growth_function, growth_function_inv, ns, sigma8, z_eq, k, normalized_pk, pk_interp, Pk_normalize, sigma2_interp, sigma2_interp_inv, filt, w0, wa)
    end
    if OmegaK < 0
        return ClosedLCDM(h, OmegaK, OmegaL, OmegaM, OmegaDM, OmegaB, OmegaR, OmegaG, OmegaN, Tcmb0,
                          Neff, m_nu, growth_function, growth_function_inv, ns, sigma8, z_eq, k, normalized_pk, pk_interp, Pk_normalize, sigma2_interp, sigma2_interp_inv, filt)
    elseif OmegaK > 0
        return OpenLCDM(h, OmegaK, OmegaL, OmegaM, OmegaDM, OmegaB, OmegaR, OmegaG, OmegaN, Tcmb0, Neff, m_nu,
                        growth_function, growth_function_inv, ns, sigma8, z_eq, k, normalized_pk, pk_interp, Pk_normalize, sigma2_interp, sigma2_interp_inv, filt)
    else
        return FlatLCDM(h, OmegaK, OmegaL, OmegaM, OmegaDM, OmegaB, OmegaR, OmegaG, OmegaN, Tcmb0, Neff, m_nu,
                        growth_function, growth_function_inv, ns, sigma8, z_eq, k, normalized_pk, pk_interp, Pk_normalize, sigma2_interp, sigma2_interp_inv, filt)
    end
end

# h, OmegaK, OmegaM, OmegaB, OmegaR, Tcmb0, OmegaG, w0, wa, Neff, m_nu, n, sigma8, tau, z_reion = 0.6766, 0, 0.30966, 0.04897, nothing, 2.7255*u.K, nothing, -1, 0, 3.046, [0, 0, 0.06] * u.eV, 0.9665, 0.8102, 0.0561, 7.82

include("default_cosmologies.jl")
include("halo.jl")
using .halo
include("peaks.jl")
using .peaks
include("concentration.jl")
using .conc_module
include("hmf.jl")
using .hmf
include("nonlinear.jl")
using .nonlinear

end # module
