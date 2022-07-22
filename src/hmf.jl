module hmf

import UnitfulAstro as ua
import Unitful as u
import SpecialFunctions: gamma, erf
import Dierckx: Spline1D
import ..cosmo: σ2, dlnσ_dlnm, dlnσ2_dlnm, ρ_m, Ω_m, AbstractCosmology, filter_m_to_r, isscalar, cumtrapz, DefaultFilter, PkFilter, SharpK, logspace
import ..peaks: lagrangianR, νpeak2, δc

export massfunc_dndm, fσ, Press1972, Sheth1999, Sheth2001, Jenkins2001, Reed2003, Warren2006, Reed2007, Peacock2007, Courtin2010, Crocce2010, Pillepich2010, Bhattacharya2011, Angulo2012_FoF, Angulo2012_Bound, Watson2013_FoF, Ishiyama2015

abstract type AbstractHMF end
Base.Broadcast.broadcastable(m::AbstractHMF) = Ref(m)

abstract type Press1972 <: AbstractHMF end
abstract type Sheth1999 <: AbstractHMF end
abstract type Sheth2001 <: AbstractHMF end
abstract type Jenkins2001 <: AbstractHMF end
abstract type Reed2003 <: AbstractHMF end
abstract type Warren2006 <: AbstractHMF end
abstract type Reed2007 <: AbstractHMF end
abstract type Peacock2007 <: AbstractHMF end
abstract type Courtin2010 <: AbstractHMF end
abstract type Crocce2010 <: AbstractHMF end
abstract type Pillepich2010 <: AbstractHMF end
abstract type Bhattacharya2011 <: AbstractHMF end
abstract type Angulo2012_FoF <: AbstractHMF end
abstract type Angulo2012_Bound <: AbstractHMF end
abstract type Watson2013_FoF <: AbstractHMF end
abstract type Ishiyama2015 <: AbstractHMF end

# wdm
abstract type Benson2013 <: AbstractHMF end

const DefaultHMF = Sheth2001 # the default HMF model to use.

#### fitting functions

fσ(cc::Type{Press1972},nu2) = @. sqrt(2/π)*sqrt(nu2)*exp(-0.5*nu2)
fσ(cc::Type{Sheth1999},nu2;a=0.707, p=0.3, A=0.3222) = (nu=sqrt.(nu2); nuprime=nu .* a; @. A * (1 + (1/nuprime^p)) * sqrt(nuprime/2) * exp(-nuprime^2/2) / sqrt(π))
fσ(cc::Type{Sheth2001},nu2;a=0.707, p=0.3, A=0.3222) = @. A * sqrt(2*a/π) * sqrt(nu2) * exp(-(a*nu2/2)) * (1 + (1/(a*nu2))^p)
fσ(cc::Type{Courtin2010},nu2;a=0.695,p=0.1,A=0.348) = fσ(Sheth2001,nu2;a=a,p=p,A=A)
fσ(cc::Type{Bhattacharya2011},nu2,z::Real;A_a=0.333,A_b=0.11,a_a=0.788,a_b=0.01,p=0.807,q=1.795) = (A = A_a*(1+z)^-A_b; a=a_a*(1+z)^-a_b; νfν = fσ(Sheth2001,nu2;A=A,a=a,p=p); @. νfν * (sqrt(a) * sqrt(nu2))^(q-1))
# Jenkins is much lower than all the others. Think it is dN / dlnσ^-1 ; might need to multiply by dln(σ^-1)_dlnm. Confusing. Think there is also a term for the power spectrum slope. 
fσ(cc::Type{Jenkins2001},sigma2;A=0.315,b=0.61,c=3.8) = @. A * exp(-abs(-log(sigma2)/2 + b)^c)
fσ(cc::Type{Warren2006},sigma2;A=0.7234,b=1.625,c=0.2538,d=1.1982,e=1) = @. A*((e/sqrt(sigma2))^b + c) * exp(-d/sigma2)
fσ(cc::Type{Watson2013_FoF},sigma2;A=0.282,b=2.163,c=1,d=1.21,e=1.406) = fσ(Warren2006,sigma2;A=A,b=b,c=c,d=d,e=e)
fσ(cc::Type{Crocce2010},sigma2,z::Real;A_a=0.58,A_b=0.13,b_a=1.37,b_b=0.15,c_a=0.3,c_b=0.084,d_a=1.036,d_b=0.024,e=1) = (A=A_a*(1+z)^-A_b; b=b_a*(1+z)^-b_b; c=c_a*(1+z)^-c_b; d=d_a*(1+z)^-d_b; fσ(Warren2006,sigma2;A=A,b=b,c=c,d=d,e=e) )
fσ(cc::Type{Pillepich2010},sigma2;A=0.6853,b=1.868,c=0.3324,d=1.2266,e=1) = fσ(Warren2006,sigma2;A=A,b=b,c=c,d=d,e=e)
fσ(cc::Type{Ishiyama2015},sigma2;A=0.193,b=1.550,c=1,d=1.186,e=2.184) = fσ(Warren2006,sigma2;A=A,b=b,c=c,d=d,e=e)
fσ(cc::Type{Reed2003},nu2,sigma2;a=0.707,p=0.3,A=0.3222,c=0.7) = (@assert length(nu2)==length(sigma2); νfν = fσ(Sheth2001,nu2;a=a, p=p, A=A); sigma=sqrt.(sigma2); @. νfν * exp(-c/(sigma*cosh(2*sigma)^5)))
fσ(cc::Type{Reed2007},nu2::Real,sigma2::Real,n_eff::Real; A=0.3222, p=0.3, c=1.08, a=0.764 ) = (lnsigma = -log(sigma2) / 2; nu = sqrt(nu2); G_1 = exp(-(lnsigma-0.4)^2 / 0.72); G_2 = exp(-(lnsigma-0.75)^2 / 0.08); return A * sqrt(2a/π) * (1 + (1/(a*nu2))^p + 0.6*G_1 + 0.4*G_2) * nu * exp(-c*a*nu2/2 - 0.03*nu^0.6/(n_eff + 3)^2))
function fσ(cc::Type{Reed2007},nu2,sigma2,n_eff; A=0.3222, p=0.3, c=1.08, a=0.764)
    a=a/c
    @assert length(nu2) == length(sigma2) == length(n_eff)
    result = similar(nu2)
    @inbounds for i in eachindex(result)
        lnsigma = -log(sigma2[i]) / 2
        nu = sqrt(nu2[i])
        G_1 = exp(-(lnsigma-0.4)^2 / 0.72)
        G_2 = exp(-(lnsigma-0.75)^2 / 0.08)
        result[i] = A * sqrt(2a/π) * (1 + (1/(a*nu2[i]))^p + 0.6*G_1 + 0.4*G_2) * nu * exp(-c*a*nu2[i]/2 - 0.03*nu^0.6/(n_eff[i] + 3)^2)
    end
    return result
end
fσ(cc::Type{Peacock2007},nu2::Real;a=1.529,b=0.704,c=0.412) = (nu = sqrt(nu2); d = 1 + a*nu^b; return nu * exp(-c * nu2) * (2 * c * d * nu + b * a * nu^(b-1)) / d^2)
function fσ(cc::Type{Peacock2007},nu2;a=1.529,b=0.704,c=0.412) 
    result = similar(nu2)
    @inbounds for i in eachindex(nu2)
        nu = sqrt(nu2[i])
        d = 1 + a*nu^b
        result[i] = nu * exp(-c * nu2[i]) * (2 * c * d * nu + b * a * nu^(b-1)) / d^2
    end
    return result
end
fσ(cc::Type{Angulo2012_FoF},sigma2;A=0.201,b=1.7,c=1.172,d=2.08) = @. A * ((d/sqrt(sigma2))^b+1) * exp(-c/sigma2)
fσ(cc::Type{Angulo2012_Bound},sigma2;A=0.265,b=1.9,c=1.4,d=1.675) = fσ(Angulo2012_FoF,sigma2;A=A,b=b,c=c,d=d)

# function fσ_interpolator(cc::Type{Benson2013}, sigma2::AbstractVector, barrier::AbstractVector)
#     @assert length(sigma2) == length(barrier)
#     result = similar(sigma2)
#     sigma2_diff = diff(sigma2)
#     idx = eachindex(sigma2,barrier)
#     len = length(sigma2)
#     for i in idx
#         sigma2[i] == 0 && (result[i] = 0; continue)
#         i == 1 && (result[i] = 2/sigma2_diff[i] * (1 - erf(barrier[i]/sqrt(2*sigma2[i]))); continue)
#         i == len && (result[i] = 2 / sigma2_diff[i-1] * (1 - erf(barrier[i] / sqrt(2 * sigma2[i])) - 2 * sum( (1 - erf( (barrier[i] - barrier[j]) / sqrt( sigma2[i] - sigma2[j]))) * result[j] for j in 1:len-1)); continue)
                     
#         tmp_idx = view(idx,1:i-1)
#         result[i] = 2 / sigma2_diff[i-1] * (1 - erf(barrier[i] / sqrt(2 * sigma2[i])) - sum( (1 - erf( (barrier[i] - barrier[j]) / sqrt(2 * (sigma2[i]-sigma2[j])))) * result[j] * (sigma2_diff[j] + sigma2_diff[j+1]) / 2 for j in tmp_idx) )
#     end
#     # return Spline1D(sigma2,result; k=4, bc="error", s=0.0)
#     return result
# end

function fσ_benson2013(c::AbstractCosmology, k::AbstractVector, Pk::AbstractVector, mx::Real, gx::Real, z::Real=0.0, Mmin::Real=1e5, Mmax::Real=1e16; filt::Type{<:PkFilter}=SharpK)
    rho_mean = ρ_m(ua.Msun / ua.Mpc^3, c, z) |> u.ustrip
    g_z = c.growth_function(z)
    M_arr = reverse(logspace(Mmin, Mmax, length=300))
    R_arr = filter_m_to_r.(M_arr,rho_mean,filt)
    R_min, R_max = extrema(R_arr)# minimum(R_arr), maximum(R_arr)
    sigma2_arr = σ2(R_arr,k,Pk; growth_factor=g_z, filt=filt)
    sigma2_min, sigma2_max = extrema(sigma2_arr)
    sigma2_M = Spline1D(sigma2_arr,M_arr; k=4, bc="error", s=0.0)

    sigma2 = range(sigma2_min,sigma2_max,length=1000)
    dsigma2 = sigma2.step.hi
    delta_c = δc(c,sigma2_M(sigma2),z,mx,gx)
    barrier = @. delta_c * sqrt(0.707) * (1 + 0.5 * (sigma2 / (0.707 * delta_c^2))^0.6) # / 1.197
    result = similar(sigma2)
    idx = eachindex(sigma2)
    for i in idx
        sigma2[i] == 0 && (result[i] = 0; continue)
        i == 1 && (result[i] = 2 / sigma2[i] * (1 - erf(barrier[i]/sqrt(2*sigma2[i]))); continue)
        tmp_idx = view(idx,1:i-1)
        result[i] = 2 / dsigma2 * (1 - erf(barrier[i]/sqrt(2*sigma2[i]))) - 2 * sum( (1 - erf( (barrier[i] - barrier[j]) / sqrt(2 * (sigma2[i] - sigma2[j])))) * result[j] for j in tmp_idx)
    end
    M_arr = sigma2_M(sigma2)
    R_arr = filter_m_to_r.(M_arr,rho_mean,filt)
    return R_arr, M_arr, sigma2, result
end

####################################
### more basic functions
""" 
    press1972_fσ(nu2)
Evaluates the multiplicity function from Press & Schecter et al. 1972. """
press1972_fσ(nu2) = @. sqrt(2/π)*sqrt(nu2)*exp(-0.5*nu2)
"""
    sheth2001_fσ(nu2;a=0.707, p=0.3, A=0.3222)
Evaluates the multiplicity function from Sheth et al. 2001. a and p are fitting coefficients, and A is calculated to normalize it correctly.  """
sheth2001_fσ(nu2;a=0.707, p=0.3, A=0.3222) =
    # (A=1/(1 + (2^-p)*gamma(0.5 - p)/gamma(0.5)); # this is a normalization factor but it is harmful for courtin2010; dont really need it
    @. A * sqrt(2*a/π) * sqrt(nu2) * exp(-(a*nu2/2)) * (1 + (1/(a*nu2))^p)
courtin2010_fσ(nu2;a=0.695,p=0.1,A=0.348) = sheth2001_fσ(nu2;a=a,p=p,A=A)
bhattacharya2011_fσ(nu2,z::Real;A_a=0.333,A_b=0.11,a_a=0.788,a_b=0.01,p=0.807,q=1.795) = (A = A_a*(1+z)^-A_b; a=a_a*(1+z)^-a_b; νfν = sheth2001_fσ(nu2;A=A,a=a,p=p); @. νfν * (sqrt(a) * sqrt(nu2))^(q-1))
"""
    jenkins2001_fσ(sigma2;A=0.315,b=0.61,c=3.8)
Evaluates the multiplicity function from Jenkins 2001. A, b, and c are fitting coefficients. This function takes sigma2 rather than nu2. """
jenkins2001_fσ(sigma2;A=0.315,b=0.61,c=3.8) = @. A * exp(-abs(-log(sigma2)/2 + b)^c)
"""
    warren2006_fσ(sigma2;A=0.7234,b=1.625,c=0.2538,d=1.1982,e=1)
Evaluates the multiplicity function from Warren et al. 2006 (http://adsabs.harvard.edu/abs/2006ApJ...646..881W). This function takes sigma2 rather than nu2.
"""
warren2006_fσ(sigma2;A=0.7234,b=1.625,c=0.2538,d=1.1982,e=1) = @. A*((e/sqrt(sigma2))^b + c) * exp(-d/sigma2)
"""
    watson2013_fof_fσ(sigma2;A=0.282,b=2.163,c=1,d=1.21,e=1.406)
    watson2013_fof_fσ(cc::HMFinit;A=0.282,b=2.163,c=1,d=1.21,e=1.406)
The friends-of-friends HMF derived in Watson et al. 2013 https://ui.adsabs.harvard.edu/abs/2013MNRAS.433.1230W/abstract. Uses WMAP5 cosmology.  """
watson2013_fof_fσ(sigma2;A=0.282,b=2.163,c=1,d=1.21,e=1.406) = warren2006_fσ(sigma2;A=A,b=b,c=c,d=d,e=e)
crocce2010_fσ(sigma2,z::Real;A_a=0.58,A_b=0.13,b_a=1.37,b_b=0.15,c_a=0.3,c_b=0.084,d_a=1.036,d_b=0.024,e=1) = (A=A_a*(1+z)^-A_b; b=b_a*(1+z)^-b_b; c=c_a*(1+z)^-c_b; d=d_a*(1+z)^-d_b; warren2006_fσ(sigma2;A=A,b=b,c=c,d=d,e=e) )
pillepich2010_fσ(sigma2;A=0.6853,b=1.868,c=0.3324,d=1.2266,e=1) = warren2006_fσ(sigma2;A=A,b=b,c=c,d=d,e=e)
ishiyama2015_fσ(sigma2;A=0.193,b=1.550,c=1,d=1.186,e=2.184) = warren2006_fσ(sigma2;A=A,b=b,c=c,d=d,e=e)
reed2003_fσ(nu2,sigma2;a=0.707,p=0.3,A=0.3222,c=0.7) = (νfν = sheth2001_fσ(nu2;a=a, p=p, A=A); sigma=sqrt.(sigma2); @. νfν * exp(-c/(sigma*cosh(2*sigma)^5)))

# function reed2007_fσ(cc::HMFinit; A=0.3222, p=0.3, c=1.08, a=0.764)
#     a=a/c
#     G_1 = @. exp(-(cc.lnsigma-0.4)^2 / 0.72)
#     G_2 = @. exp(-(cc.lnsigma-0.75)^2 / 0.08)
#     return @. A * sqrt(2a/π) * (1 + (1/(a*cc.nu2))^p + 0.6*G_1 + 0.4*G_2) * cc.nu *
#                exp(-c*a*cc.nu2/2 - 0.03*cc.nu^0.6/(cc.n_eff + 3)^2)
# end

peacock2007_fσ(nu2;a=1.529,b=0.704,c=0.412) = (nu=sqrt.(nu2); d = @. 1 + a*nu^b; @. nu * exp(-c * nu2) * (2 * c * d * nu + b * a * nu^(b-1)) / d^2)

angulo2012_fof_fσ(sigma2;A=0.201,b=1.7,c=1.172,d=2.08) = @. A * ((d/sqrt(sigma2))^b+1) * exp(-c/sigma2)

angulo2012_bound_fσ(sigma2;A=0.265,b=1.9,c=1.4,d=1.675) = angulo2012_fof_fσ(sigma2;A=A,b=b,c=c,d=d)

###############################################################################################################
#### high level functions

dispatch_fσ(model::Union{Type{Press1972},Type{Sheth1999},Type{Sheth2001},Type{Courtin2010},Type{Peacock2007}},nu2,sigma2,dlns_dlnm,z) = fσ(model,nu2)
dispatch_fσ(model::Type{Bhattacharya2011},nu2,sigma2,dlns_dlnm,z) = fσ(model,nu2,z)
dispatch_fσ(model::Union{Type{Jenkins2001},Type{Warren2006},Type{Watson2013_FoF},Type{Pillepich2010},Type{Ishiyama2015},Type{Angulo2012_FoF},Type{Angulo2012_Bound}},nu2,sigma2,dlns_dlnm,z) = fσ(model,sigma2)
dispatch_fσ(model::Type{Crocce2010},nu2,sigma2,dlns_dlnm,z) = fσ(model,sigma2,z)
dispatch_fσ(model::Type{Reed2003},nu2,sigma2,dlns_dlnm,z) = fσ(model,nu2,sigma2)
dispatch_fσ(model::Type{Reed2007},nu2,sigma2,dlns_dlnm,z) = (n_eff = @. -3 * (2 * dlns_dlnm + 1); fσ(model,nu2,sigma2,n_eff))

""" 
    massfunc_dndm(M::Union{Real,AbstractArray},c::AbstractCosmology,model::Type{T}=DefaultHMF,z::Real=0.0;filt=nothing, kws...)
Calculate the halo mass function (dN/dM) in units `h^4 Msun^-1 Mpc^-3` for the masses given in `M_arr` at redshift `z` using the cosmology and interpolators in `c` using the power spectrum filter `filt` under the halo mass function model `model`. Valid `model` are: Press1972, Sheth1999, Sheth2001, Jenkins2001, Reed2003, Warren2006, Reed2007, Peacock2007, Courtin2010, Crocce2010, Pillepich2010, Bhattacharya2011, Angulo2012\\_FoF, Angulo2012\\_Bound, Watson2013\\_FoF, Ishiyama2015. `kws...` are options to be passed to `quadgk` for the `σ2` integrations. """
function massfunc_dndm(M::Union{Real,AbstractArray},c::AbstractCosmology,model::Type{T}=DefaultHMF,z::Real=0.0;filt=nothing, kws...) where T<:AbstractHMF
    filt === nothing && (filt=c.filt)
    rho_mean0 = u.ustrip(ρ_m(ua.Msun/ua.Mpc^3,c,0.0)) / c.h^2 # if we divide by c.h^2 here, we get the same result as from python HMF.
    R = filter_m_to_r.(M,rho_mean0,filt)
    sigma2 = σ2(R, c, z; filt=filt, kws...)
    nu2 = δc(c,z)^2 ./ sigma2
    dlns_dlnm = dlnσ_dlnm(R,c,z; kws...)
    fsigma = dispatch_fσ(model,nu2,sigma2,dlns_dlnm,z)
    dndm = @. rho_mean0/M^2 * fsigma * abs(dlns_dlnm)
    # dn_dlnm = M .* dndm
    # dn_dlog10m = log(10) .* dn_dlnm
    # nless = cumtrapz(
    return dndm
end
"""
    massfunc_dndm(M::Union{Real,AbstractArray}, k::AbstractVector, Pk::AbstractVector, h::Real, Ωm::Real, Ωk::Real, ΩΛ::Real, Ωb::Real, Tcmb0::Union{Real,u.Quantity}, ns::Real, Neff::Real, m_nu, w0, wa, model::Type{T}=DefaultHMF, z::Real=0.0; filt::PkFilter=DefaultFilter)
Calculate the halo mass function (dN/dM) in units `h^4 Msun^-1 Mpc^-3` using explicitly passed cosmological parameters.
"""
function massfunc_dndm(M::Union{Real,AbstractArray}, k::AbstractVector, Pk::AbstractVector, h::Real, Ωm::Real, Ωk::Real, ΩΛ::Real, Ωb::Real, Tcmb0::Union{Real,u.Quantity}, ns::Real, Neff::Real, m_nu, w0, wa, model::Type{T}=DefaultHMF, z::Real=0.0; filt::PkFilter=DefaultFilter) where T<:AbstractHMF
    # set up constants and growth factor
    rho_mean0 = u.ustrip(ρ_m(ua.Msun/ua.Mpc^3,0.0,h,Ωm)) / h^2 # if we divide by c.h^2 here, we get the same result as from python HMF.
    Tcmb0 isa u.Quantity && (Tcmb0 = u.ustrip(Tcmb0|>u.K))
    Ωγ = 4.481620089297254e-7 * Tcmb0^4 / h^2 # calculate photon density from Tcmb0
    # set up growth function
    growth_function_tmp, growth_function_inv_tmp = setup_growth(Ωm, ΩΛ, Ωk, Ωγ, Neff, m_nu, Tcmb0, w0, wa) 
    growth_function(z) = growth_function_tmp((1 ./ (1 .+ z)))
    growth_function_inv(g) = 1 ./ growth_function_inv_tmp(g) .- 1
    g_z = growth_function(z)
    delta_c = δc(z,h,Ωm,Ωk,ΩΛ,Tcmb0,m_nu,Neff,w0,wa)
    
    R = filter_m_to_r.(M,rho_mean0,filt)
    sigma2 = σ2(R, k, Pk; filt=filt, growth_factor = g_z)
    nu2 = delta_c^2 ./ sigma2
    dlns_dlnm = dlnσ_dlnm(R,k,Pk; filt=filt, growth_factor = g_z)
    fsigma = dispatch_fσ(model,nu2,sigma2,dlns_dlnm,z)
    dndm = @. rho_mean0/M^2 * fsigma * abs(dlns_dlnm)
    return dndm
end
"""
    massfunc_dndm(M::Union{Real,AbstractArray}, c::AbstractCosmology, k::AbstractVector, Pk::AbstractVector, model::Type{T}=DefaultHMF, z::Real=0.0; filt=nothing)
Calculate the halo mass function (dN/dM) in units `h^4 Msun^-1 Mpc^-3` using the cosmological parameters in `c` but using the power spectrum `Pk` sampled at wavenumbers `k`.
"""
function massfunc_dndm(M::Union{Real,AbstractArray}, c::AbstractCosmology, k::AbstractVector, Pk::AbstractVector, model::Type{T}=DefaultHMF, z::Real=0.0; filt=nothing) where T<:AbstractHMF
    filt === nothing && (filt=c.filt)
    rho_mean0 = u.ustrip(ρ_m(ua.Msun/ua.Mpc^3,c,0.0)) / c.h^2 # if we divide by c.h^2 here, we get the same result as from python HMF.
    g_z = c.growth_function(z)
    R = filter_m_to_r.(M,rho_mean0,filt)
    sigma2 = σ2(R, k, Pk; filt=filt, growth_factor=g_z)
    nu2 = δc(c,z)^2 ./ sigma2
    dlns_dlnm = dlnσ_dlnm(R,k,Pk; filt=filt, growth_factor=g_z)
    fsigma = dispatch_fσ(model,nu2,sigma2,dlns_dlnm,z)
    dndm = @. rho_mean0/M^2 * fsigma * abs(dlns_dlnm)
    return dndm
end
"""
    massfunc_dndm(M::Union{Real,AbstractArray}, c::AbstractCosmology, k::AbstractVector, Pk::AbstractVector, mx::Real, gx::Real=1.5, model::Type{T}=DefaultHMF, z::Real=0.0; filt::PkFilter=SharpK)
For WDM using the formalism of Benson 2013. Calculate the halo mass function (dN/dM) in units `h^4 Msun^-1 Mpc^-3` using the cosmological parameters in `c` but using the power spectrum `Pk` sampled at wavenumbers `k`. `mx` is the mass of the WDM particle in keV, and `gx` is the effective degrees of freedom with 1.5 being the expected value for a fermionic spin-1/2 particle.
"""
function massfunc_dndm(M::Union{Real,AbstractArray}, c::AbstractCosmology, k::AbstractVector, Pk::AbstractVector, mx::Real, gx::Real=1.5, model::Type{Benson2013}=Benson2013, z::Real=0.0; filt::Type{<:PkFilter}=SharpK)
    rho_mean0 = u.ustrip(ρ_m(ua.Msun/ua.Mpc^3,c,0.0)) / c.h^2 # if we divide by c.h^2 here, we get the same result as from python HMF.
    g_z = c.growth_function(z)
    delta_c = δc(c,M,z,mx,gx)
    
    R = filter_m_to_r.(M,rho_mean0,filt)
    dlns_dlnm = dlnσ_dlnm(R,k,Pk; filt=filt, growth_factor=g_z)
    sigma2 = σ2(R,k,Pk; filt=filt, growth_factor = g_z)

    R_,M_,sigma2_,fsigma_ = fσ_benson2013(c, k, Pk, mx, gx, z; filt=filt)
    interp = Spline1D(log10.(reverse(M_)),log10.(reverse(fsigma_)); k=1, bc="error", s=0.0)
    fsigma = exp10.(interp(log10.(M)))
    
    dndm = @. rho_mean0/M^2 * fsigma * abs(dlns_dlnm) * 2 * (delta_c / sqrt(sigma2))^-1.8 * 2.5
    return dndm
end


# cumtrapz(reverse(py_mf.m),reverse(massfunc_dndm(py_mf.m,a,"SMT")))
# cumtrapz(-reverse(py_mf.m),reverse(massfunc_dndm(py_mf.m,a,"SMT")))
# reverse(cumtrapz(-reverse(py_mf.m),reverse(massfunc_dndm(py_mf.m,a,"SMT"))))
# reverse(cumtrapz(-reverse(log.(py_mf.m)),reverse(py_mf.m.*massfunc_dndm(py_mf.m,a,"SMT"))))
end #module
