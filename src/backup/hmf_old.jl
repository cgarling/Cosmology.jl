module hmf

import UnitfulAstro as ua
import Unitful as u
import SpecialFunctions: gamma
# using StructArrays
import ..cosmo: σ2, dlnσ_dlnm, ρ_m, Ω_m, AbstractCosmology, filter_m_to_r, isscalar, cumtrapz
import ..peaks: lagrangianR, νpeak2, δc

export massfunc_dndm

# abstract type AbstractHMF end
# Base.Broadcast.broadcastable(m::AbstractHMF) = Ref(m)

struct HMFinit{T <: Union{Real,Array{<:Real}}}
    sigma::T
    sigma2::T
    lnsigma::T
    dlns_dlnm::T
    n_eff::T
    nu::T
    nu2::T
    z::Real
    Omz::Real
end
Base.Broadcast.broadcastable(m::HMFinit) = Ref(m)

# still working on this interface for the result of massfunc. 
# struct HMFresult{T <: Union{Real,Array{<:Real}}}
#     # from HMFinit
#     sigma::T
#     sigma2::T
#     lnsigma::T
#     dlns_dlnm::T
#     n_eff::T
#     nu::T
#     nu2::T
#     z::Real
#     Omz::Real
#     # other variables and results
#     M::T          # masses at which the HMF was evaluated
#     dndm::T       # The differential mass function in units of h^4 / Msun / Mpc^3
#     dn_dlnm::T    # The differential mass function in natural log of mass, in units of h^3 / Mpc^3
#     dn_dlog10m::T # The differential mass function in log10 of mass, in units of h^3 / Mpc^3
#     # cmf::T      # cumulative mass function
#     ngreater::T   # the number of halos above M in h^3 / Mpc^3
#     nless::T      # the number of halos below M in h^3 / Mpc^3
#     ρ_greater::T  # the mass density in halos with masses above M in Msun h^2 / Mpc^3
#     ρ_less::T     # the mass density in halos with masses below M in Msun h^2 / Mpc^3
#     box_size::T   # Box side length in which to expect one halo of mass M in Mpc/h 

# end

#### fitting functions

""" 
    press1972_fσ(nu2)
Evaluates the multiplicity function from Press & Schecter et al. 1972. """
press1972_fσ(nu2) = @. sqrt(2/π)*sqrt(nu2)*exp(-0.5*nu2)
press1972_fσ(c::HMFinit) = press1972_fσ(c.nu2)

"""
    sheth2001_fσ(nu2;a=0.707, p=0.3, A=0.3222)
Evaluates the multiplicity function from Sheth et al. 2001. a and p are fitting coefficients, and A is calculated to normalize it correctly.  """
sheth2001_fσ(nu2;a=0.707, p=0.3, A=0.3222) =
    # (A=1/(1 + (2^-p)*gamma(0.5 - p)/gamma(0.5)); # this is a normalization factor but it is harmful for courtin2010; dont really need it
    @. A * sqrt(2*a/π) * sqrt(nu2) * exp(-(a*nu2/2)) * (1 + (1/(a*nu2))^p)
sheth2001_fσ(c::HMFinit;a=0.707, p=0.3, A=0.3222) = sheth2001_fσ(c.nu2;a=a,p=p,A=A)

courtin2010_fσ(nu2;a=0.695,p=0.1,A=0.348) = sheth2001_fσ(nu2;a=a,p=p,A=A)
courtin2010_fσ(c::HMFinit;a=0.695,p=0.1,A=0.348) = sheth2001_fσ(c;a=a,p=p,A=A)

bhattacharya2011_fσ(nu2,z::Real;A_a=0.333,A_b=0.11,a_a=0.788,a_b=0.01,p=0.807,q=1.795) = (A = A_a*(1+z)^-A_b; a=a_a*(1+z)^-a_b; νfν = sheth2001_fσ(nu2;A=A,a=a,p=p); @. νfν * (sqrt(a) * sqrt(nu2))^(q-1))
bhattacharya2011_fσ(c::HMFinit;A_a=0.333,A_b=0.11,a_a=0.788,a_b=0.01,p=0.807,q=1.795) = (A = A_a*(1+c.z)^-A_b; a=a_a*(1+c.z)^-a_b; νfν = sheth2001_fσ(c;A=A,a=a,p=p); @. νfν * (sqrt(a) * c.nu)^(q-1))

"""
    jenkins2001_fσ(sigma2;A=0.315,b=0.61,c=3.8)
Evaluates the multiplicity function from Jenkins 2001. A, b, and c are fitting coefficients. This function takes sigma2 rather than nu2. """
jenkins2001_fσ(sigma2;A=0.315,b=0.61,c=3.8) = @. A * exp(-abs(-log(sigma2)/2 + b)^c)
jenkins2001_fσ(cc::HMFinit) = jenkins2001_fσ(cc.sigma2)
"""
    warren2006_fσ(sigma2;A=0.7234,b=1.625,c=0.2538,d=1.1982,e=1)
Evaluates the multiplicity function from Warren et al. 2006 (http://adsabs.harvard.edu/abs/2006ApJ...646..881W). This function takes sigma2 rather than nu2.
"""
warren2006_fσ(sigma2;A=0.7234,b=1.625,c=0.2538,d=1.1982,e=1) = @. A*((e/sqrt(sigma2))^b + c) * exp(-d/sigma2)
warren2006_fσ(cc::HMFinit;A=0.7234,b=1.625,c=0.2538,d=1.1982,e=1) = @.  A*((e/cc.sigma)^b + c) * exp(-d/cc.sigma2)

"""
    watson2013_fof_fσ(sigma2;A=0.282,b=2.163,c=1,d=1.21,e=1.406)
    watson2013_fof_fσ(cc::HMFinit;A=0.282,b=2.163,c=1,d=1.21,e=1.406)
The friends-of-friends HMF derived in Watson et al. 2013 https://ui.adsabs.harvard.edu/abs/2013MNRAS.433.1230W/abstract. Uses WMAP5 cosmology.  """
watson2013_fof_fσ(sigma2;A=0.282,b=2.163,c=1,d=1.21,e=1.406) = warren2006_fσ(sigma2;A=A,b=b,c=c,d=d,e=e)
watson2013_fof_fσ(cc::HMFinit;A=0.282,b=2.163,c=1,d=1.21,e=1.406) = warren2006_fσ(cc;A=A,b=b,c=c,d=d,e=e)

crocce2010_fσ(sigma2,z::Real;A_a=0.58,A_b=0.13,b_a=1.37,b_b=0.15,c_a=0.3,c_b=0.084,d_a=1.036,d_b=0.024,e=1) = (A=A_a*(1+z)^-A_b; b=b_a*(1+z)^-b_b; c=c_a*(1+z)^-c_b; d=d_a*(1+z)^-d_b; warren2006_fσ(sigma2;A=A,b=b,c=c,d=d,e=e) )
crocce2010_fσ(cc::HMFinit;A_a=0.58,A_b=0.13,b_a=1.37,b_b=0.15,c_a=0.3,c_b=0.084,d_a=1.036,d_b=0.024,e=1) = (A=A_a*(1+cc.z)^-A_b; b=b_a*(1+cc.z)^-b_b; c=c_a*(1+cc.z)^-c_b; d=d_a*(1+cc.z)^-d_b; warren2006_fσ(cc;A=A,b=b,c=c,d=d,e=e) )

pillepich2010_fσ(sigma2;A=0.6853,b=1.868,c=0.3324,d=1.2266,e=1) = warren2006_fσ(sigma2;A=A,b=b,c=c,d=d,e=e)
pillepich2010_fσ(cc::HMFinit;A=0.6853,b=1.868,c=0.3324,d=1.2266,e=1) = warren2006_fσ(cc;A=A,b=b,c=c,d=d,e=e)

ishiyama2015_fσ(sigma2;A=0.193,b=1.550,c=1,d=1.186,e=2.184) = warren2006_fσ(sigma2;A=A,b=b,c=c,d=d,e=e)
ishiyama2015_fσ(cc::HMFinit;A=0.193,b=1.550,c=1,d=1.186,e=2.184) = warren2006_fσ(cc;A=A,b=b,c=c,d=d,e=e)

reed2003_fσ(nu2,sigma2;a=0.707,p=0.3,A=0.3222,c=0.7) = (νfν = sheth2001_fσ(nu2;a=a, p=p, A=A); sigma=sqrt.(sigma2); @. νfν * exp(-c/(sigma*cosh(2*sigma)^5)))
reed2003_fσ(cc::HMFinit;a=0.707,p=0.3,A=0.3222,c=0.7) = (νfν = sheth2001_fσ(cc.nu2;a=a, p=p, A=A); @. νfν * exp(-c/(cc.sigma*cosh(2*cc.sigma)^5)))

function reed2007_fσ(cc::HMFinit; A=0.3222, p=0.3, c=1.08, a=0.764)
    a=a/c
    G_1 = @. exp(-(cc.lnsigma-0.4)^2 / 0.72)
    G_2 = @. exp(-(cc.lnsigma-0.75)^2 / 0.08)
    return @. A * sqrt(2a/π) * (1 + (1/(a*cc.nu2))^p + 0.6*G_1 + 0.4*G_2) * cc.nu *
               exp(-c*a*cc.nu2/2 - 0.03*cc.nu^0.6/(cc.n_eff + 3)^2)
    # this reduces allocations but doesnt really help speed at all. 
    # result = similar(cc.nu2)
    # a=a/c
    # @inbounds for i in eachindex(result)
    #     G_1 = exp(-(cc.lnsigma[i]-0.4)^2 / 0.72)
    #     G_2 = exp(-(cc.lnsigma[i]-0.75)^2 / 0.08)
    #     result[i] = A * sqrt(2a/π) * (1 + (1/(a*cc.nu2[i]))^p + 0.6*G_1 + 0.4*G_2) * cc.nu[i] * exp(-c*a*cc.nu2[i]/2 - 0.03*cc.nu[i]^0.6/(cc.n_eff[i] + 3)^2)
    # end
    # return result
end

peacock2007_fσ(nu2;a=1.529,b=0.704,c=0.412) = (nu=sqrt.(nu2); d = @. 1 + a*nu^b; @. nu * exp(-c * nu2) * (2 * c * d * nu + b * a * nu^(b-1)) / d^2)
peacock2007_fσ(cc::HMFinit;a=1.529,b=0.704,c=0.412) = (d = @. 1 + a*cc.nu^b; @. cc.nu * exp(-c * cc.nu2) * (2 * c * d * cc.nu + b * a * cc.nu^(b-1)) / d^2)

angulo2012_fof_fσ(sigma2;A=0.201,b=1.7,c=1.172,d=2.08) = @. A * ((d/sqrt(sigma2))^b+1) * exp(-c/sigma2)
angulo2012_fof_fσ(cc::HMFinit;A=0.201,b=1.7,c=1.172,d=2.08) = @. A * ((d/cc.sigma)^b+1) * exp(-c/cc.sigma2)

angulo2012_bound_fσ(sigma2;A=0.265,b=1.9,c=1.4,d=1.675) = angulo2012_fof_fσ(sigma2;A=A,b=b,c=c,d=d)
angulo2012_bound_fσ(cc::HMFinit;A=0.265,b=1.9,c=1.4,d=1.675) = angulo2012_fof_fσ(cc;A=A,b=b,c=c,d=d)
#### high level functions
# requires only nu2
# nu_models=Dict([
#     ("SMT",Dict([("func",sheth2001_fσ)])),
#     ("PS",Dict([("func",press1972_fσ)]))
# ])

# # requires only sigma2
# sigma_models=Dict([
#     ("Jenkins",Dict([("func",jenkins2001_fσ)])),
#     ("Warren",Dict([("func",warren2006_fσ)])),
# ])

# # requires nu2 and sigma2
# sigma_nu_models=Dict([
#     ("Reed03",Dict([("func",reed2003_fσ)]))
# ])

hmf_models=Dict([
    ("SMT",Dict([("func",sheth2001_fσ)])),
    ("PS",Dict([("func",press1972_fσ)])),
    ("Jenkins",Dict([("func",jenkins2001_fσ)])),
    ("Warren",Dict([("func",warren2006_fσ)])),
    ("Reed03",Dict([("func",reed2003_fσ)])),
    ("Reed07",Dict([("func",reed2007_fσ)])),
    ("Watson_FoF",Dict([("func",watson2013_fof_fσ)])),
    ("Peacock",Dict([("func",peacock2007_fσ)])),
    ("Angulo_FoF",Dict([("func",angulo2012_fof_fσ)])),
    ("Angulo_Bound",Dict([("func",angulo2012_bound_fσ)])),
    ("Crocce",Dict([("func",crocce2010_fσ)])),
    ("Courtin",Dict([("func",courtin2010_fσ)])),
    ("Bhattacharya",Dict([("func",bhattacharya2011_fσ)])),
    ("Pillepich",Dict([("func",pillepich2010_fσ)])),
    ("Ishiyama",Dict([("func",ishiyama2015_fσ)]))
])

function massfunc_dndm(M::Union{Real,AbstractArray},c::AbstractCosmology,model::String="SMT",z::Real=0.0;filt=nothing, kws...)
    !("SMT" in keys(hmf_models)) && error("Invalid model, $model.")
    filt === nothing && (filt=c.filt)
    # R_L = lagrangianR.(c,M)
    rho_mean = u.ustrip(ρ_m(ua.Msun/ua.Mpc^3,c,0.0)) ./ c.h^2 # if we divide by c.h^2 here, we get the same result as from python HMF. 
    R = filter_m_to_r(M,rho_mean,filt)
    sigma2 = σ2(R, c, z; filt=filt, kws...)
    nu2 = δc(c,z)^2 ./ sigma2
    dlns_dlnm = dlnσ_dlnm(R,c,z; kws...)
    n_eff = @. -3 * (2 * dlns_dlnm + 1)
    lnsigma = @. -log(sigma2) / 2
    obj = HMFinit(sqrt.(sigma2),sigma2,lnsigma,dlns_dlnm,n_eff,sqrt.(nu2),nu2,z,Ω_m(c,z))
    fσ = hmf_models[model]["func"](obj)
    
    dndm = @. rho_mean/M^2 * fσ * abs(dlns_dlnm)
    # dn_dlnm = M .* dndm
    # dn_dlog10m = log(10) .* dn_dlnm
    # nless = cumtrapz(
    return dndm
end


# cumtrapz(reverse(py_mf.m),reverse(massfunc_dndm(py_mf.m,a,"SMT")))
# cumtrapz(-reverse(py_mf.m),reverse(massfunc_dndm(py_mf.m,a,"SMT")))
# reverse(cumtrapz(-reverse(py_mf.m),reverse(massfunc_dndm(py_mf.m,a,"SMT"))))
# reverse(cumtrapz(-reverse(log.(py_mf.m)),reverse(py_mf.m.*massfunc_dndm(py_mf.m,a,"SMT"))))
end #module
