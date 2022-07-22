module nonlinear

import Dierckx: Spline1D
import Trapz: trapz
import SpecialFunctions: sinint, cosint
import UnitfulAstro as ua
import Unitful as u
import UnitfulEquivalences as ue
import ..cosmo: eisenstein98zerobaryon, power_spectrum, AbstractCosmology, PkFilter, DefaultFilter, filter_function_k, filter_m_to_r, σ2, dlnσ2_dlnr, ρ_m, growth_integral, Ω_ν, Ω_m, smooth_array_gaussian, logspace, setup_growth
import ..peaks: δc, nonlinear_scale, lagrangianR
import ..halo: NFWmu, NFW_fundamental_parameters
import ..hmf: fσ, Sheth2001

# include("utils.jl")

export HMCode2020, Barkana2001_WDM, nonlinear_power_spectrum

abstract type NonlinearModel end
Base.Broadcast.broadcastable(m::NonlinearModel) = Ref(m)
abstract type HMCode2020 <: NonlinearModel end
abstract type HMCode2020_feedback <: NonlinearModel end
abstract type Barkana2001_WDM <: NonlinearModel end


""" 
    dewiggle_mead2020(k::AbstractVector, Pk::AbstractVector, h::Real, Om0::Real, Ob0::Real, Tcmb0::Real, ns::Real)
    dewiggle_mead2020(c::AbstractCosmology; kws...)
Perform the BAO dampening according to the HMCode2020 recipe and return the smooth power spectrum, Pdwl in Equation 15 of Mead 2021. `kmin_wiggle` is the minimum wavenumber in h/Mpc to smooth, `kmax_wiggle` is the maximum wavenumber to smooth in h/Mpc, `wiggle_sigma` is the Gaussian smoothing width in log(k), and `knorm` is the wavenumber  in h/Mpc at which to force the linear and no-wiggle power spectra to match. """
function dewiggle_mead2020(k::AbstractVector, Pk::AbstractVector, h::Real, Om0::Real, Ob0::Real, Tcmb0::Real, ns::Real, filt::Type{<:PkFilter}; kmin_wiggle::Real=5e-3,kmax_wiggle::Real=5,wiggle_sigma::Real=0.25,knorm::Real=0.03)# , growth_factor::Real=1.0)
    @assert length(k) == length(Pk)
    logk=log.(k)

    # construct interpolators and figure out the normalization point. 
    # Pk_interp = Spline1D(log10.(k),log10.(Pk);k=3,bc="zero",s=0.0)
    # Pk_func(k) = exp10.(Pk_interp(log10.(k)))
    # Pk_atnorm = Pk_func(knorm_nowiggle)
    # nowiggle_atnorm = eisenstein98zerobaryon(knorm_nowiggle, h, Om0, Ob0, Tcmb0)^2 * knorm_nowiggle^ns # * growth_factor^2
    # nowiggle_norm = Pk_atnorm / nowiggle_atnorm

    # instead of constructing interpolators, just use the nearest point in the k vector and normalize there
    _,idx = findmin(abs.(k.-knorm))
    nowiggle_atnorm = eisenstein98zerobaryon(k[idx], h, Om0, Ob0, Tcmb0)^2 * knorm^ns # * growth_factor^2
    nowiggle_norm = Pk[idx] / nowiggle_atnorm
    
    # calculate the nowiggle power spectrum. 
    nowiggle = eisenstein98zerobaryon(k, h, Om0, Ob0, Tcmb0).^2 .* k.^ns .* nowiggle_norm # .* growth_factor.^2

    # find the indices between kmin_wiggle and kmax_wiggle
    good_idx=findall((k .> kmin_wiggle) .& (k.<kmax_wiggle))
    Pk_smt=deepcopy(Pk)
    Pk_smt[good_idx] .= smooth_array_gaussian(logk[good_idx],Pk[good_idx] ./ nowiggle[good_idx], wiggle_sigma;nsig=3) .* nowiggle[good_idx]

    # now calculate equation 15
    y_arr = @. Pk * filter_function_k(k,0.0,filt)^2 * k 
    sigma2V0 = trapz(logk, y_arr) / π^2 / 6.0
    Pk_dwl = @. Pk - (1 - exp(-k^2*sigma2V0)) * (Pk - Pk_smt) 
    return Pk_dwl
end
dewiggle_mead2020(c::AbstractCosmology;kws...) = dewiggle_mead2020(c.k,c.Pk,c.h,c.Ω_m,c.Ω_b,c.Tcmb0,c.ns,c.filt;kws...)

"""
    function NFW_window(k,r_vir,conc)
Analytic Fourier transform of the NFW profile used for the HMCode2020 nonlinear model. """
function NFW_window(k,r_vir,conc)
    ks = k*r_vir/conc # Scale wavenumber

    # sine and cosine integrals
    si1 = sinint(ks)
    si2 = sinint((1+conc)*ks)
    ci1 = cosint(ks)
    ci2 = cosint((1+conc)*ks)

    p1 = cos(ks)*(ci2-ci1)
    p2 = sin(ks)*(si2-si1)
    p3 = sin(ks*conc)/(ks*(1+conc))
    
    # calculate window function and divide by NFWmu
    winnfw = p1+p2-p3
    winnfw /= NFWmu(conc)
    winnfw > 1 && (winnfw=one(winnfw))
    winnfw < 0 && (winnfw=zero(winnfw))
    return winnfw
end


"""
    nonlinear_power_spectrum(model::Type{HMCode2020},c::AbstractCosmology,z::Real=0.0; filt::Type{<:PkFilter}=DefaultFilter)
Given the linear power spectrum saved in `c`, calculate the power spectrum at redshift `z` applying the nonlinear correction of the HMCode2020 model. 
"""
function nonlinear_power_spectrum(model::Type{HMCode2020},c::AbstractCosmology,z::Real=0.0; filt::Type{<:PkFilter}=DefaultFilter)
    k, Pk, h, Om0, Ob0, Tcmb0, ns, m_nu = c.k,c.Pk,c.h,c.Ω_m,c.Ω_b,c.Tcmb0,c.ns,c.m_nu
    g_z = c.growth_function(z)
    sum(m_nu) != 0 ? f_ν = Ω_ν(c,z) / (Ω_m(c,z) + Ω_ν(c,z)) : f_ν=0.0
    # set up constants and such
    sigma8 = c.sigma8
    logk = log.(k)
    kd = 0.05699*sigma8^-1.089
    f_damp = 0.2696 * sigma8^0.9403
    nd = 2.853
    kstar = 0.05618 * sigma8^-1.013
    η = 0.1281 * sigma8^-0.3644
    B = 5.196
    F = 0.01
    z_c = 10.0
    delta_c = δc(c,z)
    rho_mean0 = u.ustrip(ρ_m(ua.Msun/ua.Mpc^3,c,0.0)) # / h^2 # if we divide by h^2 here, we get the same result as from python HMF.

    # get the dewiggle'd power spectrum, Equation 15
    Pk_dwl = dewiggle_mead2020(c)

    # set up the two-halo term
    r_nl = nonlinear_scale(c,z) # the non-linear scale in Mpc/h
    n_eff = -dlnσ2_dlnr(r_nl,k,Pk) - 3
    α = 1.875 * 1.603^n_eff
    Pk_2H = @. Pk_dwl * (1 - (f_damp * (k/kd)^nd / (1 + (k/kd)^nd)))

    ########################################
    # set up the 1 halo term
    M_arr=logspace(1e5,1e18,length=100) # these are M200c for Bullock concentration
    
    # get the bullock 2001 concentrations with the changed scaling
    # R = lagrangianR.(c,M_arr * F)
    R = filter_m_to_r.(M_arr.*F,rho_mean0,filt)
    sigma2_arr = σ2(R, k, Pk; growth_factor=g_z,filt=filt)
    
    z_f = c.growth_function_inv(@. delta_c / sqrt(sigma2_arr) * g_z) 
    z_f[z_f .< z] .= z
    conc = @. B * (1+z_f) / (1+z)
    # calculate the dark energy rescaling
    gi_norm = growth_integral(c,0.0)
    conc_correct = c.growth_function(z_c) * growth_integral(c,z,gi_norm) / (growth_integral(c,z_c,gi_norm) * c.growth_function(z))
    @. conc *= conc_correct

    # get the NFW rs and rhos
    rs,rhos = NFW_fundamental_parameters(c,M_arr,conc,z,"vir")
    r_vir = conc .* rs ./ 1e3 # convert from kpc to Mpc 

    # set up the halo mass function
    R .= filter_m_to_r.(M_arr,rho_mean0,filt)
    sigma2_arr .= σ2(R, k, Pk; growth_factor=g_z,filt=filt)
    ν2_arr = delta_c^2 ./ sigma2_arr
    ν_arr = sqrt.(ν2_arr)
    fsigma = fσ.(Sheth2001,ν2_arr)

    integrand1 = (M_arr .* NFW_window.(k * ν_arr'.^η, ones(eltype(k),length(k)) * r_vir', ones(eltype(k),length(k)) * conc')'.^2 .* fsigma ./ rho_mean0 ./ ν_arr .* (1 - f_ν)^2)' # using matrix operations to avoid allocations
    Pk_1H = [trapz(ν_arr, view(integrand1,i,:)) for i in eachindex(k)]
    @. Pk_1H = Pk_1H * (k/kstar)^4 / (1+(k/kstar)^4) 
    
    Δ²2H = @. 4π * (k/(2π))^3 * Pk_2H
    Δ²1H = @. 4π * (k/(2π))^3 * Pk_1H
    result = @. ((Δ²1H^α + Δ²2H^α)^(1/α))
    @. result = result * (2π/k)^3 / (4π) * g_z^2
    return result
end


# now with all options
# nonlinear_power_spectrum(HMCode2020, c.k, c.Pk, c.h, c.Ω_m, c.Ω_k, c.Ω_Λ, c.Ω_b, c.Tcmb0, c.ns, c.Neff, c.m_nu, -1, 0, 0.0; filt=TopHat)
"""
    nonlinear_power_spectrum(model::Type{HMCode2020}, k::AbstractArray, Pk::AbstractArray, h::Real, Ωm::Real, Ωk::Real, ΩΛ::Real, Ωb::Real, Tcmb0::Union{Real,u.Quantity}, ns::Real, Neff::Real, m_nu, w0, wa, z::Real=0.0; filt::Type{<:PkFilter}=DefaultFilter)

Takes a present-day (`z`=0) linear power spectrum `Pk` sampled at wavenumbers `k` for the given cosmological parameters and evolves it to redshift `z` applying the HMCode2020 nonlinear corrections. """
# """If providing growth functions (through `growth_function` and `growth_function_inv`), `growth_function` must take a redshift `z::Real` and return the linear growth factor, while `growth_function_inv` must take a linear growth factor `D` and return the redshift `z`; both, or neither, must be provided; providing only one is unsupported. """
function nonlinear_power_spectrum(model::Type{HMCode2020}, k::AbstractArray, Pk::AbstractArray, h::Real, Ωm::Real, Ωk::Real, ΩΛ::Real, Ωb::Real, Tcmb0::Union{Real,u.Quantity}, ns::Real, Neff::Real, m_nu, w0, wa, z::Real=0.0; filt::Type{<:PkFilter}=DefaultFilter)
    Tcmb0 isa u.Quantity && (Tcmb0 = u.ustrip(Tcmb0|>u.K))
    Ωγ = 4.481620089297254e-7 * Tcmb0^4 / h^2 # calculate photon density from Tcmb0
    # set up growth function
    growth_function_tmp, growth_function_inv_tmp = setup_growth(Ωm, ΩΛ, Ωk, Ωγ, Neff, m_nu, Tcmb0, w0, wa) 
    growth_function(z) = growth_function_tmp((1 ./ (1 .+ z)))
    growth_function_inv(g) = 1 ./ growth_function_inv_tmp(g) .- 1
    # set up constants and such
    g_z = growth_function(z)
    sum(m_nu) != 0 ? f_ν = Ω_ν(c,z) / (Ω_m(c,z) + Ω_ν(c,z)) : f_ν=0.0
    sigma8 = sqrt(σ2(8.0,k,Pk;filt=filt,growth_factor=1.0))
    logk = log.(k)
    kd = 0.05699*sigma8^-1.089
    f_damp = 0.2696 * sigma8^0.9403
    nd = 2.853
    kstar = 0.05618 * sigma8^-1.013
    η = 0.1281 * sigma8^-0.3644
    B = 5.196
    F = 0.01
    z_c = 10.0
    delta_c = δc(z,h,Ωm,Ωk,ΩΛ,Tcmb0,m_nu,Neff,w0,wa)
    rho_mean0 = u.ustrip(ρ_m(ua.Msun/ua.Mpc^3,z,h,Ωm)) # / h^2 # if we divide by h^2 here, we get the same result as from python HMF.

    # set up interpolants
    sigma2_r = logspace(1/maximum(k),1/minimum(k),length=100)
    s2_arr = σ2(sigma2_r, k, Pk; filt=filt, growth_factor=1.0)
    sigma2_interp_tmp = Spline1D(log10.(sigma2_r), log10.(s2_arr); k=4, bc="zero",s=0.0)
    sigma2_interp(x) = exp10.(sigma2_interp_tmp(log10.(x)))
    sigma2_interp_inv_tmp = Spline1D(log10.(reverse(s2_arr)), log10.(reverse(sigma2_r)); k=4, bc="error",s=0.0) # inverse of above
    sigma2_interp_inv(x) = exp10.(sigma2_interp_inv_tmp(log10.(x)))

    # get the dewiggle'd power spectrum, Equation 15
    Pk_dwl = dewiggle_mead2020(k, Pk, h, Ωm, Ωb, Tcmb0, ns, filt)

    # set up the two-halo term
    r_nl = sigma2_interp_inv( (delta_c / g_z)^2) # the non-linear scale in Mpc
    n_eff = -dlnσ2_dlnr(r_nl,k,Pk) - 3
    α = 1.875 * 1.603^n_eff
    Pk_2H = @. Pk_dwl * (1 - (f_damp * (k/kd)^nd / (1 + (k/kd)^nd)))

    ########################################
    # set up the 1 halo term
    M_arr=logspace(1e5,1e18,length=100) # these are M200c for Bullock concentration
    
    # get the bullock 2001 concentrations with the changed scaling
    R = filter_m_to_r.(M_arr.*F,rho_mean0,filt)
    sigma2_arr = σ2(R, k, Pk; growth_factor=g_z,filt=filt)

    z_f = growth_function_inv(@. delta_c / sqrt(sigma2_arr) * g_z) 
    z_f[z_f .< z] .= z
    conc = @. B * (1+z_f) / (1+z)
    # calculate the dark energy rescaling
    gi_norm = growth_integral(z, h, Ωm, Ωk, ΩΛ, Tcmb0, w0, wa)
    conc_correct = growth_function(z_c) * growth_integral(z, h, Ωm, Ωk, ΩΛ, Tcmb0, w0, wa, gi_norm) /
        (growth_integral(z_c, h, Ωm, Ωk, ΩΛ, Tcmb0, w0, wa, gi_norm) * growth_function(z))
    @. conc *= conc_correct

    # get the NFW rs and rhos
    rs, rhos = NFW_fundamental_parameters(M_arr, conc, z, "vir", h, Ωm, Ωk, ΩΛ, Tcmb0, m_nu, Neff, w0, wa)
    r_vir = conc .* rs ./ 1e3 # convert from kpc to Mpc

    # set up the halo mass function
    R .= filter_m_to_r.(M_arr,rho_mean0,filt)
    sigma2_arr .= σ2(R, k, Pk; growth_factor=g_z,filt=filt)
    ν2_arr = delta_c^2 ./ sigma2_arr
    ν_arr = sqrt.(ν2_arr)
    fsigma = fσ.(Sheth2001,ν2_arr)

    # integrand1=[Vector{Float64}(undef,length(M_arr)) for i in eachindex(k)]
    # for i in eachindex(k), j in eachindex(M_arr)
    #     integrand1[i][j] = M_arr[j] * NFW_window(k[i] * ν_arr[j]^η, r_vir[j], conc[j])^2 * fsigma[j] / rho_mean0 / ν_arr[j] * (1 - f_ν)^2
    # end
    # Pk_1H = [trapz(ν_arr, integrand1[i]) for i in eachindex(integrand1)]

    # (M_arr .* NFW_window.(k[i] .* ν_arr.^η, r_vir, conc).^2 .* fsigma ./ rho_mean0 ./ ν_arr .* (1 - f_ν)^2 for i in eachindex(k))
    # Pk_1H = collect((trapz(ν_arr, @. M_arr * NFW_window(k[i] * ν_arr^η, r_vir, conc)^2 * fsigma / rho_mean0 / ν_arr * (1 - f_ν)^2) for i in eachindex(k)))
    
    integrand1 = (M_arr .* NFW_window.(k * ν_arr'.^η, ones(eltype(k),length(k)) * r_vir', ones(eltype(k),length(k)) * conc')'.^2 .* fsigma ./ rho_mean0 ./ ν_arr .* (1 - f_ν)^2)' # using matrix operations to avoid allocations
    Pk_1H = [trapz(ν_arr, view(integrand1,i,:)) for i in eachindex(k)]
    # Pk_1H = [ trapz(ν_arr, @. M_arr * NFW_window(k[i] * ν_arr^η, r_vir, conc)^2 * fsigma / rho_mean0 / ν_arr * (1 - f_ν)^2) for i in eachindex(k) ]
    # Pk_1H = [ trapz(ν_arr, M_arr .* NFW_window.(k[i] .* ν_arr.^η, r_vir, conc).^2 .* fsigma ./ rho_mean0 ./ ν_arr .* (1 - f_ν)^2) for i in eachindex(k) ]
    @. Pk_1H = Pk_1H * (k/kstar)^4 / (1+(k/kstar)^4) 
    
    Δ²2H = @. 4π * (k/(2π))^3 * Pk_2H
    Δ²1H = @. 4π * (k/(2π))^3 * Pk_1H
    result = @. ((Δ²1H^α + Δ²2H^α)^(1/α))
    @. result = result * (2π/k)^3 / (4π) * g_z^2
    return result
end


################################################################################################
######## warm dark matter models ###############################################################

function nonlinear_power_spectrum(model::Type{Barkana2001_WDM},c::AbstractCosmology, m_x::Real, g_x::Real, z::Real=0.0; ϵ=0.361, η=5, ν=1.2)
    Ωdm = c.Ω_m - c.Ω_b
    R0c = 0.201 * ( Ωdm * c.h^2 / 0.15)^0.15 * (g_x / 1.5)^-0.29 * m_x^-1.15 # in Mpc, Equation 4
    g_z = c.growth_function(z)
    return @. c.Pk * (1 + (ϵ*c.k*R0c)^2ν)^(-η/ν) * g_z^2
end







end # module
