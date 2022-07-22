abstract type TransferModel end
Base.Broadcast.broadcastable(m::TransferModel) = Ref(m)
struct Sugiyama1995 <: TransferModel end
struct Eisenstein1998 <: TransferModel end
struct Eisenstein1998_ZB <: TransferModel end
struct Eisenstein1999 <: TransferModel end
struct CAMB <: TransferModel end
const DefaultTransferModel = Eisenstein1999

abstract type PkFilter end
Base.Broadcast.broadcastable(m::PkFilter) = Ref(m)
struct TopHat <: PkFilter end
struct Gaussian <: PkFilter end
struct SharpK <: PkFilter end
# struct SharpK_Benson2013 <: PkFilter end # Andrew Benson's sharp-k filter for WDM
const DefaultFilter = TopHat

"""
    transfer_function(k, h, Om0, Ob0, Onu0, Neff, Tcmb0, model::Type{TransferModel})
    transfer_function(k,c::AbstractCosmology,model::Type{<:TransferModel}=DefaultTransferModel)
The transfer function transforms the spectrum of primordial fluctuations into the linear power spectrum of the matter density fluctuations. The primordial power spectrum is usually described as a power law, leading to a power spectrum P(k) = T(k)^2 k^{ns} where P(k) is the matter power spectrum, T(k) is the transfer function, and :math:`n_s` is  the tilt of the primordial power spectrum. 

Parameters
-------------------------------------------------------------------------------------------
k: array_like
	The wavenumber k (in comoving h/Mpc).
h: Real
	The Hubble constant in units of 100 km/s/Mpc.
Om0: Real
	The matter density in units of the critical density at z = 0.
Ob0: Real
	The baryon density in units of the critical density at z = 0.
Onu0: Real
	The neutrino density in units of the critical density at z = 0.
Neff: Real
	The effective number of neutrino species.
Tcmb0: Real (in Kelvin) or Unitful Quantity
	The temperature of the CMB at z = 0.

Returns
-------------------------------------------------------------------------------------------
Tk: array_like
	The transfer function; has the same dimensions as ``k`` in units of (Mpc/h)^3. """
transfer_function(k, h, Om0, Ob0, Onu0, Neff, Tcmb0, model::Type{Sugiyama1995}) = sugiyama95(k,h,Om0,Ob0,Tcmb0)
transfer_function(k, h, Om0, Ob0, Onu0, Neff, Tcmb0, model::Type{Eisenstein1998}) = eisenstein98(k,h,Om0,Ob0,Tcmb0)
transfer_function(k, h, Om0, Ob0, Onu0, Neff, Tcmb0, model::Type{Eisenstein1999}) = eisenstein99(k,h,Om0,Ob0,Onu0,Neff,Tcmb0)
transfer_function(k, h, Om0, Ob0, Onu0, Neff, Tcmb0, model::Type{Eisenstein1998_ZB}) = eisenstein98zerobaryon(k, h, Om0, Ob0, Tcmb0)
transfer_function(k,c::AbstractCosmology,model::Type{<:TransferModel}=DefaultTransferModel) = transfer_function(k,c.h,c.Ω_m,c.Ω_b,c.Ω_n,c.Neff,c.Tcmb0,model)
################################################################################################
# models

function sugiyama95(k,h,Om0,Ob0,Tcmb0)
    k *= h    # Convert kh from h/Mpc to 1/Mpc
    Tcmb0 isa u.Quantity && (Tcmb0 = Tcmb0 |> u.K |> u.ustrip)

    # this explicit loop is 13% faster but harder to read than the broadcasted code below.
    # not a function of k, so eval outside loop
    q_factor = (Tcmb0 / 2.7)^2 / (Om0 * h^2 * exp(-Ob0 * (1.0 + sqrt(2 * h) / Om0))) 
    if isscalar(k)
        q = k * q_factor# (Tcmb0 / 2.7)^2 / (Om0 * h^2 * exp(-Ob0 * (1.0 + sqrt(2 * h) / Om0)))
        q < 1e-9 && return 1.0
        return log(1.0 + 2.34 * q) / (2.34 * q) *
            (1.0 + 3.89 * q + (16.1 * q)^2 + (5.46 * q)^3 + (6.71 * q)^4)^-0.25
    else
        # Tk=one.(k) #don't allocate Tk, just mutate the input k array in place and return it.
        # doesn't mutate outer k; I think because we are redefining k=k*h above? 
        for i in eachindex(k)
            q = k[i] * q_factor # (Tcmb0 / 2.7)^2 / (Om0 * h^2 * exp(-Ob0 * (1.0 + sqrt(2 * h) / Om0)))
            q<1e-9 ? k[i]=1.0 : k[i] = log(1.0 + 2.34 * q) / (2.34 * q) *
                (1.0 + 3.89 * q + (16.1 * q)^2 + (5.46 * q)^3 + (6.71 * q)^4)^-0.25
        end
        return k
    end
    # q = @. k * (Tcmb0 / 2.7)^2 / (Om0 * h^2 * exp(-Ob0 * (1.0 + sqrt(2 * h) / Om0)))
	
    # Tk = @. log(1.0 + 2.34 * q) / (2.34 * q) *
    #     (1.0 + 3.89 * q + (16.1 * q)^2 + (5.46 * q)^3 + (6.71 * q)^4)^-0.25

    # # Numerically, very small values of q lead to issues with T become zero rather than one.
    # if length(k) == 1
    #     q < 1e-9 && (Tk=1.0)
    # else
    #     Tk[q.<1e-9].=1.
    # end
    return Tk
end

function eisenstein98(k, h, Om0, Ob0, Tcmb0)
    Ob0 < 1e-20 && (return eisenstein98zerobaryon(k,h,Om0,Ob0,Tcmb0))
    k *= h    # Convert kh from h/Mpc to 1/Mpc
    Tcmb0 isa u.Quantity && (Tcmb0 = Tcmb0 |> u.K |> u.ustrip)
    # Define shorter expressions
    omc = Om0 - Ob0
    ombom0 = Ob0 / Om0
    h2 = h^2
    om0h2 = Om0 * h2
    ombh2 = Ob0 * h2
    theta2p7 = Tcmb0 / 2.7
    theta2p72 = theta2p7^2
    theta2p74 = theta2p72^2

    # Equation 2
    zeq = 2.50e4 * om0h2 / theta2p74

    # Equation 3
    keq = 7.46e-2 * om0h2 / theta2p72

    # Equation 4
    b1d = 0.313 * om0h2^-0.419 * (1.0 + 0.607 * om0h2^0.674)
    b2d = 0.238 * om0h2^0.223
    zd = 1291.0 * om0h2^0.251 / (1.0 + 0.659 * om0h2^0.828) * (1.0 + b1d * ombh2^b2d)

    # Equation 5
    Rd = 31.5 * ombh2 / theta2p74 / (zd / 1e3)
    Req = 31.5 * ombh2 / theta2p74 / (zeq / 1e3)

    # Equation 6
    s = 2.0 / 3.0 / keq * sqrt(6.0 / Req) * log((sqrt(1.0 + Rd) +
	sqrt(Rd + Req)) / (1.0 + sqrt(Req)))

    # Equation 7
    ksilk = 1.6 * ombh2^0.52 * om0h2^0.73 * (1.0 + (10.4 * om0h2)^-0.95)

    # Equation 10
    q = @. k / 13.41 / keq

    # Equation 11
    a1 = (46.9 * om0h2)^0.670 * (1.0 + (32.1 * om0h2)^-0.532)
    a2 = (12.0 * om0h2)^0.424 * (1.0 + (45.0 * om0h2)^-0.582)
    ac = a1^(-ombom0) * a2^(-ombom0^3)

    # Equation 12
    b1 = 0.944 / (1.0 + (458.0 * om0h2)^-0.708)
    b2 = (0.395 * om0h2)^-0.0266
    bc = 1.0 / (1.0 + b1 * ((omc / Om0)^b2 - 1.0))

    # Equation 15
    y = (1.0 + zeq) / (1.0 + zd)
    Gy = y * (-6.0 * sqrt(1.0 + y) + (2.0 + 3.0 * y) *
	log((sqrt(1.0 + y) + 1.0) / (sqrt(1.0 + y) - 1.0)))

    # Equation 14
    ab = 2.07 * keq * s * (1.0 + Rd)^(-3.0 / 4.0) * Gy

    # Get CDM part of transfer function

    # Equation 18
    f = @. 1.0 / (1.0 + (k * s / 5.4)^4)

    # Equation 20
    C = @. 14.2 / ac + 386.0 / (1.0 + 69.9 * q^1.08)

    # Equation 19
    T0t = @. log(ℯ + 1.8 * bc * q) / (log(ℯ + 1.8 * bc * q) + C * q^2)

    # Equation 17
    C1bc = @. 14.2 + 386.0 / (1.0 + 69.9 * q^1.08)
    T0t1bc = @. log(ℯ + 1.8 * bc * q) / (log(ℯ + 1.8 * bc * q) + C1bc * q^2)
    Tc = @. f * T0t1bc + (1.0 - f) * T0t

    # Get baryon part of transfer function

    # Equation 24
    bb = 0.5 + ombom0 + (3.0 - 2.0 * ombom0) * sqrt((17.2 * om0h2)^2 + 1.0)

    # Equation 23
    bnode = 8.41 * om0h2^0.435

    # Equation 22
    st = @. s / (1.0 + (bnode / k / s)^3)^(1.0 / 3.0)

    # Equation 21
    C11 = @. 14.2 + 386.0 / (1.0 + 69.9 * q^1.08)
    T0t11 = @. log(ℯ + 1.8 * q) / (log(ℯ + 1.8 * q) + C11 * q * q)
    Tb = @. (T0t11 / (1.0 + (k * s / 5.2)^2) + ab / (1.0 + (bb / k / s)^3) * exp(-(k / ksilk)^1.4)) *
	sin(k * st) / (k * st)

    # Total transfer function
    Tk = @. ombom0 * Tb + omc / Om0 * Tc
    return Tk
end
function eisenstein98zerobaryon(k, h, Om0, Ob0, Tcmb0)
    Tcmb0 isa u.Quantity && (Tcmb0 = Tcmb0 |> u.K |> u.ustrip)
    k *= h    # Convert kh from h/Mpc to 1/Mpc
    ombom0 = Ob0 / Om0
    h2 = h^2
    om0h2 = Om0 * h2
    ombh2 = Ob0 * h2
    theta2p7 = Tcmb0 / 2.7
    # Equation 26
    s = 44.5 * log(9.83 / om0h2) / sqrt(1.0 + 10.0 * ombh2^0.75)
    # Equation 31
    alphaGamma = 1.0 - 0.328 * log(431.0 * om0h2) * ombom0 + 0.38 * log(22.3 * om0h2) * ombom0^2
    # Equation 30
    Gamma = @. Om0 * h * (alphaGamma + (1.0 - alphaGamma) / (1.0 + (0.43 * k * s)^4))
    # Equation 28
    q = @. k * theta2p7 * theta2p7 / Gamma
    # Equation 29
    C0 = @. 14.2 + 731.0 / (1.0 + 62.5 * q)
    L0 = @. log(2.0 * ℯ + 1.8 * q)
    Tk = @. L0 / (L0 + C0 * q * q)
    return Tk
end

""" Only valid if cosmology has massive neutrinos."""
function eisenstein99(k, h, Om0, Ob0, Onu0, Neff, Tcmb0)
    k *= h    # Convert kh from h/Mpc to 1/Mpc
    Tcmb0 isa u.Quantity && (Tcmb0 = Tcmb0 |> u.K |> u.ustrip)
    Odm0 = Om0 - Ob0
    Obh2 = Ob0 * h^2
    Omh2 = (Om0+Onu0) * h^2
    f_b = Ob0 / (Om0+Onu0)
    f_dm = Odm0 / (Om0+Onu0)
    f_nu = Onu0 / (Om0+Onu0)
    f_cb = f_dm + f_b
    theta_cmb = Tcmb0 / 2.7

    z_eq = 2.5e4 * Omh2 * theta_cmb^-4 #really 1+z
    k_eq = 7.46e-2 * Omh2 * theta_cmb^-2  #units Mpc^-1 (no h!)

    z_drag_b1 = 0.313 * Omh2^-0.419 * (1.0 + 0.607 * Omh2^0.674) #equation 2
    z_drag_b2 = 0.238 * Omh2^0.223
    z_drag = (1291.0 * Omh2^0.251 / (1.0 + 0.659 * Omh2^0.828) *
        (1.0 + z_drag_b1 * Obh2^z_drag_b2))

    y_d = (1.0 + z_eq) / (1.0 + z_drag) #equation 3

    sound_horizon = 44.5 * log(9.83/Omh2) / sqrt( 1. + 10. * Obh2^0.75 ) * h # in Mpc

    p_dm = 0.25 * (5.0 - sqrt( 1.0 + 24.0 * f_dm) ) #Equation 11
    p_b = 0.25 * (5.0 - sqrt( 1.0 + 24.0 * f_b) ) #Equation 11
    p_nu = 0.25 * (5.0 - sqrt( 1.0 + 24.0 * f_nu) ) #Equation 11
    p_cb = 0.25 * (5.0 - sqrt( 1.0 + 24.0 * f_cb) ) #Equation 11
    alpha_nu = f_dm/f_cb * ( ( 5.0 - 2.0 * (p_dm + p_cb) ) / (5.0 - 4.0 * p_cb) ) *
        ( 1.0 - 0.553 * (f_nu+f_b) + 0.126 * (f_nu + f_b)^3. ) /
        (1.0 - 0.193 * sqrt( f_nu * Neff) + 0.169 * f_nu * Neff^0.2 ) *
        (1.0 + y_d)^(p_cb-p_dm) *
        ( 1.0 + ( (p_dm - p_cb) / 2.0 ) * ( 1.0 + (1.0 / ( (3.0 - 4.0 * p_dm) *
        (7.0 - 4.0 * p_cb ) ) ) ) / (1.0 + y_d) ) #equation 15
    beta_c = 1.0 / (1.0 - 0.949 * (f_nu + f_b) ) #equation 21 

    q = @. k * theta_cmb^2 / Omh2 

    gamma_eff = @. Omh2 * ( sqrt(alpha_nu) + (1. - sqrt(alpha_nu) ) /
        (1.0 + (0.43 * k * sound_horizon / h )^4 ) ) #equation 16
    q_eff = @. k * theta_cmb^2 / gamma_eff #equation 17
    C = @. 14.4 + (325.0 / (1.0 + 60.5 * q_eff^1.11 ) ) #equation 20 
    L = @. log( ℯ + 1.84 * beta_c * sqrt( alpha_nu ) * q_eff ) #equation 19
    q_nu = @. 3.92 * q * sqrt( Neff ) / f_nu
    T_sup = @. L / ( L + C * q_eff^2. )
    B_k = @. 1.0 + ( (1.2 * f_nu^0.64 * Neff^(0.3 + 0.6 * f_nu) ) /
        (q_nu^-1.6 + q_nu^0.8 ) )
    return @. T_sup * B_k

end

"""
    power_spectrum(k::Union{Real,AbstractArray},T::Union{Real,AbstractArray}, ns::Real, growth_factor::Real=1.0)
    power_spectrum(k::Union{Real,AbstractArray},c::AbstractCosmology,model::Type{<:TransferModel}=DefaultTransferModel, z::Real=0.0)
Constructs the power spectrum from wavenumber `k` (in h/Mpc) and the transfer function `T` scaled by the growth factor. Returns the power spectrum in (Mpc/h)^3. This result is unnormalized and must be properly normalized by σ8, unless calculated from a cosmology object, in which case the saved Pk_normalize constant is used. """
power_spectrum(k::Real,T::Real, ns::Real, growth_factor::Real=1.0) = T^2 * k^ns * growth_factor^2
power_spectrum(k::AbstractArray,T::AbstractArray, ns::Real, growth_factor::Real=1.0) = [T[i]^2 * k[i]^ns * growth_factor^2 for i in eachindex(k,T)]
power_spectrum(k::Real,c::AbstractCosmology,model::Type{<:TransferModel}=DefaultTransferModel, z::Real=0.0) = transfer_function(k,c,model)^2 * k^c.ns * c.growth_function(z)^2 * c.Pk_normalize
power_spectrum(k::AbstractArray,c::AbstractCosmology,model::Type{<:TransferModel}=DefaultTransferModel, z::Real=0.0) = (T=transfer_function(k,c,model); ns=c.ns; norm=c.Pk_normalize; gg=c.growth_function; [norm * T[i]^2 * k[i]^ns * gg(z)^2 for i in eachindex(k,T)])


############################################################################################################################# routines related to mass variance ##############################################################################

"""
    filter_function_k(k,R,filt::Type{<:PkFilter}=DefaultFilter)
The Fourier transform of certain filter functions, used mainly to compute matter variance.
k is wavenumber in h/Mpc, R is filter scale in comoving Mpc/h, filt is a type specifying the type of filter; either TopHat (in real space), SharpK (a tophat in Fourier space) or Gaussian (a Gaussian in real and Fourier space). """
filter_function_k(k,R,filt::Type{TopHat}) = (x = k * R; x < 1e-6 ? (return 1.0) : ( x > 1e10 ? (return 0.0) : (return 3.0 / x^3 * (sin(x) - x * cos(x))) ) )
filter_function_k(k,R,filt::Type{Gaussian}) = (x = k * R; return exp(-x^2 / 2))
filter_function_k(k,R,filt::Type{SharpK}) = (x = k * R; return heaviside(1.0-x, 1.0) )
# filter_function_k(k,R,filt::Type{SharpK_Benson2013}) = (x = k * R / 2.5; return heaviside(1.0-x, 1.0) )
filter_function_k(k,R) = filter_function_k(k,R,DefaultFilter)

"""
    filter_function_r(r,R,filt::Type{<:PkFilter}=DefaultFilter)
The filter functions in real space. `r` is the radius of evaluation in comoving Mpc/h, `R` is the filter scale in comoving Mpc/h, filt is a type specifying the type of filter to use; either TopHat (in real space), SharpK (a tophat in Fourier space) or Gaussian (a Gaussian in real and Fourier space).  """
filter_function_r(r,R,filt::Type{TopHat}) = heaviside(R-r, 0.5)
filter_function_r(r,R,filt::Type{Gaussian}) = exp(-r^2 / 2 / R^2) / (2π)^1.5 / R^3
filter_function_r(r,R,filt::Type{SharpK}) = (x = r/R; x<1e-6 && return inv(6π^2*R^3); (sin(x) - x * cos(x)) / (2 * π^2 * r^3) )
filter_function_r(r,R) = filter_function_r(r,R,DefaultFilter)

"""
    filter_m_to_r(M,ρ_mean,filt::Type{<:PkFilter}=DefaultFilter)
Convert a mass (in Msun) to a radius (in Mpc) for a given filter and mean density of the Universe (as calculated by ρ_m). For TopHat filter this is identical to lagrangianR from `cosmo.peaks`. If using a ShapK filter, there is an additional parameter `c` which is typically 2.5. """
filter_m_to_r(M,ρ_mean,filt::Type{TopHat};kws...) = cbrt(3 * M / (4 * π * ρ_mean))
filter_m_to_r(M,ρ_mean,filt::Type{Gaussian};kws...) = cbrt(M / ρ_mean) / sqrt(2π)
filter_m_to_r(M,ρ_mean,filt::Type{SharpK};c=2.5) = 1/c * cbrt(3 * M / (4 * π * ρ_mean))
filter_m_to_r(M,ρ_mean;c=2.5) = filter_m_to_r(M,ρ_mean,DefaultFilter;c=c)

"""
    filter_r_to_M(R,ρ_mean,filt::Type{<:PkFilter}=DefaultFilter)
Convert radius (in Mpc) to a mass (in Msun) for a given filter and mean density of the Universe (as calculated by ρ_m). For TopHat filter this is identical to lagrangianM from `cosmo.peaks`. If using a SharpK filter, there is an additional parameter `c` which is typically 2.5. """
filter_r_to_m(R,ρ_mean,filt::Type{TopHat};kws...) = 4 * π * R^3 * ρ_mean / 3
filter_r_to_m(R,ρ_mean,filt::Type{Gaussian};kws...) = (2π)^1.5 * R^3 * ρ_mean
filter_r_to_m(R,ρ_mean,filt::Type{SharpK};c=2.5) = 4 * π * (c * R)^3 * ρ_mean / 3
filter_r_to_m(R,ρ_mean;c=2.5) = filter_r_to_m(R,ρ_mean,DefaultFilter;c=c)

"""
    dw_dlnkr(k,R,filt::Type{<:PkFilter}=DefaultFilter)
The logarithmic derivative of the window function, dW / dln(k*R). 
k is the wavenumber in h/Mpc and R is the filter scale in comoving Mpc/h. """
dw_dlnkr(k,R,filt::Type{TopHat}) = (x=k*R; x<1e-6 || x>1e10 ? (return 0.0) : (return (9 * x * cos(x) + 3 * (x^2 - 3) * sin(x)) / x^3) )
dw_dlnkr(k,R,filt::Type{Gaussian}) = (x2=(k*R)^2; -x2 * exp(-x2 / 2) )
dw_dlnkr(k,R,filt::Type{SharpK}) = (k==inv(R) ? (return 1.0) : (return 0.0))
dw_dlnkr(k,R) = dw_dlnkr(k,R,DefaultFilter)

"""
    dw_dkr(k,R,filt="tophat")
The linear derivative of the window function, dW / dln(k*R). """
dw_dkr(k,R,filt::Type{<:PkFilter}) = dw_dlnkr(k,R,filt) * k * R
dw_dkr(k,R) = dw_dlnkr(k,R,DefaultFilter) * k * R

"""
    σ2(R::Union{Real,AbstractArray}, Pk_func; growth_factor=1.0, k_min=9.99e-3, k_max=1.99e4, j=0, filt::Type{<:PkFilter}=DefaultFilter, kws...)
    σ2(R::Union{Real,AbstractArray},c::AbstractCosmology,z=0; j=0, filt=nothing, kws...)
    σ2(R::Union{Real,AbstractArray},k::AbstractVector, Pk::AbstractVector; growth_factor=1.0, j=0, filt::Type{<:PkFilter}=DefaultFilter)
The rms variance of the linear density field on a scale R in comoving Mpc/h. R can be scalar or array-like. 
Pk\\_func is a function that takes a wavenumber k in h/Mpc and returns the present-day power spectrum in (Mpc/h)^3. growth\\_factor gives the linear growth factor at the time σ2 is to be evaluated. k\\_min sets the minimum wavenumber for integration, and k_max sets the maximum. j is the order of the moment, 0 for the common mass variance. filt is a type that sets which filter to use to compute σ2; either TopHat (in real space), Gaussian (in real and k-space), or SharpK (a tophat in k-space). If no filter is provided when giving an AbstractCosmology object, it will use the filter used to construct the interpolant in the object and will be fast. 
kws... are passed to quadgk. By default maxevals=100 is already set for intergration. """
function σ2(R::Real, Pk_func; growth_factor=1.0, k_min=9.99e-3, k_max=1.99e4, j=0, filt::Type{<:PkFilter}=DefaultFilter, kws...)
    if filt == TopHat && j > 0
	throw(DomainError(filt,"Higher-order moments of sigma are not well-defined for tophat filter. Choose gaussian filter instead."))
        return
    end

    function log_integrand(k)
        ret = Pk_func(k) * filter_function_k(k,R,filt)^2 * k^3
        j > 0 && (ret*=k^(2*j))
        return ret
    end

    result = quadgk(x->log_integrand(exp(x)), log(k_min), log(k_max); maxevals=100, kws...)[1] / π^2 / 2.0
    return result*growth_factor^2
end
σ2(R::AbstractArray, Pk_func; kws...) = [σ2(i,Pk_func; kws...) for i in R]
function σ2(R::Real,c::AbstractCosmology,z=0; j=0, filt=nothing, kws...)
    if j == 0 && (filt==c.filt || filt===nothing) # use the precomputed interpolator
        return c.σ2_interp(R) * c.growth_function(z)^2
    else
        filt===nothing && (filt=c.filt)
        return σ2(R,c.Pk_interp; growth_factor=c.growth_function(z), k_min=minimum(c.k), k_max=maximum(c.k), j=j, filt=filt, kws...)
    end
end
function σ2(R::AbstractArray,c::AbstractCosmology,z=0; j=0, filt=nothing, kws...)
    if j == 0 && (filt==c.filt || filt===nothing) # use the precomputed interpolator
        return c.σ2_interp(R) .* c.growth_function(z).^2
    else
        filt===nothing && (filt=c.filt)
        return [σ2(i,c.Pk_interp; growth_factor=c.growth_function(z), k_min=minimum(c.k), k_max=maximum(c.k), j=j, filt=filt, kws...) for i in R]
    end
end
function σ2(R::Real,k::AbstractVector, Pk::AbstractVector; growth_factor=1.0, j=0, filt::Type{<:PkFilter}=DefaultFilter)
    if filt == TopHat && j > 0
	throw(DomainError(filt,"Higher-order moments of sigma are not well-defined for tophat filter. Choose gaussian filter instead."))
        return
    end
    logk = log.(k)
    integrand = @. Pk * filter_function_k(k,R,filt)^2 * k^3
    j > 0 && (integrand*=k^(2*j))
    return trapz(log.(k), integrand) / π^2 / 2.0 * growth_factor^2
end
σ2(R::AbstractArray, k::AbstractVector, Pk::AbstractVector; kws...) = [σ2(i,k,Pk; kws...) for i in R]

"""
    dlnσ2_dlnr(R::Union{Real,AbstractArray}, Pk_func; s2=nothing, growth_factor=1.0, k_min=9.99e-3, k_max=1.99e4, j=0, filt::Type{<:PkFilter}=DefaultFilter, kws...)
    dlnσ2_dlnr(R::Union{Real,AbstractArray}, c::AbstractCosmology, z::Union{Real,AbstractArray}=0.0; s2=nothing, j=0, filt="tophat", kws...)
    dlnσ2_dlnr(R::Union{Real,AbstractArray},k::AbstractVector,Pk::AbstractVector; s2=nothing, growth_factor=1.0, j=0, filt::Type{<:PkFilter}=DefaultFilter)
The logarithmic derivative of σ2 with respect to lnR. Options are as in σ2, except for kws s2=nothing; if s2 is not nothing, this is taken to be the value of σ2 at radius R so it doesn't need to be recomputed. """
function dlnσ2_dlnr(R::Real, Pk_func; s2=nothing, growth_factor=1.0, k_min=9.99e-3, k_max=1.99e4, j=0, filt::Type{<:PkFilter}=DefaultFilter, kws...)
    s2===nothing && ( s2 = σ2(R, Pk_func; growth_factor, k_min, k_max, j, filt, kws...) )
    filt === SharpK && (return -Pk_func(1/R) / (2 * π^2 * s2 * R^3)) # simple form for SharpK
    log_integrand(k) = Pk_func(k) * filter_function_k(k,R,filt) * dw_dlnkr(k,R,filt) * k^3
    result = quadgk(x->log_integrand(exp(x)), log(k_min), log(k_max); maxevals=100, kws...)[1] / π^2 / s2
    return result*growth_factor^2
end
dlnσ2_dlnr(R::AbstractArray, Pk_func; s2=nothing, kws...) = (s2===nothing ? (return [dlnσ2_dlnr(i, Pk_func; s2=nothing, kws...) for i in R]) : (@assert length(R) == length(s2); return [dlnσ2_dlnr(R[i], Pk_func; s2=s2[i], kws...) for i in eachindex(R,s2)] ) )
function dlnσ2_dlnr(R::Real, c::AbstractCosmology, z::Real=0.0; s2=nothing, j=0, filt=nothing, kws...)
    if j==0 && (filt==c.filt || filt===nothing)
        # return derivative(c.σ2_loginterp,log10(R)) # .* c.growth_function(z).^2
        return derivative(c.σ2_interp.sigma2_interp_tmp,log10(R)) # .* c.growth_function(z).^2
    else
        filt===nothing && (filt="tophat")
        return dlnσ2_dlnr(R,c.Pk_interp;s2=s2,growth_factor=c.growth_function(z),k_min=minimum(c.k),k_max=maximum(c.k),j=j,filt=filt,kws...)
    end
end
function dlnσ2_dlnr(R::AbstractArray, c::AbstractCosmology, z::Union{Real,AbstractArray}=0.0; s2=nothing, j=0, filt=nothing, kws...)
    if j==0 && (filt==c.filt || filt===nothing)
        # return derivative(c.σ2_loginterp,log10.(R)) # .* c.growth_function(z).^2
        return derivative(c.σ2_interp.sigma2_interp_tmp,log10.(R)) # .* c.growth_function(z).^2
    else
        filt===nothing && (filt=DefaultFilter)
        return dlnσ2_dlnr(R,c.Pk_interp;s2=s2,growth_factor=c.growth_function(z),k_min=minimum(c.k),k_max=maximum(c.k),j=j,filt=filt,kws...)
    end
end
function dlnσ2_dlnr(R::Real,k::AbstractVector,Pk::AbstractVector; s2=nothing, growth_factor=1.0, j=0, filt::Type{<:PkFilter}=DefaultFilter)
    s2===nothing && ( s2 = σ2(R, k, Pk; growth_factor, j, filt) )
    filt === SharpK && (return -Pk[findmin( abs.(k.-(1/R)))[2]] / (2 * π^2 * s2 * R^3)) # simple form for SharpK
    return trapz(log.(k), @. Pk * filter_function_k(k,R,filt) * dw_dlnkr(k,R,filt) * k^3) ./ π^2 ./ s2 .* growth_factor^2
end
dlnσ2_dlnr(R::AbstractArray,k::AbstractVector,Pk::AbstractVector;s2=nothing,kws...) = (s2===nothing ? (return [dlnσ2_dlnr(i,k,Pk; s2=nothing, kws...) for i in R]) : (@assert length(R) == length(s2); return [dlnσ2_dlnr(R[i],k,Pk; s2=s2[i], kws...) for i in eachindex(R,s2)]) )
"""
    dlnσ_dlnr(R::Union{Real,AbstractArray}, Pk_func; kws...)
    dlnσ_dlnr(R::Union{Real,AbstractArray}, c::AbstractCosmology, z::Union{Real,AbstractArray}=0.0; kws...)
    dlnσ_dlnr(R::Union{Real,AbstractArray}, k::AbstractVector, Pk::AbstractVector; kws...)
A convenience function that computes `dlnσ2_dlnr` and divides it by 2 to get `dlnσ_dlnr` """
dlnσ_dlnr(R::Union{Real,AbstractArray}, Pk_func; kws...) = dlnσ2_dlnr(R,Pk_func;kws...) ./ 2
dlnσ_dlnr(R::Union{Real,AbstractArray}, c::AbstractCosmology, z::Union{Real,AbstractArray}=0.0; kws...) = dlnσ2_dlnr(R,c,z;kws...) ./ 2
dlnσ_dlnr(R::Union{Real,AbstractArray}, k::AbstractVector, Pk::AbstractVector; kws...) = dlnσ2_dlnr(R,k,Pk;kws...) ./ 2

"""
    dlnr_dlnm(R)
The derivative of log radius with log mass. For the usual m ∝ r^3 mass assignment, this is just 1/3. """
dlnr_dlnm(R) = 1/3

""" 
    dlnσ2_dlnm(R::Union{Real,AbstractArray}, Pk_func; s2=nothing, growth_factor=1.0, k_min=9.99e-3, k_max=1.99e4, j=0, filt::Type{<:PkFilter}=DefaultFilter, kws...) 
    dlnσ2_dlnm(R::Union{Real,AbstractArray},c::AbstractCosmology,z::Union{Real,AbstractArray}=0.0; kws...)
    dlnσ2_dlnm(R::Union{Real,AbstractArray}, k::AbstractVector, Pk::AbstractVector; kws...)
The logarithmic derivative of σ2 with respect to lnM. Options are as in σ2, except for kws s2=nothing; if s2 is not nothing, this is taken to be the value of σ2 at radius R so it doesn't need to be recomputed."""
dlnσ2_dlnm(R::Real, Pk_func; s2=nothing, growth_factor=1.0, k_min=9.99e-3, k_max=1.99e4, j=0, filt::Type{<:PkFilter}=DefaultFilter, kws...) = dlnσ2_dlnr(R, Pk_func; s2, growth_factor, k_min, k_max, j, filt, kws...) * dlnr_dlnm(R)
dlnσ2_dlnm(R::AbstractArray, Pk_func; kws...) = [dlnσ2_dlnm(i,Pk_func; kws...) for i in R]
dlnσ2_dlnm(R::Real,c::AbstractCosmology,z::Real=0.0; kws...) = dlnσ2_dlnr(R,c,z;kws...) * dlnr_dlnm(R)
dlnσ2_dlnm(R::AbstractArray,c::AbstractCosmology,z::Union{Real,AbstractArray}=0.0; kws...) = dlnσ2_dlnr(R,c,z;kws...) .* dlnr_dlnm.(R)
dlnσ2_dlnm(R::Union{Real,AbstractArray}, k::AbstractVector, Pk::AbstractVector; kws...) = dlnσ2_dlnr(R,k,Pk;kws...) .* dlnr_dlnm.(R)

"""
    dlnσ_dlnm(R::Union{Real,AbstractArray}, Pk_func; kws...)
    dlnσ_dlnm(R::Union{Real,AbstractArray},c::AbstractCosmology,z::Union{Real,AbstractArray}=0.0; kws...)
    dlnσ_dlnm(R::Union{Real,AbstractArray}, k::AbstractVector, Pk::AbstractVector; kws...)
A convenience function that computes `dlnσ2_dlnm` and divides it by 2 to get `dlnσ_dlnm` """
dlnσ_dlnm(R::Union{Real,AbstractArray}, Pk_func; kws...) = dlnσ2_dlnm(R,Pk_func;kws...) ./ 2
dlnσ_dlnm(R::Union{Real,AbstractArray},c::AbstractCosmology,z::Union{Real,AbstractArray}=0.0; kws...) = dlnσ2_dlnm(R,c,z;kws...) ./ 2
dlnσ_dlnm(R::Union{Real,AbstractArray}, k::AbstractVector, Pk::AbstractVector; kws...) = dlnσ2_dlnm(R,k,Pk;kws...) ./ 2
############################################################################################################################# routines related to correlation functions ######################################################################
""" 
    ξmm(R::Union{Real,AbstractArray}, Pk_func; growth_factor=1.0, k_min=9.99e-3, k_max=1.99e4, kws...)
    ξmm(R::Union{Real,AbstractArray}, Pk_func; kws...)
    ξmm(R::Real, c::AbstractCosmology, z::Real=0.0; kws...)
Matter-matter autocorrelation. """
function ξmm(R::Real, Pk_func; growth_factor=1.0, k_min=9.99e-3, k_max=1.99e4, kws...)
    func(k) = k * Pk_func(k) * sin(k * R) / (2 * π^2 * R)
    # ξ, err = quadde(func,k_min,k_max; kws...) # quaddeo( func, R, 0.0, k_min, k_max) # technically, quaddeo should be optimal here but it tries to make calls outside k_min and k_max so it errors if the Pk_func is an interpolation.
    # ξ, err = quaddeo( func, R, 0.0, k_min, k_max; rtol=1e-2, atol=1e-2, kws...)
    ξ, err = quaddeo( func, R, 0.0, 0.0, Inf; rtol=1e-8, kws...)
    # func(k) = k^3 * Pk_func(k) * sin(k * R) / (2 * π^2 * k * R)
    # ξ, err = quadgk(x->func(exp(x)),log(k_min),log(k_max))
    return ξ * growth_factor^2
end
ξmm(R::AbstractArray, Pk_func; kws...) = [ξmm(i,Pk_func;kws...) for i in R]
ξmm(R::Real, c::AbstractCosmology, z::Real=0.0; kws...) = ξmm(R,c.Pk_interp; growth_factor=c.growth_function(z), k_min=minimum(c.k), k_max=maximum(c.k), kws...)

# function dξmm_dr(R::Real, Pk_func; growth_factor=1.0, k_min=9.99e-3, k_max=1.99e4, ret_interp=false, kws...)
#     func(k) = k * Pk_func(k) * (k * cos(k*R) - sin(k*R)/R) / (2 * π^2 * R)
#     # ξ, err = quadde(func,k_min,k_max; kws...)
#     ξ, err = quaddeo( func, R, 0.0, 0.0, Inf; kws...)
#     # ξ, err = quadgk(func,k_min,k_max; kws...)
#     return ξ * growth_factor^2
# end
# dξmm_dr(R::AbstractArray, Pk_func; kws...) = [dξmm_dr(i,Pk_func;kws...) for i in R]
# dξmm_dr(R::Real, c::AbstractCosmology, z::Real=0.0; kws...) = dξmm_dr(R,c.Pk_interp; growth_factor=c.growth_function(z), k_min=minimum(c.k), k_max=maximum(c.k), kws...)

############################################################################################################################# routines related to peak height ##############################################################################


