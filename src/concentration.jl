module conc_module

import Roots: find_zero, Bisection
import Dierckx: Spline1D, Spline2D, derivative
import Unitful as u
import UnitfulAstro as ua
import HDF5: h5open, close, read, write
# include("peaks.jl")
import ..peaks: δc, νpeak2, lagrangianR, nonlinear_mass, mass_from_νpeak
import ..cosmo: AbstractCosmology, H, isscalar, dlnσ_dlnr, TopHat, σ2
import ..halo: NFWmu, convert_mdef_NFW, NFW_fundamental_parameters, density_threshold


export concentration


const INVALID_CONCENTRATION = -1.0


function concentration(c::AbstractCosmology, M::Union{Real,AbstractArray,u.Quantity,Array{<:u.Quantity}}, z::Union{Real,AbstractArray},mdef::String, model::String="diemer19"; range_warning::Bool=true, range_return::Bool=false, return_mask::Bool=false, kws...)
    @assert (length(M)==length(z) || isscalar(z) || isscalar(M))    
    !(model in keys(models)) && error("Unknown model $model.")
    mdefs_model = models[model]["mdefs"]
    func = models[model]["func"]
    # Now check whether the definition the user has requested is the native definition of the model.
    # If yes, we just return the model concentration. If not, the problem is much harder. Without 
    # knowing the concentration, we do not know what mass in the model definition corresponds to 
    # the input mass M. Thus, we need to find both M and c iteratively.
    if mdef in mdefs_model
        args=[c,M,z]
        length(mdefs_model)>1 && push!(args,mdef)
        # conc,mask=func(c,M,z;kws...)
        conc,mask=func(args...;kws...)
        range_return ? (return conc,mask) : (return conc)

    else
	# ---------------------------------------------------------------------------------------------
	# This equation is zero for a mass MDelta (in the mass definition of the c-M model) when the
	# corresponding mass in the user's mass definition is M_desired.
        function eq(MDelta,M_desired,mdef_model,mdef,func,args;kws...)
            cDelta, mask_tmp = func(args[1],MDelta,args[2:end]...;kws...)
            cDelta < 0.0 && return NaN
            rs,rhos=NFW_fundamental_parameters(args[1], MDelta, cDelta, args[2], mdef_model) #(c::AbstractCosmology, M, conc, z, mdef)
            M_new, R_new, c_new = convert_mdef_NFW(rs,rhos,args[1],args[2],mdef)
            return M_new - M_desired
        end
	guess_factors = [2.0, 5.0, 10.0, 100.0, 10000.0]
        # for the iterative process, it is easiest if M and z are of the same length.
        if isscalar(M) && isscalar(z)
            M,z,mask,conc = [M],[z],BitArray([false]),[0.0]
            scalar_vals = true
        elseif length(M) > length(z)
            z,mask,conc=zero(M).+z, BitArray(zeros(length(M))), zero(M)
            scalar_vals = false
        elseif length(z) > length(M)
            M,mask,conc=zero(z).+M, BitArray(zeros(length(M))), zero(z)
            scalar_vals = false
        else
            mask, conc = BitArray(zeros(length(M))), zero(z)
            scalar_vals = false
        end
        # length(M)>1 ? (mask = similar(M,BitArray)) : mask=false
	# conc = zero(M)
        mdef_model = mdefs_model[1]
	# To a good approximation, the relation M2 / M1 = Delta1 / Delta2. We use this mass
	# as a guess around which to look for the solution.
        delta_ratio = density_threshold.(c,z,mdef) ./ density_threshold.(c,z,mdef_model)
        M_guess = M .* delta_ratio
        for i in eachindex(M,z)
            j=1
            MDelta = nothing
            length(mdefs_model)>1 ? (args = [c,z[i],mdef_model]) : (args = [c,z[i]])
            while MDelta === nothing && j < length(guess_factors)
                M_min = M_guess[i] / guess_factors[j]
		M_max = M_guess[i] * guess_factors[j]

                eq_min = eq(M_min,M[i],mdef_model,mdef,func,args;kws...)
                eq_max = eq(M_max,M[i],mdef_model,mdef,func,args;kws...)
		(isnan(eq_min) || isnan(eq_max)) && break						
		if eq_min * eq_max < 0.0
		    # MDelta = scipy.optimize.brentq(eq, M_min, M_max, args = args_solver)
                    # find_zero(f,(zmin,zmax), Bisection(); kws...)
                    MDelta = find_zero(x->eq(x,M[i],mdef_model,mdef,func,args;kws...),(M_min,M_max);rtol=1e-3) 
		else
		    j += 1
                end
            end
	    if MDelta === nothing || MDelta < 0.1
		range_warning && println("Could not find concentration")#for model $model, mass , mdef %s.' % (model, M_array[i], mdef)
		c[i] = INVALID_CONCENTRATION
		mask[i] = false
            else
                args=[c,MDelta,z[i]]
                length(mdefs_model)>1 && push!(args,mdef_model)
                cDelta, mask_element = func(args...;kws...)
                rs,rhos=NFW_fundamental_parameters(c, MDelta, cDelta, z[i], mdef_model) #(c::AbstractCosmology, M, conc, z, mdef)
                M_new, R_new, c_new = convert_mdef_NFW(rs,rhos,c,z[i],mdef)
                conc[i] = c_new
                mask[i] = mask_element
            end
        end
            # range_warning && (0 in mask)
        scalar_vals && (conc=conc[1]; mask=mask[1])
        range_return ? (return conc,mask) : (return conc)
    end
end



"""
    bullock2001_conc(c::AbstractCosmology, M200c::Union{AbstractArray,Real,u.Quantity,Array{<:u.Quantity}}, z::Union{Real,AbstractArray}; kws...)
"""
function bullock2001_conc(c::AbstractCosmology, M200c::Union{Real,AbstractArray}, z::Union{Real,AbstractArray}; kws...)
    @assert (length(M200c)==length(z) || isscalar(z) || isscalar(M200c))    
    K=3.85
    F=0.01
    # The math works out such that we are looking for the redshift where the growth factor is
    # equal to the peak height of a halo with mass F * M. Utilize an inversion of the growth factor
    # as a function of z to get the z as a function of growth_factor. 
    R = lagrangianR.(c,M200c * F)
    growth_fac = sqrt.(νpeak2(R,c,0.0; s2=nothing, j=0, filt=TopHat, kws...))
    # H0 = H.(c,z) .|> u.ustrip# c.h * 100
    # conc = @. K * (u.ustrip(H(c,c.growth_function_inv(growth_fac))) / H0)^(2/3)
    # new implementation following the text in Mead2021; seems to give better results
    # not sure what I was doing with the H factors tbh
    conc = K .* (1 .+ c.growth_function_inv(growth_fac)) ./ (1 + z)
    if isscalar(M200c) && isscalar(z)
        mask=true
    else
        mask=BitArray(ones(maximum([length(M200c),length(z)])))
    end
    return conc, mask
end
# bullock2001_conc(c::AbstractCosmology, M200c::u.Quantity, z::Real; kws...) = bullock2001_conc(c,M200c|>ua.Msun|>u.ustrip,z;kws...)
bullock2001_conc(c::AbstractCosmology, M200c::Union{u.Quantity, Array{<:u.Quantity}}, z::Union{Real,AbstractArray}; kws...) = bullock2001_conc(c,M200c.|>ua.Msun.|>u.ustrip,z;kws...)

"""
    duffy2008_conc(c::AbstractCosmology, M::Union{AbstractArray,Real,u.Quantity,Array{<:u.Quantity}}, z::Union{AbstractArray,Real}, mdef::String)
Concentration model of Duffy et al. 2008. A power-law fit calibrated for WMAP5 cosmology. Valid `mdef` are '200c', 'vir', and '200m'. """
function duffy2008_conc(c::AbstractCosmology, M::Union{AbstractArray,Real}, z::Union{AbstractArray,Real}, mdef::String)
    @assert (length(M)==length(z) || isscalar(z) || isscalar(M))
    if mdef == "200c"
	A = 5.71
	B = -0.084
	C = -0.47
    elseif mdef == "vir"
	A = 7.85
	B = -0.081
	C = -0.71
    elseif mdef == "200m"
	A = 10.14
	B = -0.081
	C = -1.01
    else
	error("Invalid mass definition for Duffy et al. 2008 model, $mdef")
    end
    c = @. A * (M / 2E12)^B * (1.0 + z)^C
    mask = (M .>= 1E11) .& (M .<= 1E15) .& (z .<= 2.0)
    return c, mask
end
duffy2008_conc(c::AbstractCosmology,M::Union{u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real}, mdef::String) = duffy2008_conc(c,M.|>ua.Msun.|>u.ustrip,z,mdef)


"""
    klypin2011_conc(c::AbstractCosmology,Mvir::Union{AbstractArray,Real,u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real})
Concentration model of Klypin et al. 2011. This power-law fit was calibrated for the WMAP7 cosmology of the Bolshoi simulation. Note that this model relies on concentrations that were measured from circular velocities, rather than from a fit to the actual density profiles. Klypin et al. 2011 also give fits at particular redshifts other than zero. However, there is no clear procedure to interpolate between redshifts, particularly since the z = 0 relation has a different functional form than the high-z relations. Thus, we only implement the `z=0` relation here. """
function klypin2011_conc(c::AbstractCosmology,Mvir::Union{AbstractArray,Real},z::Union{AbstractArray,Real})
    cvir = @. 9.6 * (Mvir / 1e12)^-0.075
    mask = @. (Mvir > 3e10) & (Mvir .< 5E14) & (z < 0.01)
    return cvir, mask
end
klypin2011_conc(c::AbstractCosmology,Mvir::Union{u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real}) = klypin2001_conc(c,Mvir.|>ua.Msun.|>u.ustrip,z)

"""
    prada2012_conc(c::AbstractCosmology,M200c::Union{AbstractArray,Real,u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real})
    Concentration model of Prada et al. 2012. This model predicts concentration based on the c-ν relation. The model was calibrated on the Bolshoi and Multidark simulations, but is in principle applicable to any cosmology. The implementation follows equations 12 to 22 in Prada et al. 2012. This function uses the exact values for peak height rather than their approximation. """
function prada2012_conc(c::AbstractCosmology, M200c::Union{AbstractArray,Real},z::Union{AbstractArray,Real}; kws...)
    @assert (length(M200c)==length(z) || isscalar(z) || isscalar(M200c))
    cmin(x) = 3.681 + (5.033 - 3.681) * (1.0 / π * atan(6.948 * (x - 0.424)) + 0.5)
    smin(x) = 1.047 + (1.646 - 1.047) * (1.0 / π * atan(7.386 * (x - 0.526)) + 0.5)
    R = lagrangianR.(c,M200c)
    peak_height = sqrt.(νpeak2(R,c,z; s2=nothing, j=0, filt=TopHat, kws...))
    a = @. 1 / (1 + z)
    x = @. cbrt(c.Ω_Λ / c.Ω_m) * a
    B0 = @. cmin(x) / cmin(1.393)
    B1 = @. smin(x) / smin(1.393)
    temp_sig = @. 1.686 / peak_height
    temp_sigp = @. temp_sig * B1
    temp_C = @. 2.881 * ((temp_sigp / 1.257)^1.022 + 1) * exp(0.06 / temp_sigp^2)
    c200c = @. B0 * temp_C

    # the body of the do block becomes an anonymous function, which is passed as the first argument to broadcast as
    # broadcast(body,c,z,peak_height) and so will handle various cases.
    # c200c = 
    # broadcast(c,z,M200c) do c,z,M200c
    # broadcast(c,z,a,x,peak_height) do c,z,a,x,peak_height
    #     # R = lagrangianR(c,M200c)
    #     # peak_height = sqrt(νpeak2(R,c,z; s2=nothing, j=0, filt="tophat", kws...))
    #     # a = 1.0 / (1.0 + z)
    #     # x = cbrt(c.Ω_Λ / c.Ω_m) * a
    #     B0 = cmin(x) / cmin(1.393)
    #     B1 = smin(x) / smin(1.393)
    #     temp_sig = 1.686 / peak_height
    #     temp_sigp = temp_sig * B1
    #     temp_C = 2.881 * ((temp_sigp / 1.257)^1.022 + 1) * exp(0.06 / temp_sigp^2)
    #     c200c = B0 * temp_C
    # end
    if isscalar(M200c) && isscalar(z)
        mask=true
    else
        mask=BitArray(ones(length(c200c)))
    end
    return c200c, mask
end
prada2012_conc(c::AbstractCosmology,M200c::Union{u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real}; kws...) = prada2012_conc(c,M200c.|>ua.Msun.|>u.ustrip,z;kws...)


"""
    bhattacharya2013_conc(c::AbstractCosmology,M::Union{AbstractArray,Real,u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real},mdef::String; kws...)
Concentration model of Bhattacharya et al. 2013. This power-law fit in c-ν was calibrated for a WMAP7 cosmology. Valid `mdef` are '200c', 'vir', and '200m'. """
function bhattacharya2013_conc(c::AbstractCosmology,M::Union{Real,AbstractArray},z::Union{AbstractArray,Real},mdef::String; kws...)
    @assert (length(M)==length(z) || isscalar(z) || isscalar(M))
    D = c.growth_function(z)
    # Note that peak height in the B13 paper is defined wrt. the mass definition in question, so 
    # we can just use M to evaluate nu.
    R = lagrangianR.(c,M)
    nu = sqrt.(νpeak2(R,c,z; s2=nothing, j=0, filt=TopHat, kws...))

    if mdef == "200c"
	c_fit = @. 5.9 * D^0.54 * nu^-0.35
    elseif mdef == "vir"
	c_fit = @. 7.7 * D^0.90 * nu^-0.29
    elseif mdef == "200m"
	c_fit = @. 9.0 * D^1.15 * nu^-0.29
    else
        error("Invalid mass definition for Bhattacharya et al. 2013 model, $mdef")
    end
				
    if isscalar(M) && isscalar(z)
        z<=0.5 && M<2e12 ? mask=false : (M>2e15 ? mask=false : mask=true)
    else
        M_min = 2e12
        mask = BitArray(ones(maximum([length(M),length(z)])))
        mask[((z .< 0.5) .& ((M .< M_min) .| (M .> 2e15)))] .= false
        mask[((z .> 0.5) .& (z .<= 1.5) .& ((M .< M_min) .| (M .> 2e14)))] .= false
        mask[((z .> 1.5) .& ((M .< M_min) .| (M .> 1e14)))] .= false
    end
    
    return c_fit, mask
end
bhattacharya2013_conc(c::AbstractCosmology,M::Union{u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real},mdef::String; kws...) = bhattacharya2013_conc(c,M.|>ua.Msun.|>u.ustrip,z,mdef;kws...)


"""
    dutton2014_conc(c::AbstractCosmology,M::Union{AbstractArray,Real,u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real},mdef::String)
Concentration model of Dutton et al. 2014. This power-law fit was calibrated for the Planck 2013 cosmology. Valid `mdef` are '200c' and 'vir'. """
function dutton2014_conc(c::AbstractCosmology,M::Union{Real,AbstractArray},z::Union{AbstractArray,Real},mdef::String)
    @assert (length(M)==length(z) || isscalar(z) || isscalar(M))
    if mdef == "200c"
	a = @. 0.520 + (0.905 - 0.520) * exp(-0.617 * z^1.21)
	b = @. -0.101 + 0.026 * z
    elseif mdef == "vir"
	a = @. 0.537 + (1.025 - 0.537) * exp(-0.718 * z^1.08)
	b = @. -0.097 + 0.024 * z
    else
	error("Invalid mass definition for Dutton & Maccio 2014 model, $mdef")
    end	
    logc = @. a + b * log10(M / 1E12)
    c = exp10.(logc)
    mask = (M .>= 1E10) .& (z .<= 5.0)
    return c, mask
end
dutton2014_conc(c::AbstractCosmology,M::Union{u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real},mdef::String) = dutton2014_conc(c,M.|>ua.Msun.|>u.ustrip,z,mdef)

"""
    child2018_conc(c::AbstractCosmology,M200c::Union{AbstractArray,Real,u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real};halo_sample="individual_all")
Concentration model of Child et al. 2018. Valid definitions for `halo_sample` are 'individual_all', 'individual_relaxed', 'stacked_nfw', and 'stacked_einasto'. 
"""
function child2018_conc(c::AbstractCosmology,M200c::Union{Real,AbstractArray},z::Union{AbstractArray,Real};halo_sample::String="individual_all")
    if halo_sample == "individual_all"
	m = -0.10
	A = 3.44
	b = 430.49
	c0 = 3.19
    elseif halo_sample == "individual_relaxed"
	m = -0.09
	A = 2.88
	b = 1644.53
	c0 = 3.54
    elseif halo_sample == "stacked_nfw"
	m = -0.07
	A = 4.61
	b = 638.65
	c0 = 3.59
    elseif halo_sample == "stacked_einasto"
	m = -0.01
	A = 63.2
	b = 431.48
	c0 = 3.36
    else
	error("Unknown halo sample for child18 concentration model, $halo_sample")
    end
    mask = (M200c .>= 2.1e11) .& (z .>= 0.0) .& (z .<= 4.0)
    
    Mstar = nonlinear_mass(c,z)
    M_MT = @. M200c / (Mstar * b)
    c200c = @. c0 + A * (M_MT^m * (1.0 + M_MT)^-m - 1.0)

    return c200c, mask
end
child2018_conc(c::AbstractCosmology,M200c::Union{u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real};halo_sample::String="individual_all") = child2018_conc(c,M200c.|>ua.Msun.|>u.ustrip,z,halo_sample)

#############################################################################################
###### support for Diemer et al. 2019

# this will generate the interpolation table for Diemer 2019. 
function generate_diemer19_table(n_n=40,n_c=80)
    n = range(-4,0.0,length=n_n)
    c = logspace(1e-1,1e3,length=n_c)
    mu = NFWmu.(c)
    # lhs = [@. log10( c / mu^((5 + n) / 6)) for c in c, mu in mu]
    lhs = @. log10( c / mu^((5 + n') / 6)) #column major order
    # At very low concentration and shallow slopes, the LHS begins to rise again. This will cause
    # issues with the inversion. We set those parts of the curve to the minimum concentration of 
    # a given n bin.
    mask_ascending = BitArray(ones(Bool,size(lhs))) # construct mask
    mask_ascending[1:end-1, :] = (mapslices(diff,lhs;dims=1) .> 0.0) # this is equivalent to (np.diff(lhs, axis = 0) > 0.0)
    # mask_ascending[:,1:end-1] = (mapslices(diff,lhs;dims=2) .> 0.0) # this is equivalent to (np.diff(lhs, axis = 0) > 0.0)
    # the mapslices takes slices along the first dimension from lhs and passes it to diff
    
    # Create a table of c as a function of G and n. First, use the absolute min and max of G as 
    # the table range
    G_min = minimum(lhs)
    G_max = maximum(lhs)
    G = range(G_min, G_max, length=n_c)
		
    gc_table = ones(Float64,(n_c, n_n)) .* -10.0
    mins = zero(n)
    maxs = zero(n)
    for i in 1:n_n
        # We interpolate only the ascending values to get c(G)
        mask_ = mask_ascending[:, i]
        lhs_ = lhs[mask_, i]
        mins[i] = minimum(lhs_)
        maxs[i] = maximum(lhs_)
        interp = Spline1D(lhs_,log10.(c[mask_]); k=3, bc="error",s=0.0)
		
        # Not all G exist for all n
        mask = @. (G >= mins[i]) & (G <= maxs[i])
        res = interp(G[mask])
        gc_table[mask, i] .= res
        
        mask_low = (G .< mins[i])
        gc_table[mask_low, i] .= minimum(res)
        mask_high = (G .> maxs[i])
        gc_table[mask_high, i] .= maximum(res)

    end
    return gc_table, G, n, mins, maxs
    # file=h5open("data/diemer19_interp_table.hdf5","w")
    # write(file,"gc_table",gc_table)
    # write(file,"G",collect(G))
    # write(file,"n",collect(n))
    # write(file,"mins",mins)
    # write(file,"maxs",maxs)
    # close(file)
end

# loads the saved interpolation table and loads it into the module's namespace
function load_diemer19_table()
    filename = joinpath(pkgdir(conc_module), "data/diemer19_interp_table.hdf5")
    file=h5open(filename,"r")
    n = read(file["n"])
    G = read(file["G"])
    mins = read(file["mins"])
    maxs = read(file["maxs"])
    gc_table = read(file["gc_table"])
    interp1 = Spline1D(n,mins;bc="error",k=3,s=0.0)
    interp2 = Spline1D(n,maxs;bc="error",k=3,s=0.0)
    interp3 = Spline2D(G,n,gc_table;kx=3,ky=3,s=0.0)
    close(file)
    return interp1, interp2, interp3
end
const diemer19_interp_Gmin, diemer19_interp_Gmax, diemer19_interp_Gc = load_diemer19_table()


# general function for concentration modules of the type used in diemer2019
function _diemer2019_general(c::AbstractCosmology,M200c::Union{Real,AbstractArray},z::Union{AbstractArray,Real},Κ,a0,a1,b0,b1,c_α; kws...)
    R_L = lagrangianR.(c,M200c)
    nu = sqrt.(νpeak2(R_L,c,z; kws...))

    n_eff = -2 .* dlnσ_dlnr(R_L.*Κ, c, z; kws...) .- 3
    alpha_eff = -derivative(c.growth_function.growth_function_tmp, 1. / (1 .+ z)) .* (1 .+ z) ./ c.growth_function(z)
    # Compute input parameters and the right-hand side of the c-M equation. We use interpolation
    # tables to find the concentration at which the equation gives that RHS.    A_n = @. a0 * (1+a1) * (n_eff+3)
    A_n = @. a0 * (1 + a1 * (n_eff + 3))
    B_n = @. b0 * (1 + b1 * (n_eff+3))
    C_α = @. 1 - c_α * (1 - alpha_eff)
    rhs = @. log10(A_n / nu * (1 + nu^2 / B_n) )
    # Mask out values of rhs for which do G not exist. The interpolator will still work because it
    # is based on a full grid of G values, but we mask out the invalid array elements.
    mask = @. (rhs >= diemer19_interp_Gmin(n_eff)) & (rhs <= diemer19_interp_Gmax(n_eff))
    c = exp10.(diemer19_interp_Gc.(rhs, n_eff)) .* C_α
    if isscalar(mask)
        mask === false && (c=INVALID_CONCENTRATION)
    else
        c[.!mask] .= INVALID_CONCENTRATION
    end
    return c,mask
end

"""
    diemer2019_conc(c::AbstractCosmology,M200c::Union{AbstractArray,Real,u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real},statistic::String="median"; kws...)
Concentration model of Diemer & Joyce 2019. Valid `statistic`s are 'mean' and 'median'.
"""
function diemer2019_conc(c::AbstractCosmology,M200c::Union{Real,AbstractArray},z::Union{AbstractArray,Real};statistic::String="median", kws...)
    if statistic == "median"
        Κ,a0,a1,b0,b1,c_α = 0.41, 2.45, 1.82, 3.2, 2.3, 0.21
    elseif statistic == "mean"
        Κ,a0,a1,b0,b1,c_α = 0.42, 2.37, 1.74, 3.39, 1.82, 0.20
    else
        error("Statistic $statistic is not implemented in diemer2019_conc.")
    end
    return _diemer2019_general(c,M200c,z,Κ,a0,a1,b0,b1,c_α; kws...)
end
diemer2019_conc(c::AbstractCosmology,M200c::Union{u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real};kws...) = diemer2019_conc(c,M200c.|>ua.Msun.|>u.ustrip,z;kws...)



"""
    ishiyama_conc(c::AbstractCosmology,M::Union{AbstractArray,Real,u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real},mdef::String,c_type::String="fit_all"; kws...)
Concentration model of Ishiyama et al. 2020. A recalibration of the Diemer 2019 model based on the Uchuu simulation. Valid mass definitions `mdef` are '200c' and 'vir' for all samples, with the addition of '500c' for 'fit_all' and 'fit_relaxed' `c_type`. `c_type` indicates the measurement technique or sample used for the fit. Valid `c_type` are 'fit_all', 'fit_relaxed', 'vmax_all', and 'vmax_relaxed'.   
"""
function ishiyama2020_conc(c::AbstractCosmology,M::Union{Real,AbstractArray},z::Union{AbstractArray,Real},mdef::String;c_type::String="fit_all", kws...)
    if c_type == "fit_all"
	if mdef == "200c"
            Κ,a0,a1,b0,b1,c_α = 1.19, 2.54, 1.33, 4.04, 1.21, 0.22
	elseif mdef == "vir"
            Κ,a0,a1,b0,b1,c_α = 1.64, 2.67, 1.23, 3.92, 1.30, -0.19
	elseif mdef == "500c"
            Κ,a0,a1,b0,b1,c_α = 1.83, 1.95, 1.17, 3.57, 0.91, 0.26
        else
	    error("Invalid mdef $mdef for ishiyama20 model, allowed are 200c and vir for all, additionally 500c for 'fit' and 'fit_relaxed'.")
        end
    elseif c_type == "fit_relaxed"
	if mdef == "200c"
            Κ,a0,a1,b0,b1,c_α = 0.60, 2.14, 2.63, 1.69, 6.36, 0.37
	elseif mdef == "vir"
            Κ,a0,a1,b0,b1,c_α = 1.22, 2.52, 1.87, 2.13, 4.19, -0.017
	elseif mdef == "500c"
            Κ,a0,a1,b0,b1,c_α = 0.38, 1.44, 3.41, 2.86, 2.99, 0.42
        else
	    error("Invalid mdef $mdef for ishiyama20 model, allowed are 200c and vir for all, additionally 500c for 'fit' and 'fit_relaxed'.")
        end
    elseif c_type == "vmax_all"
	if mdef == "200c"
            Κ,a0,a1,b0,b1,c_α = 1.10, 2.30, 1.64, 1.72, 3.60, 0.32
        elseif mdef == "vir"
            Κ,a0,a1,b0,b1,c_α = 0.76, 2.34, 1.82, 1.83, 3.52, -0.18
        else
	    error("Invalid mdef $mdef for ishiyama20 model, allowed are 200c and vir for all, additionally 500c for 'fit' and 'fit_relaxed'.")
        end
    elseif c_type == "vmax_relaxed"
	if mdef == "200c"
            Κ,a0,a1,b0,b1,c_α = 1.79, 2.15, 2.06, 0.88, 9.24, 0.51
        elseif mdef == "vir"
            Κ,a0,a1,b0,b1,c_α = 2.40, 2.27, 1.80, 0.56, 13.24, 0.079
        else
	    error("Invalid mdef $mdef for ishiyama20 model, allowed are 200c and vir for all, additionally 500c for 'fit' and 'fit_relaxed'.")
        end
    end
    return _diemer2019_general(c,M,z,Κ,a0,a1,b0,b1,c_α; kws...)
end
ishiyama2020_conc(c::AbstractCosmology,M::Union{u.Quantity,Array{<:u.Quantity}},z::Union{AbstractArray,Real},mdef::String; kws...) = ishiyama2020_conc(c,M.|>ua.Msun.|>u.ustrip,z,mdef,kws...)



const models=Dict([("bullock01",Dict([("mdefs",["200c"]),("func",bullock2001_conc)])),
             ("duffy08",Dict([("mdefs",["200c","vir","200m"]),("func",duffy2008_conc)])),
             ("klypin11",Dict([("mdefs",["vir"]),("func",klypin2011_conc)])),
             ("prada12",Dict([("mdefs",["200c"]),("func",prada2012_conc)])),
             ("bhattacharya13",Dict([("mdefs",["200c","200m","vir"]),("func",bhattacharya2013_conc)])),
             ("dutton14",Dict([("mdefs",["200c","vir"]),("func",dutton2014_conc)])),
             ("child18",Dict([("mdefs",["200c"]),("func",child2018_conc)])),
             ("diemer19",Dict([("mdefs",["200c"]),("func",diemer2019_conc)])),
             ("ishiyama20",Dict([("mdefs",["200c","200m","500c"]),("func",ishiyama2020_conc)])),
             ])



end # module
