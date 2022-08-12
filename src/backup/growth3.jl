# import OrdinaryDiffEq: ODEProblem, solve, Tsit5
# import Dierckx: Spline1D # the odesolution from OrdinaryDiffEq is annoying so we're going to construct a new one
"""
    setup_growth(c::AbstractCosmology,a1::Real=1e-5,a2::Real=1.0)
    setup_growth(Ωm::Real,ΩΛ::Real,Ωk::Real,Ωg::Real,Neff::Real,m_nu::Union{Real,AbstractArray},Tcmb0::Real,w0::Real,wa::Real,a1::Real=1e-5,a2::Real=1.0)
Function that computes the linear growth factor as Equation 11 in Linder 2003, modelled after the implementation in Colossus and agrees to a few percent at the limit of Colossus' redshift range. This version is much better than the previous one, but still diverges from the Colossus result at high redshift; think it's fine for now. Although it goes negative at high redshift, need to fix that ...
"""
function setup_growth(c::AbstractCosmology,a1::Real=1e-5,a2::Real=1.0) # Union{FlatWCDM,OpenWCDM,ClosedWCDM}
    # Ok = Ω_k(c,0)# 1 - Ωm - ΩΛ
    z1=1/a1-1
    u0 = [1.0; 0.0]
    tspan = (a1, a2)
    :wa in keys(typeof(c)) ? wa=c.wa : wa=0
    :w0 in keys(typeof(c)) ? w0=c.w0 : w0=-1
    prob = ODEProblem(growth!, u0, tspan, [c.Ω_m, c.Ω_Λ, c.Ω_k, c.Ω_γ, c.Neff, c.m_nu, c.Tcmb0, w0, wa])
    sol = solve(prob, Tsit5())
    tmp_sol = transpose(reduce(hcat,sol.u))
    growth_f = tmp_sol[:,1] .* sol.t
    # z_arr= @. 1.0/sol.t - 1.0
    # interp1 = Spline1D(reverse(z_arr), reverse(growth_f./maximum(growth_f)); k=4, bc="error", s=0.0)
    # interp2 = Spline1D(growth_f./maximum(growth_f),z_arr; k=4, bc="nearest", s=0.0) #invert the growth function
    interp1 = Spline1D(sol.t, growth_f./maximum(growth_f); k=4, bc="error", s=0.0)
    interp2 = Spline1D(growth_f./maximum(growth_f),sol.t; k=4, bc="nearest", s=0.0)
    return [interp1,interp2]
end
function setup_growth(Ωm::Real,ΩΛ::Real,Ωk::Real,Ωg::Real,Neff::Real,m_nu::Union{Real,AbstractArray},Tcmb0::Real,w0::Real,wa::Real,a1::Real=1e-5,a2::Real=1.0)
    z1=1/a1-1
    u0 = [1.0; 0.0]
    tspan = (a1, a2)
    prob = ODEProblem(growth!, u0, tspan, [Ωm, ΩΛ, Ωk, Ωg, Neff, m_nu, Tcmb0, w0, wa])
    sol = solve(prob, Tsit5())
    # tmp_sol = transpose(reduce(hcat,sol.u))
    # growth_f = tmp_sol[:,1] .* sol.t
    # z_arr= @. 1.0/sol.t - 1.0
    
    # a_arr = range(a1,a2,length=1000)
    # tmp_sol = transpose(reduce(hcat,sol(a_arr).u))
    # growth_f = tmp_sol[:,1] .* a_arr
    # z_arr = @. 1/a_arr - 1
    # interp1 = Spline1D(reverse(z_arr), reverse(growth_f./maximum(growth_f)); k=4, bc="error", s=0.0)
    # interp2 = Spline1D(growth_f./maximum(growth_f),z_arr; k=4, bc="nearest", s=0.0)#invert the growth function
    tmp_sol = transpose(reduce(hcat,sol.u))
    growth_f = tmp_sol[:,1] .* sol.t
    interp1 = Spline1D(sol.t, growth_f./maximum(growth_f); k=4, bc="error", s=0.0)
    interp2 = Spline1D(growth_f./maximum(growth_f),sol.t; k=4, bc="nearest", s=0.0)
    return [interp1,interp2]
end

"""
Differential equation to solve for the growth function. Following Equation 11 in Linder 2003. 
"""
function growth!(du, U, p, a)
    G, δG = U
    Ωm0, ΩΛ0, Ωk0, Ωg0, Neff, m_nu, Tcmb0, w0, wa0 = p
    z=1/a-1
    Or = Ωg0 * (1+nu_relative_density(m_nu,Neff,Int(round(Neff)),u.ustrip(T_nu(Tcmb0,z))))
    wa = wa0 * z / (1.0 + z)
    ade = exp((1 - 3 * (w0 + wa)) * log(a) + 3 * wa * (a - 1))
    E = sqrt(Or + (Ωm0 + Ωk0 * a) * a + ΩΛ0 * ade)

    zp1 = 1+z
    de_density_scale = zp1^(3 * (1 + w0 + wa)) * exp(-3 * wa * z / zp1)
    X = (Ωm0 * (1+z)^3 / E^2) / (ΩΛ0 * de_density_scale / E^2)
    t1 = -(3.5 - 1.5 * (wa+w0) / (1.0 + X)) / a
    t2 = -1.5 * (1.0 - (wa+w0)) / (1.0 + X) / a^2

    du[1] = δG
    du[2] = t1*δG + t2 * G
end


##### integral form for ΛCDM ################################################
""" 
     growth_integral(c::Union{FlatLCDM,OpenLCDM,ClosedLCDM},z::Real,normalize=1)
Convenience implementation for calculating the linear growth factor for lambda CDM cosmologies using Equation 8 in Eisenstein1999. Generally less accurate and takes longer to run than the ODE solution. Not used for any internal calculations. """
function growth_integral(c::Union{FlatLCDM,OpenLCDM,ClosedLCDM},z::Real,normalize=1)
    Om, Ok, Ol, h, Tcmb0 = c.Ω_m, c.Ω_k, c.Ω_Λ, c.h, c.Tcmb0
    zeq = 2.5e4 * Om * h^2 * (Tcmb0/2.7)^-4
    function Ez_D(z1)
	ai = (1.0 + z1)
        E = sqrt( Om * ai^3 + Ok * ai^2 + 
            Ol * de_density_scale(c,z1) )
    end
    return 5 * Om / 2 * (1+zeq) * Ez_D(z) * quadgk(x-> (1+x) / Ez_D(x)^3, z, Inf)[1] / normalize
end
function growth_integral(z::Real, h::Real, Ωm::Real, Ωk::Real, ΩΛ::Real, Tcmb0::Union{Real,u.Quantity}, w0::Real=-1.0, wa::Real=0.0, normalize=1)
    Tcmb0 isa u.Quantity && (Tcmb0=u.ustrip(Tcmb0|>u.K))
    zeq = 2.5e4 * Ωm * h^2 * (Tcmb0/2.7)^-4
    function Ez_D(z1)
	ai = (1.0 + z1)
        E = sqrt( Ωm * ai^3 + Ωk * ai^2 + 
            ΩΛ * de_density_scale(z1,w0,wa) )
    end
    return 5 * Ωm / 2 * (1+zeq) * Ez_D(z) * quadgk(x-> (1+x) / Ez_D(x)^3, z, Inf)[1] / normalize
end
