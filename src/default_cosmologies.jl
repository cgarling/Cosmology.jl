# OmegaL = 1 - Ω_m - Ω_γ0(h,nothing,Tcmb0) * (1 + nu_relative_density(m_ν,N_eff,u.ustrip(u.K,T_nu(Tcmb0)), n_nu(N_eff)))
# Ω_γ0(h,OmegaG::Nothing,Tcmb0::Nothing) = zero(h)
# Ω_γ0(h,OmegaG::Number,Tcmb0::Nothing) = OmegaG
# Ω_γ0(h,OmegaG::Nothing,Tcmb0::Number) = 4.481620089297254e-7 * Tcmb0^4 / h^2
# function Ω_γ0(h,OmegaG::Number,Tcmb0::Number)
#     if !isapprox(OmegaG,4.481620089297254e-7 * Tcmb0^4 / h^2)
#         @warn "Both Tcmb0 and OmegaG are provided but they do not agree; OmegaG will take precedence."
#     end
#     return OmegaG
# end
# Ω_L(Ω_m,h,Tcmb0,m_ν,N_eff) = 1 - Ω_m - Ω_γ0(h,nothing,Tcmb0) * (1 + Cosmology.nu_relative_density(m_ν,N_eff,u.ustrip(u.K,T_nu(Tcmb0)), Cosmology.n_nu(N_eff)))

# Constructing these does not meaningfully affect load time
# FlatLCDM(h::Real, Ω_Λ::Real, Ω_m::Real, Ω_b::Real, Tcmb0::Real, Neff::Real, m_ν)
""" The full Planck 18 cosmology (best fit with BAO; column 6). """
const Planck18 = FlatLCDM(0.6766, 0.6888463055445425, 0.30966, 0.04897, 2.7255, 3.046, (0.0,0.0,0.06))
""" The full Planck 15 cosmology (TT,TE,EE+lowP+lensing+ext; column 6 in Table 4 of Planck 2015 paper). """
const Planck15 = FlatLCDM(0.6774,0.6896098315260938,0.3089,0.04869,2.7255,3.046,(0.0,0.0,0.06))
""" The conservative Planck 15 cosmology (TT+lowP+lensing; column 2 in Table 4 of Planck 2015 paper). """ 
const Planck15_only = FlatLCDM(0.6781,0.6905129065283897,0.308,0.04868,2.7255,3.046,(0.0,0.0,0.06))
""" The full Planck 13 cosmology. """
const Planck13 = FlatLCDM(0.6774,0.6913898315260937,0.30712,0.048252,2.7255,3.046,(0.0,0.0,0.06))
""" The full WMAP 9 cosmology. """
const WMAP9 = FlatLCDM(0.6932,0.7134367388797183,0.2865,0.04628,2.725,3.04,(0.0,))
""" The full WMAP 7 cosmology. """
const WMAP7 = FlatLCDM(0.704,0.7279386649578158,0.272,0.0455,2.725,3.04,(0.0,))
""" The full WMAP 5 cosmology. """
const WMAP5 = FlatLCDM(0.702,0.7229383149725506,0.277,0.0459,2.725,3.04,(0.0,))

# Ω_L(Ω_m,h,Tcmb0,m_ν,N_eff)
# Ω_L(0.30712,0.6774,2.7255,(0.0,0.0,0.06),3.046) # Planck13
# Ω_L(0.2865,0.6932,2.725,(0.0,),3.04) # WMAP9
# d=WMAP9; Ω_m(d) + Ω_Λ(d) + Ω_k(d) + Ω_r(d)
# Ω_L(0.272,0.704,2.725,(0.0,),3.04) # WMAP7
# d=WMAP7; Ω_m(d) + Ω_Λ(d) + Ω_k(d) + Ω_r(d)
# Ω_L(0.277,0.702,2.725,(0.0,),3.04) # WMAP5
