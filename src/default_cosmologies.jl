# we can delay evaluation of these by making them expressions and then eval'ing them when we want to use them. doesn't seem to help import time much at the moment.

# Planck18 = cosmology(h = 0.6766, OmegaK = 0,OmegaM = 0.30966,OmegaB = 0.04897,OmegaR = nothing,Tcmb0 = 2.7255*u.K,OmegaG = nothing,w0 = -1,wa = 0,Neff = 3.046,m_nu=[0, 0, 0.06] * u.eV,ns=0.9665,sigma8=0.8102,tau=0.0561,z_reion=7.82)

# Planck15 = cosmology(h = 0.6774, OmegaK = 0,OmegaM = 0.3075,OmegaB = 0.04860,OmegaR = nothing,Tcmb0 = 2.7255*u.K,OmegaG = nothing,w0 = -1,wa = 0,Neff = 3.046,m_nu=[0, 0, 0.06] * u.eV,ns=0.9667,sigma8=0.8159,tau=0.066,z_reion=8.8)

# Planck13 = cosmology(h = 0.6774, OmegaK = 0,OmegaM = 0.30712,OmegaB = 0.048252,OmegaR = nothing,Tcmb0 = 2.7255*u.K,OmegaG = nothing,w0 = -1,wa = 0,Neff = 3.046,m_nu=[0, 0, 0.06] * u.eV,ns=0.9611,sigma8=0.8288,tau=0.0952,z_reion=11.52)

# WMAP9 = cosmology(h = 0.6932, OmegaK = 0,OmegaM = 0.2865,OmegaB = 0.04628,OmegaR = nothing,Tcmb0 = 2.725*u.K,OmegaG = nothing,w0 = -1,wa = 0,Neff = 3.04,m_nu=0.,ns=0.9608,sigma8=0.820,tau=0.081,z_reion=10.1)

# WMAP7 = cosmology(h = 0.704, OmegaK = 0,OmegaM = 0.272,OmegaB = 0.0455,OmegaR = nothing,Tcmb0 = 2.725*u.K,OmegaG = nothing,w0 = -1,wa = 0,Neff = 3.04,m_nu=0.,ns=0.967,sigma8=0.810,tau=0.085,z_reion=10.3)

# WMAP5 = cosmology(h = 0.702, OmegaK = 0,OmegaM = 0.277,OmegaB = 0.0459,OmegaR = nothing,Tcmb0 = 2.725*u.K,OmegaG = nothing,w0 = -1,wa = 0,Neff = 3.04,m_nu=0.,ns=0.967,sigma8=0.817,tau=0.088,z_reion=11.3)
"""
    Planck18(;kws...)
Constructor function for the full Planck 18 cosmology (best fit with BAO; column 6). You can specify valid kws to override defaults or specify other options for the `cosmology` constructor. """
function Planck18(;kws...)
    # construct the dictionary with the default options for this cosmology
    default = Dict(:h=>0.6766, :OmegaM=>0.30966, :OmegaK=>0.0, :OmegaB=>0.04897,:sigma8=>0.8102,:ns=>0.9665,:OmegaR => nothing,:Tcmb0 => 2.7255*u.K,:OmegaG => nothing,:w0 => -1,:wa => 0,:Neff => 3.046,:m_nu=>[0, 0, 0.06] * u.eV, :tau=>0.0561,:z_reion=>7.82)
    # alter / add to the default dict for the provided kws.
    for i in eachindex(kws)
        default[i] = kws[i]
    end
    return cosmology(;default...)
end
"""
    Planck15(;kws...)
Constructor function for the full Planck 15 cosmology (TT,TE,EE+lowP+lensing+ext; column 6 in Table 4 of Planck 2015 paper.). You can specify valid kws to override defaults or specify other options for the `cosmology` constructor. """
function Planck15(;kws...)
    # construct the dictionary with the default options for this cosmology
    default = Dict(:h=>0.6774, :OmegaM=>0.3089, :OmegaK=>0.0, :OmegaB=>0.04869,:sigma8=>0.8159,:ns=>0.9667,:OmegaR => nothing,:Tcmb0 => 2.7255*u.K,:OmegaG => nothing,:w0 => -1,:wa => 0,:Neff => 3.046,:m_nu=>[0, 0, 0.06] * u.eV, :tau=>0.066,:z_reion=>8.8)
    # alter / add to the default dict for the provided kws.
    for i in eachindex(kws)
        default[i] = kws[i]
    end
    return cosmology(;default...)
end
"""
    Planck15_only(;kws...)
Constructor function for the conservative Planck 15 cosmology (TT+lowP+lensing; column 2 in Table 4 of Planck 2015 paper.). You can specify valid kws to override defaults or specify other options for the `cosmology` constructor. """
function Planck15_only(;kws...)
    # construct the dictionary with the default options for this cosmology
    default = Dict(:h=>0.6781, :OmegaM=>0.308, :OmegaK=>0.0, :OmegaB=>0.04868,:sigma8=>0.8149,:ns=>0.9677,:OmegaR => nothing,:Tcmb0 => 2.7255*u.K,:OmegaG => nothing,:w0 => -1,:wa => 0,:Neff => 3.046,:m_nu=>[0, 0, 0.06] * u.eV, :tau=>0.066,:z_reion=>8.8)
    # alter / add to the default dict for the provided kws.
    for i in eachindex(kws)
        default[i] = kws[i]
    end
    return cosmology(;default...)
end
"""
    Planck13(;kws...)
Constructor function for the full Planck 13 cosmology. You can specify valid kws to override defaults or specify other options for the `cosmology` constructor. """
function Planck13(;kws...)
    # construct the dictionary with the default options for this cosmology
    default = Dict(:h=>0.6774, :OmegaM=>0.30712, :OmegaK=>0.0, :OmegaB=>0.048252,:sigma8=>0.8288,:ns=>0.9611,:OmegaR => nothing,:Tcmb0 => 2.7255*u.K,:OmegaG => nothing,:w0 => -1,:wa => 0,:Neff => 3.046,:m_nu=>[0, 0, 0.06] * u.eV, :tau=>0.0952,:z_reion=>11.52)
    # alter / add to the default dict for the provided kws.
    for i in eachindex(kws)
        default[i] = kws[i]
    end
    return cosmology(;default...)
end
"""
    WMAP9(;kws...)
Constructor function for the full WMAP 9 cosmology. You can specify valid kws to override defaults or specify other options for the `cosmology` constructor. """
function WMAP9(;kws...)
    # construct the dictionary with the default options for this cosmology
    default = Dict(:h=>0.6932, :OmegaM=>0.2865, :OmegaK=>0.0, :OmegaB=>0.04628,:sigma8=>0.820,:ns=>0.9608,:OmegaR => nothing,:Tcmb0 => 2.725*u.K,:OmegaG => nothing,:w0 => -1,:wa => 0,:Neff => 3.04,:m_nu=>0.0, :tau=>0.081,:z_reion=>10.1)
    # alter / add to the default dict for the provided kws.
    for i in eachindex(kws)
        default[i] = kws[i]
    end
    return cosmology(;default...)
end
"""
    WMAP7(;kws...)
Constructor function for the full WMAP 7 cosmology. You can specify valid kws to override defaults or specify other options for the `cosmology` constructor. """
function WMAP7(;kws...)
    # construct the dictionary with the default options for this cosmology
    default = Dict(:h=>0.704, :OmegaM=>0.272, :OmegaK=>0.0, :OmegaB=>0.0455,:sigma8=>0.810,:ns=>0.967,:OmegaR => nothing,:Tcmb0 => 2.725*u.K,:OmegaG => nothing,:w0 => -1,:wa => 0,:Neff => 3.04,:m_nu=>0.0, :tau=>0.085,:z_reion=>10.3)
    # alter / add to the default dict for the provided kws.
    for i in eachindex(kws)
        default[i] = kws[i]
    end
    return cosmology(;default...)
end
"""
    WMAP5(;kws...)
Constructor function for the full WMAP 5 cosmology. You can specify valid kws to override defaults or specify other options for the `cosmology` constructor. """
function WMAP5(;kws...)
    # construct the dictionary with the default options for this cosmology
    default = Dict(:h=>0.702, :OmegaM=>0.277, :OmegaK=>0.0, :OmegaB=>0.0459,:sigma8=>0.817,:ns=>0.967,:OmegaR => nothing,:Tcmb0 => 2.725*u.K,:OmegaG => nothing,:w0 => -1,:wa => 0,:Neff => 3.04,:m_nu=>0.0, :tau=>0.088,:z_reion=>11.3)
    # alter / add to the default dict for the provided kws.
    for i in eachindex(kws)
        default[i] = kws[i]
    end
    return cosmology(;default...)
end
