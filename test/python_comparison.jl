using cosmo
using cosmo.halo
using cosmo.hmf
using utils
using PyCall
import PyPlot as plt

### power spectrum stuff
# colossus_cosmo = pyimport("colossus.cosmology.cosmology")
# colossus_obj = colossus_cosmo.setCosmology("planck18")

# cosmo_obj = cosmology(transfer_model="eisenstein98")
# k_arr = logspace(1e-3,1e3,length=1000)
# doing this with additional comparison to HMF.
# Pk1 = cosmo_obj.Pk_interp(k_arr)
# Pk2 = colossus_obj.matterPowerSpectrum(k_arr,z=0.0,model="eisenstein98")

# fig=plt.figure()
# ax1=plt.subplot(111)
# ax1.plot(k_arr,Pk1,label="cosmo E98")
# ax1.plot(k_arr,Pk2,label="Colossus E98",ls="--")
# ax1.set_xscale("log")
# ax1.set_yscale("log")
# ax1.legend()
# plt.show()

## hmf stuff
colossus_cosmo = pyimport("colossus.cosmology.cosmology")
colossus_obj = colossus_cosmo.setCosmology("planck15")
astropy = pyimport("astropy.cosmology")
astropy_planck15 = astropy.Planck15
py_hmf = pyimport("hmf")
camb = pyimport("camb")
cosmo_obj = cosmo.Planck15(;transfer_model=Eisenstein1998)# cosmology(transfer_model="eisenstein98")
z_mf1=0.0
# z_mf2=1.0

py_mf = py_hmf.hmf.MassFunction(Mmin=5.,Mmax=20.,dlog10m=.01, hmf_model=py_hmf.fitting_functions.SMT,hmf_params=Dict([("cosmo",astropy_planck15)]),delta_h=200.0,delta_wrt="crit",delta_c=δc(cosmo_obj,z_mf1),filter_model=py_hmf.filters.TopHat,transfer_model=py_hmf.transfer_models.EH_BAO,transfer_params=nothing,n=cosmo_obj.ns, sigma_8=cosmo_obj.sigma8,takahashi=false,lnk_min=-18.420680743952367,lnk_max=9.9034875525361272,dlnk=0.05,growth_model=py_hmf.growth_factor.GrowthFactor,growth_params=nothing,z=z_mf1)
# py_mf = hmf.hmf.MassFunction(Mmin=5.,Mmax=20.,dlog10m=.01, hmf_model=hmf.fitting_functions.Jenkins,hmf_params={"cosmo":Planck15},delta_h=200.0,delta_wrt="crit",delta_c=1.686,filter_model=hmf.filters.TopHat,transfer_model=hmf.transfer_models.EH_BAO,transfer_params=None,n=0.9665, sigma_8=0.8102,takahashi=False,lnk_min=-18.420680743952367,lnk_max=9.9034875525361272,dlnk=0.05,growth_model=hmf.growth_factor.GrowthFactor,growth_params=None,z=0.0)

camb_params=camb.CAMBparams(H0=cosmo_obj.h*100,ombh2=cosmo_obj.Ω_b*cosmo_obj.h^2,omch2=cosmo_obj.Ω_dm * cosmo_obj.h^2, omk=cosmo_obj.Ω_k, nnu=cosmo_obj.Neff, standard_neutrino_neff=cosmo_obj.Neff, num_massive_neutrinos=length(cosmo_obj.m_nu[cosmo_obj.m_nu.!=0]), m_nu=sum(cosmo_obj.m_nu),TCMB=cosmo_obj.Tcmb0,InitPower=camb.initialpower.InitialPowerLaw(ns=cosmo_obj.ns),WantTransfer=true, Transfer=camb.model.TransferParams(high_precision=true, accurate_massive_neutrinos=true,kmin=1e-4,kmax=1e3,PK_num_redshifts=1,PK_redshifts=[z_mf1],npoints=1000),WantCMB=false,WantScalars=false,WantDerivedParameters=false,Want_cl_2D_array=false,Want_CMB_lensing=false,DoLensing=false)
camb_result = camb.get_results(camb_params)
camb_sigma8 = camb_result.get_sigma8()[1]
# camb_k, camb_z, camb_Pk = camb_result.get_matter_power_spectrum(minkh=1e-10,maxkh=1e10,npoints=1000)
# run again to renormalize sigma8. Key is the As constant in the call to InitialPowerLaw
# actually, just renormalize the initial result
# camb_params=camb.CAMBparams(H0=cosmo_obj.h*100,ombh2=cosmo_obj.Ω_b*cosmo_obj.h^2,omch2=cosmo_obj.Ω_dm * cosmo_obj.h^2, omk=cosmo_obj.Ω_k, nnu=cosmo_obj.Neff, standard_neutrino_neff=cosmo_obj.Neff, num_massive_neutrinos=length(cosmo_obj.m_nu[cosmo_obj.m_nu.!=0]), m_nu=sum(cosmo_obj.m_nu),TCMB=cosmo_obj.Tcmb0,InitPower=camb.initialpower.InitialPowerLaw(As=camb_params.InitPower.As*(cosmo_obj.sigma8/camb_sigma8)^2,ns=cosmo_obj.ns),WantTransfer=true, Transfer=camb.model.TransferParams(high_precision=true, accurate_massive_neutrinos=false,kmin=1e-4,kmax=1e3,PK_num_redshifts=1,PK_redshifts=[z_mf1],npoints=1000),WantCMB=false,WantScalars=false,WantDerivedParameters=false,Want_cl_2D_array=false,Want_CMB_lensing=false,DoLensing=false)
# camb_result = camb.get_results(camb_params)
# camb_sigma8 = camb_result.get_sigma8()[1]
camb_k, camb_z, camb_Pk = camb_result.get_matter_power_spectrum(minkh=1e-4,maxkh=1e3,npoints=1000)
camb_Pk = @. camb_Pk * cosmo_obj.sigma8^2 / camb_sigma8^2

k_arr = py_mf.k
Pk1 = cosmo_obj.Pk_interp(k_arr)
T1 = transfer_function(k_arr,cosmo_obj,Eisenstein1998)
Pk2 = colossus_obj.matterPowerSpectrum(k_arr,z=0.0,model="eisenstein98")
T2 = colossus_obj.transferFunction(k_arr)
T3 = py_mf.transfer_function ./ maximum(py_mf.transfer_function)

fig,ax1=plt.subplots()
ax1.plot(k_arr,T3,label="HMF Transfer")
ax1.plot(k_arr,T1,label="My P(k)", ls="--")
ax1.plot(k_arr,T2,label="Colossus P(k)", ls="dotted")
ax1.set_xscale("log")
ax1.legend()
plt.show()

fig,ax1=plt.subplots()
ax1.plot(py_mf.k,py_mf.power,label="HMF P(k)")
ax1.plot(k_arr,Pk1,label="My P(k)", ls="--")
ax1.plot(k_arr,Pk2,label="Colossus P(k)", ls="dotted")
ax1.plot(camb_k,camb_Pk',label="CAMB", c="k")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.legend()
plt.show()

my_mf = massfunc_dndm(py_mf.m,cosmo_obj,Sheth2001,z_mf1)

fig=plt.figure()
ax1=plt.subplot(111)
ax1.plot(py_mf.m,py_mf.dndm,label="Python HMF")
ax1.plot(py_mf.m,my_mf,label="My HMF same z",ls="--")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.legend()
plt.show()

# we get agreement with python hmf if we divide rho_mean by cosmo.h^2 in hmf.jl
fig,ax1=plt.subplots()
ax1.plot(py_mf.m,(py_mf.dndm .- (my_mf)) ./ ((my_mf.+py_mf.dndm)./2),label="Fractional Difference")
ax1.set_xscale("log")
plt.show()

# growth factors
z_arr=collect(0:0.001:499)
fig=plt.figure()
ax1=plt.subplot(111)
ax1.plot(z_arr,colossus_obj.growthFactor(z_arr),label="Colossus")
ax1.plot(z_arr,cosmo_obj.growth_function(z_arr),ls="--",label="Mine")
ax1.set_ylabel("D")
ax1.set_xlabel("z")
ax1.set_xscale("log")
ax1.legend()
plt.show()

fig=plt.figure()
ax1=plt.subplot(111)
ax1.plot(z_arr,colossus_obj.growthFactor(z_arr) ./ cosmo_obj.growth_function(z_arr) .- 1)
ax1.set_xlabel("z")
ax1.set_ylabel("(Colossus growth / my growth) - 1")
ax1.set_xscale("log")
plt.show()


# look at HMcode2020 nonlinear model
camb_params2=camb.CAMBparams(H0=cosmo_obj.h*100,ombh2=cosmo_obj.Ω_b*cosmo_obj.h^2,omch2=cosmo_obj.Ω_dm * cosmo_obj.h^2, omk=cosmo_obj.Ω_k, nnu=cosmo_obj.Neff, standard_neutrino_neff=cosmo_obj.Neff, num_massive_neutrinos=length(cosmo_obj.m_nu[cosmo_obj.m_nu.!=0]), m_nu=sum(cosmo_obj.m_nu),TCMB=cosmo_obj.Tcmb0,InitPower=camb.initialpower.InitialPowerLaw(ns=cosmo_obj.ns),WantTransfer=true, Transfer=camb.model.TransferParams(high_precision=true, accurate_massive_neutrinos=true,kmin=1e-4,kmax=1e3,PK_num_redshifts=1,PK_redshifts=[z_mf1],npoints=1000),WantCMB=false,WantScalars=false,WantDerivedParameters=false,Want_cl_2D_array=false,Want_CMB_lensing=false,DoLensing=false,NonLinear="NonLinear_pk",NonLinearModel=camb.nonlinear.Halofit(halofit_version="mead2020"))
camb_result2 = camb.get_results(camb_params2)
camb_sigma82 = camb_result2.get_sigma8()[1]
# camb_k, camb_z, camb_Pk = camb_result.get_matter_power_spectrum(minkh=1e-10,maxkh=1e10,npoints=1000)
# run again to renormalize sigma8. Key is the As constant in the call to InitialPowerLaw
# actually, just renormalize the initial result
# camb_params=camb.CAMBparams(H0=cosmo_obj.h*100,ombh2=cosmo_obj.Ω_b*cosmo_obj.h^2,omch2=cosmo_obj.Ω_dm * cosmo_obj.h^2, omk=cosmo_obj.Ω_k, nnu=cosmo_obj.Neff, standard_neutrino_neff=cosmo_obj.Neff, num_massive_neutrinos=length(cosmo_obj.m_nu[cosmo_obj.m_nu.!=0]), m_nu=sum(cosmo_obj.m_nu),TCMB=cosmo_obj.Tcmb0,InitPower=camb.initialpower.InitialPowerLaw(As=camb_params.InitPower.As*(cosmo_obj.sigma8/camb_sigma8)^2,ns=cosmo_obj.ns),WantTransfer=true, Transfer=camb.model.TransferParams(high_precision=true, accurate_massive_neutrinos=false,kmin=1e-4,kmax=1e3,PK_num_redshifts=1,PK_redshifts=[z_mf1],npoints=1000),WantCMB=false,WantScalars=false,WantDerivedParameters=false,Want_cl_2D_array=false,Want_CMB_lensing=false,DoLensing=false)
# camb_result = camb.get_results(camb_params)
# camb_sigma8 = camb_result.get_sigma8()[1]
camb_k2, camb_z2, camb_Pk2 = camb_result2.get_matter_power_spectrum(minkh=1e-4,maxkh=1e3,npoints=1000)
camb_Pk2 = @. camb_Pk2 * cosmo_obj.sigma8^2 / camb_sigma82^2

fig=plt.figure()
ax1=plt.subplot(111)
ax1.plot(camb_k,camb_Pk',label="Linear CAMB")
ax1.plot(camb_k2,camb_Pk2',label="Nonlinear CAMB + Mead",ls="--")
# ax1.plot(camb_k2,camb_Pk' ./ camb_Pk2',label="Nonlinear Difference",ls="--")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.legend()
ax1.set_ylabel("P(k)")
ax1.set_xlabel("k")
plt.show(fig)

                             
