using Cosmology
# using Cosmology.halo
using Test
import Unitful as u
using Unitful
using UnitfulAstro
import QuadGK: quadgk

# values from http://icosmos.co.uk/

const dist_rtol = 1e-6
const age_rtol = 2e-4
const density_rtol=1e-4
# Integrating a unitful function would require UnitfulIntegration.jl.  Without using it, we
# strip the units away from the integrand function
integrand(c, z) = 4pi*u.ustrip(comoving_volume_element(c, z))

@testset verbose = true "cosmo" begin
    @testset "FlatLCDM" begin
        c = cosmology(h=0.7, OmegaM=0.3, OmegaG=0)
        @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1651.9144029437339u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 625.3444228406352u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
        @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3303.8288058874678u"Mpc" rtol = dist_rtol
        @test comoving_volume(c,1,rtol=dist_rtol) ≈ 151.05712532061932u"Gpc^3" rtol = dist_rtol
        @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ u.ustrip(comoving_volume(c, 2.5))
        @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6607.6576117749355u"Mpc" rtol = dist_rtol
        @test distmod(c,1,rtol=dist_rtol) ≈ 44.10023765554372 rtol = dist_rtol
        @test age(c,0,rtol=age_rtol) ≈ 13.467105300439439u"Gyr" rtol = age_rtol
        @test age(c,1,rtol=age_rtol) ≈ 5.751689848348662u"Gyr" rtol = age_rtol
        @test lookback_time(c,1,rtol=age_rtol) ≈ 7.715337003613595u"Gyr" rtol = age_rtol
        @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0) rtol = age_rtol
        # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
        # @test ρ_c(c,0) ≈ cosmo.constants.RHO_C_Z0_CGS * c.h^2 * u.g / u.cm^3 rtol = density_rtol
        # @test nu_relative_density(c, 0) ≈ 26.65068992843329
    end

    @testset "OpenLCDM" begin
        c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaG=0)
        @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1619.9586273816346u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 598.9117950429865u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
        @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3209.7837001394314u"Mpc" rtol = dist_rtol
        @test comoving_volume(c,1,rtol=dist_rtol) ≈ 140.0855965772839u"Gpc^3" rtol = dist_rtol
        @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ u.ustrip(comoving_volume(c, 2.5))
        @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6479.834509526539u"Mpc" rtol = dist_rtol
        @test distmod(c,1,rtol=dist_rtol) ≈ 44.057819572163766 rtol = dist_rtol
        @test age(c,0,rtol=age_rtol) ≈ 13.061816787881945u"Gyr" rtol = age_rtol
        @test age(c,1,rtol=age_rtol) ≈ 5.545627298881875u"Gyr" rtol = age_rtol
        @test lookback_time(c,1,rtol=age_rtol) ≈ 7.5161110299320155u"Gyr" rtol = age_rtol
        @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0) rtol = age_rtol
        # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
        # @test nu_relative_density(c, 0) ≈ 26.65068992843329

    end

    @testset "ClosedLCDM" begin
        c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaG=0)
        @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1686.5271861439733u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 655.6019184358607u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
        @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3408.9370198733986u"Mpc" rtol = dist_rtol
        @test comoving_volume(c,1,rtol=dist_rtol) ≈ 163.84786905899605u"Gpc^3" rtol = dist_rtol
        @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ u.ustrip(comoving_volume(c, 2.5))
        @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6746.108744575893u"Mpc" rtol = dist_rtol
        @test distmod(c,1,rtol=dist_rtol) ≈ 44.14526668781513 rtol = dist_rtol
        @test age(c,0,rtol=age_rtol) ≈ 13.92262503761163u"Gyr" rtol = age_rtol
        @test age(c,1,rtol=age_rtol) ≈ 5.985745756567641u"Gyr" rtol = age_rtol
        @test lookback_time(c,1,rtol=age_rtol) ≈ 7.936800843174262u"Gyr" rtol = age_rtol
        @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0) rtol = age_rtol
        # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
        # @test nu_relative_density(c,0) ≈ 26.65068992843329
    end

    @testset "FlatWCDM" begin
        c = cosmology(h=0.7, OmegaM=0.3, OmegaG=0, w0=-0.9, wa=0.1)
        @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1612.0582924897212u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 607.6801988608612u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
        @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3224.1165849794425u"Mpc" rtol = dist_rtol
        @test comoving_volume(c,1,rtol=dist_rtol) ≈ 140.38503745162228u"Gpc^3" rtol = dist_rtol
        @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ u.ustrip(comoving_volume(c, 2.5))
        @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6448.233169958885u"Mpc" rtol = dist_rtol
        @test distmod(c,1,rtol=dist_rtol) ≈ 44.04720366646934 rtol = dist_rtol
        @test age(c,0,rtol=age_rtol) ≈ 13.189290392581498u"Gyr" rtol = age_rtol
        @test age(c,1,rtol=age_rtol) ≈ 5.645434355117561u"Gyr" rtol = age_rtol
        @test lookback_time(c,1,rtol=age_rtol) ≈ 7.543777588964826u"Gyr" rtol = age_rtol
        @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0) rtol = age_rtol
        # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
        # @test nu_relative_density(c,0) ≈ 26.65068992843329
    end

    @testset "OpenWCDM" begin
        c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaG=0, w0=-0.9, wa=0.1)
        @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1588.0178528134486u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 585.4929230997847u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
        @test comoving_radial_dist(c,rtol=dist_rtol,1) ≈ 3147.622246929028u"Mpc" rtol = dist_rtol
        @test comoving_volume(c,1,rtol=dist_rtol) ≈ 132.04651936026573u"Gpc^3" rtol = dist_rtol
        @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ u.ustrip(comoving_volume(c, 2.5))
        @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6352.071411253794u"Mpc" rtol = dist_rtol
        @test distmod(c,1,rtol=dist_rtol) ≈ 44.01457685935837 rtol = dist_rtol
        @test age(c,0,rtol=age_rtol) ≈ 12.846672864196098u"Gyr" rtol = age_rtol
        @test age(c,1,rtol=age_rtol) ≈ 5.464935753468967u"Gyr" rtol = age_rtol
        @test lookback_time(c,1,rtol=age_rtol) ≈ 7.3816586516400555u"Gyr" rtol = age_rtol
        @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0) rtol = age_rtol
        # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
        # @test nu_relative_density(c,0) ≈ 26.65068992843329
    end

    @testset "ClosedWCDM" begin
        c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaG=0, w0=-0.9, wa=0.1)
        @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1637.5992334783048u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 632.5829291313883u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
        @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3307.9930096427556u"Mpc" rtol = dist_rtol
        @test comoving_volume(c,1,rtol=dist_rtol) ≈ 149.83003658984978u"Gpc^3" rtol = dist_rtol
        @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ u.ustrip(comoving_volume(c, 2.5))
        @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6550.396933913219u"Mpc" rtol = dist_rtol
        @test distmod(c,1,rtol=dist_rtol) ≈ 44.08133808849712 rtol = dist_rtol
        @test age(c,0,rtol=age_rtol) ≈ 13.567944716196108u"Gyr" rtol = age_rtol
        @test age(c,1,rtol=age_rtol) ≈ 5.847135987422405u"Gyr" rtol = age_rtol
        @test lookback_time(c,1,rtol=age_rtol) ≈ 7.720730290879207u"Gyr" rtol = age_rtol
        @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0) rtol = age_rtol
        # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
        # @test nu_relative_density(c,0) ≈ 26.65068992843329
    end

    @testset "Non-Float64" begin
        # Test that FlatLCDM works with non-Float64 (BigFloat in this example)
        c = cosmology(h=0.7, OmegaM=big(0.3), OmegaG=0)
        @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1651.9144029437341u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 625.344422840635148u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
        @test comoving_volume_element(c, big(1.41)) ≈ 3.40308795871941707e10u"Mpc^3" rtol = dist_rtol
        # Test that FlatWCDM works with non-Float64 (BigFloat in this example)
        c = cosmology(h=big(0.7), OmegaM=0.3, OmegaG=0, w0=-0.9, wa=0.1)
        @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1612.05829248972126693u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 607.6801988608u"Mpc" rtol = dist_rtol
        @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
        @test comoving_volume_element(c, big(1.41)) ≈ 3.13786257398184e10u"Mpc^3" rtol = dist_rtol
        # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
        # @test nu_relative_density(c,0) ≈ 26.65068992843329
    end

    @testset "Units" begin
        c = cosmology(h=0.9, OmegaM=0.5, OmegaG=0)
        for u in (u"m", u"pc", u"ly")
            @test unit(luminosity_dist(u, c, 1)) == u
            @test unit(angular_diameter_dist(u, c, 2)) == u
            @test unit(hubble_dist(u,c,1)) == u
            @test unit(comoving_radial_dist(u,c,1)) == u
            @test unit(comoving_transverse_dist(u,c,1)) == u
            @test unit(angular_diameter_dist(u,c,1)) == u
            @test unit(sound_horizon(u,c)) == u
        end
        for u in (u"s", u"yr")
            @test unit(age(u, c, 3)) == u
            @test unit(lookback_time(u, c, 4)) == u
            @test unit(hubble_time(u, c, 4)) == u
        end
        for u in (u"m^3", u"km^3", u"ly^3")
            @test unit(comoving_volume(u, c, 1)) == u
            @test unit(comoving_volume_element(u, c, 1)) == u
        end
        for u in (u"°C", u"°F")
            @test unit(T_nu(u, c, 1)) == u
            @test unit(T_cmb(u, c, 1)) == u
        end
        for u in (u"kg/m^3", u"g/m^3")
            @test unit(ρ_c(u,c,1)) == u
            @test unit(ρ_m(u,c,1)) == u
            @test unit(ρ_b(u,c,1)) == u
            @test unit(ρ_dm(u,c,1)) == u
            @test unit(ρ_Λ(u,c,1)) == u
            @test unit(ρ_γ(u,c,1)) == u
            @test unit(ρ_r(u,c,1)) == u
        end
    end

    @testset "Utilities" begin
        c = cosmology(h = 0.7,N_eff=3.046,Tcmb0=2.75,OmegaK=0)
        @test hubble_time(c, 0) ≈ Cosmology.hubble_time0(c)
        @test hubble_dist(c, 0) ≈ Cosmology.hubble_dist0(c)
        @test H(c, 0) ≈ 70u"km/s/Mpc"
        @test n_nu(c) == 3
        @test T_nu(c,0) ≈ 1.9628561026349225 * u.K
        @test T_nu(c,1.) ≈ 3.925712205269845 * u.K
        @test T_cmb(c,0) == 2.75 * u.K
        @test T_cmb(c,1.) ≈ 5.5 * u.K
        @test z_at_value(c, scale_factor, 0.8) ≈ 0.25
        @test Cosmology.isscalar(1) == true
        @test Cosmology.isscalar([1,2]) == false
        @test Cosmology.isscalar(1*u.s) == true
        @test length(Cosmology.logspace(1e-2,1e2,length=100)) == 100
        @test (x=Cosmology.logspace(1e-2,1e2,step=0.1); log10(x[2])-log10(x[1])) ≈ 0.1
    end

    # @testset "Halo" begin
    #     c = cosmology(h=0.6766,OmegaK=0,OmegaM=0.30966,OmegaB=0.04897,OmegaR=nothing,w0=-1,wa=0,Neff=3.046,m_nu=[0, 0, 0.06] * u.eV,ns=0.9665,sigma8=0.8102,tau=0.0561,z_reion=7.82,a_min=1e-5,transfer_model=Eisenstein1998,k_min=1e-20,k_max=1e20,dlogk=0.01)
    #     b = NFW(c=c,M=1e12,conc=10,z=0.0,mdef="200c")
    #     @test b.rs ≈ 21.10080548234697
    #     @test b.rhos ≈ 5.689255015396299e6
    #     @test Δvir(c,0.) ≈ 102.45879591120845
    #     @test ρ(1.0,b) ≈ 1.094299746287284e8
    #     @test ∇ρ(1.0,b) ≈ -1.1933277882436149e8
    #     @test ρmean(1.0,b) ≈ 1.692707127584995e8
    #     @test ∇ρmean(1.0,b) ≈ -1.7952221438931662e8
    #     @test enclosed_mass(1.0,b) ≈ 7.090395035600148e8
    #     @test all(Σ([1.0,b.rs,50.0],b) .≈ [6.608123771664459e8, 8.003190894622941e7, 2.453549690288129e7])
    #     @test all(projected_mass([1.0,b.rs,50.0],b) .≈
    #         [2.449752480502437e9, 2.0610687805610764e11, 4.688096750221652e11])
    #     @test Φ(1.0,b) ≈ -133761.4882774559 # * (u.km/u.s)^2 #might not include units in the future
    #     @test ∇Φ(1.0,b) ≈ 9.882822055975421e-14
    #     @test ∇∇Φ(1.0,b) ≈ -1.9397068199885201e-31
    #     @test ∇²Φ(1.0,b) ≈ 6.2116383260549736e-30
    #     @test circular_velocity(1.0,b) ≈ 55.22246142648133
    #     @test escape_velocity(1.0,b) ≈ 517.226233436503
    # end
end
