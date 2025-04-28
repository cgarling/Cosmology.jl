using Cosmology
using Test
# import Unitful as u # This breaks Julia 1.5 compatibility
using Unitful
const u = Unitful
using UnitfulAstro
import QuadGK: quadgk
using Documenter

# the doctest for z_at_value will throw deprecation warnings so we need
# to set depwarn=no in CI; equivalent for local testing is
# import Pkg
# import Cosmology
# Pkg.test("Cosmology";julia_args=["--depwarn=no"])
DocMeta.setdocmeta!(Cosmology, :DocTestSetup, :(using Cosmology; import Unitful; import UnitfulAstro); recursive=true)
doctest(Cosmology,manual=false)

# values from http://icosmos.co.uk/

const dist_rtol = 1e-6
const age_rtol = 2e-4
const density_rtol=1e-4
# Integrating a unitful function would require UnitfulIntegration.jl.  Without using it, we
# strip the units away from the integrand function; for comoving_volume
integrand(c, z) = 4pi*u.ustrip(comoving_volume_element(c, z))

@testset "FlatLCDM" begin
    c = cosmology(h=0.7, OmegaM=0.3, OmegaG=0, Tcmb0=0) # Tcmb0=0 is not necessary but suppresses a warning
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
    @test Cosmology.w0(c) == -1
    @test Cosmology.wa(c) == 0
    @test Cosmology.w(c,2.0) == -1
    @test Cosmology.w(2.0,Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.w(c,2.0)
    @test Cosmology.a2E(0.8,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.a2E(c,0.8)
    @test Cosmology.E(2.0,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.E(c,2.0)

    # Planck18 has neutrinos, baryons, dark matter, dark energy, and non-zero CMB temperature
    let c=Cosmology.Planck18 
        @test ρ_c(c,2.0) == ρ_c(2.0,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c))
        @test ρ_m(c,2.0) == ρ_m(2.0,Cosmology.h(c),Ω_m(c))
        @test ρ_b(c,2.0) == ρ_b(2.0,Cosmology.h(c),Ω_b(c))
        @test ρ_dm(c,2.0) == ρ_dm(2.0,Cosmology.h(c),Ω_dm(c))
        @test ρ_Λ(c,2.0) == ρ_Λ(2.0,Cosmology.h(c),Ω_Λ(c),Cosmology.w0(c),Cosmology.wa(c))
        @test ρ_γ(c,2.0) == ρ_γ(2.0,Cosmology.h(c),Ω_γ(c))
        @test ρ_ν(c,2.0) ≈ ρ_ν(2.0,Cosmology.h(c),u.ustrip(u.K,T_cmb(c)),Cosmology.Neff(c),Cosmology.m_nu(c))
        @test ρ_r(c,2.0) ≈ ρ_r(2.0,Cosmology.h(c),u.ustrip(u.K,T_cmb(c)),Cosmology.Neff(c),Cosmology.m_nu(c))
        @test Ω_m(c,2.0) == Ω_m(2.0,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c))
        @test Ω_b(c,2.0) == Ω_b(2.0,Cosmology.h(c),Ω_m(c),Ω_b(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c))
        @test Ω_dm(c,2.0) == Ω_dm(2.0,Cosmology.h(c),Ω_m(c),Ω_b(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c))

    end
    # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
    # @test ρ_c(c,0) ≈ cosmo.constants.RHO_C_Z0_CGS * c.h^2 * u.g / u.cm^3 rtol = density_rtol
    # @test nu_relative_density(c, 0) ≈ 26.65068992843329
end

@testset "OpenLCDM" begin
    c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaG=0, Tcmb0=0) # Tcmb0=0 is not necessary but suppresses a warning
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
    @test Cosmology.w0(c) == -1
    @test Cosmology.wa(c) == 0
    @test Cosmology.w(c,2.0) == -1
    @test Cosmology.w(2.0,Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.w(c,2.0)
    @test Cosmology.a2E(0.8,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.a2E(c,0.8)
    @test Cosmology.E(2.0,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.E(c,2.0)

    # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
    # @test nu_relative_density(c, 0) ≈ 26.65068992843329

end

@testset "ClosedLCDM" begin
    c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaG=0, Tcmb0=0) # Tcmb0=0 is not necessary but suppresses a warning
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
    @test Cosmology.w0(c) == -1
    @test Cosmology.wa(c) == 0
    @test Cosmology.w(c,2.0) == -1
    @test Cosmology.w(2.0,Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.w(c,2.0)
    @test Cosmology.a2E(0.8,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.a2E(c,0.8)
    @test Cosmology.E(2.0,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.E(c,2.0)

    # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
    # @test nu_relative_density(c,0) ≈ 26.65068992843329
end

@testset "FlatWCDM" begin
    c = cosmology(h=0.7, OmegaM=0.3, OmegaG=0, Tcmb0=0, w0=-0.9, wa=0.1) # Tcmb0=0 is not necessary but suppresses a warning
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
    @test Cosmology.w0(c) == -0.9
    @test Cosmology.wa(c) == 0.1
    @test Cosmology.w(c,2.0) ≈ -0.8333333333333334
    @test Cosmology.w(2.0,Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.w(c,2.0)
    @test Cosmology.de_density_scale(c,2.0) == Cosmology.de_density_scale(2.0,Cosmology.w0(c),Cosmology.wa(c))
    @test Cosmology.a2E(0.8,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.a2E(c,0.8)
    @test Cosmology.E(2.0,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.E(c,2.0)
    # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
    # @test nu_relative_density(c,0) ≈ 26.65068992843329
end

@testset "OpenWCDM" begin
    c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaG=0, Tcmb0=0, w0=-0.9, wa=0.1) # Tcmb0=0 is not necessary but suppresses a warning
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
    @test Cosmology.w0(c) == -0.9
    @test Cosmology.wa(c) == 0.1
    @test Cosmology.w(c,2.0) ≈ -0.8333333333333334
    @test Cosmology.w(2.0,Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.w(c,2.0)
    @test Cosmology.de_density_scale(c,2.0) == Cosmology.de_density_scale(2.0,Cosmology.w0(c),Cosmology.wa(c))
    @test Cosmology.a2E(0.8,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.a2E(c,0.8)
    @test Cosmology.E(2.0,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.E(c,2.0)
    # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
    # @test nu_relative_density(c,0) ≈ 26.65068992843329
end

@testset "ClosedWCDM" begin
    c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaG=0, Tcmb0=0, w0=-0.9, wa=0.1) # Tcmb0=0 is not necessary but suppresses a warning
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
    @test Cosmology.w0(c) == -0.9
    @test Cosmology.wa(c) == 0.1
    @test Cosmology.w(c,2.0) ≈ -0.8333333333333334
    @test Cosmology.w(2.0,Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.w(c,2.0)
    @test Cosmology.de_density_scale(c,2.0) == Cosmology.de_density_scale(2.0,Cosmology.w0(c),Cosmology.wa(c))
    @test Cosmology.a2E(0.8,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.a2E(c,0.8)
    @test Cosmology.E(2.0,Cosmology.h(c),Ω_m(c),Ω_k(c),Ω_Λ(c),u.ustrip(u.K,T_cmb(c)),Cosmology.m_nu(c),Cosmology.Neff(c),Cosmology.w0(c),Cosmology.wa(c)) == Cosmology.E(c,2.0)
    # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
    # @test nu_relative_density(c,0) ≈ 26.65068992843329
end

@testset "Non-Float64" begin
    # Test that FlatLCDM works with non-Float64 (BigFloat in this example)
    c = cosmology(h=0.7, OmegaM=big(0.3), OmegaG=0, Tcmb0=0) # Tcmb0=0 is not necessary but suppresses a warning
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1651.9144029437341u"Mpc" rtol = dist_rtol
    @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 625.344422840635148u"Mpc" rtol = dist_rtol
    @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
    @test comoving_volume_element(c, big(1.41)) ≈ 3.40308795871941707e10u"Mpc^3" rtol = dist_rtol
    # Test that FlatWCDM works with non-Float64 (BigFloat in this example)
    c = cosmology(h=big(0.7), OmegaM=0.3, OmegaG=0, Tcmb0=0, w0=-0.9, wa=0.1) # Tcmb0=0 is not necessary but suppresses a warning
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1612.05829248972126693u"Mpc" rtol = dist_rtol
    @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 607.6801988608u"Mpc" rtol = dist_rtol
    @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
    @test comoving_volume_element(c, big(1.41)) ≈ 3.13786257398184e10u"Mpc^3" rtol = dist_rtol
    # @test ρ_c(c,0) ≈ 9.21671792415115e-30 * u.g / u.cm^3 rtol = density_rtol
    # @test nu_relative_density(c,0) ≈ 26.65068992843329
end

@testset "Units" begin
    c = cosmology(h=0.9, OmegaM=0.5, OmegaG=0, Tcmb0=0) # Tcmb0=0 is not necessary but suppresses a warning
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
    c = cosmology(h = 0.7, N_eff=3.046, Tcmb0=2.75, OmegaK=0)
    @test hubble_time(c, 0) ≈ Cosmology.hubble_time0(c)
    @test hubble_dist(c, 0) ≈ Cosmology.hubble_dist0(c)
    @test H(c, 0) ≈ 70u"km/s/Mpc"
    @test Cosmology.n_nu(c) == 3
    @test T_nu(c,0) ≈ 1.9628561026349225 * u.K
    @test T_nu(c,1.) ≈ 3.925712205269845 * u.K
    @test T_cmb(c,0) == 2.75 * u.K
    @test T_cmb(c,1.) ≈ 5.5 * u.K
    @test z_at_value(c, T_cmb, 5.5*u.K) ≈ 1
end

