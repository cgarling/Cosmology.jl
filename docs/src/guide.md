# [Using Cosmology.jl](@id guide)

## Constructing Cosmology Instances
To create an instance of a cosmology that you can use to do computations, you can use the convenience constructor [`cosmology`](@ref Cosmology.cosmology). By default this constructor will use cosmological parameters from the Planck 2018 results; no arguments are necessary unless you want to override the defaults.

```@example guide
using Cosmology
c = cosmology()
```

Alternatively, you can use one of the [pre-constructed instances](@ref default_cosmologies), like [`Planck18`](@ref Cosmology.Planck18) or [`WMAP9`](@ref Cosmology.WMAP9). These are not exported and so must be imported as

```@example guide
import Cosmology: Planck18
Planck18
```
```@example guide
import Cosmology: WMAP9
WMAP9
```

!!! tip
    [`cosmology`](@ref Cosmology.cosmology) is type-unstable because it returns different concrete subtypes of [`AbstractCosmology`](@ref Cosmology.AbstractCosmology) depending on the parameters you give it. This results in a significant performance hit; for example, `c=cosmology()` takes ~300 ns, while `FlatLCDM(c.h,c.Ω_Λ,c.Ω_m,c.Ω_b,c.Tcmb0,c.Neff,c.m_nu)` takes ~1 ns. If you want to create many cosmology instances very quickly, it is recommended that you use the [base constructors](@ref concrete_types). However, these basic constructors do not include much in the way of argument checking; for example, see the section below on neutrino handling.

```@example guide; ansicolor = true
import BenchmarkTools: @benchmark
import Cosmology: FlatLCDM
@benchmark cosmology(h=0.6766,OmegaK=0.0,OmegaM=0.30966,OmegaB=0.04897,OmegaG=nothing,Tcmb0=2.7255,w0=-1,wa=0,N_eff=3.046,m_ν=(0.0,0.0,0.06))
```
```@example guide; ansicolor = true
@benchmark FlatLCDM($0.6766, $0.6888463055445425, $0.30966, $0.04897, $2.7255, $3.046, $(0.0,0.0,0.06))
```

## Calling Methods
With a cosmological type constructed, we can use it to call [methods](@ref methods). For example, to calculate the critical density at a redshift of 1,

```@example guide
ρ_c(c,1.0)
```

We can transform this to different units as well,

```@example guide
import Unitful as u
import UnitfulAstro as ua
ρ_c(ua.Msun/ua.kpc^3,c,1.0)
```

Methods with return values that have attached units should also have this conversion interface. See [the section on methods](@ref methods) for a full list of defined methods.

## Neutrinos
This package supports massive neutrinos by assuming that the total mass in neutrinos is split between the species. For greatest efficiency you should pass `m_ν` to `cosmology` or one of the basic constructors as an `NTuple{N,T}`; that is, a tuple where all the elements are of the same concrete type. Internally they are stored in eV but with the units stripped, so a good argument would be `m_ν=(0.0,0.0,0.06)` for 3 neutrino species with one massive species with mass 0.06 eV. You can provide neutrino masses to [`cosmology`](@ref Cosmology.cosmology) and the [basic constructors](@ref concrete_types) with energy or mass units as well;

```@example guide
cosmology(m_ν=(0.0,0.0,0.06) .* u.eV)
```
```@example guide
FlatLCDM(0.6766, 0.6888463055445425, 0.30966, 0.04897, 2.7255, 3.046, (0.0,0.0,0.06) .* u.eV)
```

There is currently some weirdness with how massless neutrinos are dealt with; for now, I recommend that you make the length of the `m_ν` iterable you provide equal to the result of `Cosmology.n_nu(N_eff)` for the effective number of neutrino species you choose. For example, if `N_eff=3.046` then `Cosmology.n_nu(N_eff)=3` and your `m_ν` iterable should have length 3. This is done for you in the [pre-constructed instances](@ref default_cosmologies). A warning will be issued by [`cosmology`](@ref Cosmology.cosmology) if `length(m_ν) != Cosmology.n_nu(N_eff)` and `!iszero(N_eff)`, but none of the basic constructors contain such checks. The implementation of massive neutrinos is open to change.

!!! tip
    Inclusion of massive neutrinos is expensive. For example, for the default massive neutrino parameters `c=cosmology()`, the evaluation of `E(c, 0.8)` takes 114.613 ns, while `E( cosmology(m_ν=(0.0,),N_eff=3.046), 0.8)` takes 6.986 ns and `E( cosmology(m_ν=(0.0,),N_eff=0), 0.8)` takes 6.095 ns. This makes a significant difference in methods that involve integrations (e.g., [`comoving_radial_dist`](@ref)). If speed is a concern, consider if you can neglect neutrinos for your calculation.

```@example guide
c_massivenu = cosmology()
c_masslessnu = cosmology(m_ν=(0.0,),N_eff=3.046)
c_nonu = cosmology(m_ν=(0.0,),N_eff=0.0)
nothing # hide
```
```@example guide; ansicolor = true
@benchmark E($c_massivenu,$0.8)
```
```@example guide; ansicolor = true
@benchmark E($c_masslessnu,$0.8)
```
```@example guide; ansicolor = true
@benchmark E($c_nonu,$0.8)
```

## Integrated Packages
Cosmology.jl is envisioned as the base for a collection of packages for things like cosmological power spectra, growth factors, halo mass functions, etc. As they become available, some examples will be shown here and on the [integrated packages](@ref integrated_packages) page.