# [Defined Types](@id types)
## Abstract Types
The highest level type defined in this package is the abstract type [`Cosmology.AbstractCosmology`](@ref Cosmology.AbstractCosmology), which all other cosmology objects are subtyped from. Below these are [`Cosmology.AbstractFlatCosmology`](@ref Cosmology.AbstractFlatCosmology), [`Cosmology.AbstractOpenCosmology`](@ref Cosmology.AbstractOpenCosmology), and [`Cosmology.AbstractClosedCosmology`](@ref Cosmology.AbstractClosedCosmology).
```@docs
Cosmology.AbstractCosmology
Cosmology.AbstractFlatCosmology
Cosmology.AbstractOpenCosmology
Cosmology.AbstractClosedCosmology
```

## [Concrete Types](@id concrete_types)
The following concrete types, representing specific types of cosmologies, are currently implemented:
```@docs
Cosmology.FlatLCDM
Cosmology.OpenLCDM
Cosmology.ClosedLCDM
Cosmology.FlatWCDM
Cosmology.OpenWCDM
Cosmology.ClosedWCDM
```
 
## Convenience Constructors
```@docs
Cosmology.cosmology
Cosmology.WCDM
```

## [Pre-Constructed Instances](@id default_cosmologies)
These are constant instances that implement published cosmological results.
!!! note
    These instances are not exported, so you must import them explicitly; e.g,

```@example
import Cosmology: Planck18
Planck18
```

```@docs
Cosmology.Planck18
Cosmology.Planck15
Cosmology.Planck15_only
Cosmology.Planck13
Cosmology.WMAP9
Cosmology.WMAP7
Cosmology.WMAP5
```

## Retrieving Parameters
The parameters that define these types are accessed via the following unexported, internal methods.

```@docs
Cosmology.h
Cosmology.partype
Cosmology.m_nu
Cosmology.Neff
Cosmology.n_nu
Cosmology.w0
Cosmology.wa
```

Other parameters that have more complicated associated methods (e.g., parameters with redshift evolution like [`T_cmb`](@ref)) are given in the [Methods](@ref methods) section.