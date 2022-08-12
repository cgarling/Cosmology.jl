# [Defined Types](@id types)
## Abstract Types
The highest level type defined in this package is the abstract type [`Cosmology.AbstractCosmology`](@ref Cosmology.AbstractCosmology), which all other cosmology objects are subtyped from. Below these are [`Cosmology.AbstractFlatCosmology`](@ref Cosmology.AbstractFlatCosmology), [`Cosmology.AbstractOpenCosmology`](@ref Cosmology.AbstractOpenCosmology), and [`Cosmology.AbstractClosedCosmology`](@ref Cosmology.AbstractClosedCosmology).
```@docs
Cosmology.AbstractCosmology
Cosmology.AbstractFlatCosmology
Cosmology.AbstractOpenCosmology
Cosmology.AbstractClosedCosmology
```

## Concrete Types
The following concrete types, representing specific types of cosmologies, are currently implemented:
 - `FlatLCDM`, `OpenLCDM`, `ClosedLCDM`, `FlatWCDM`, `OpenWCDM`,
 
## Constructors

## Retrieving Parameters
The parameters that define these types are accessed via the following unexported, internal methods

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