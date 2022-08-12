# [Methods](@id methods)
## Public Methods
The following methods are part of our publicly exported API.

### Temperatures

```@docs
T_cmb
T_nu
```

### Hubble Factors

```@docs
E
H
hubble_time
hubble_dist
```

### Lengths and Volumes

```@docs
comoving_radial_dist
comoving_transverse_dist
angular_diameter_dist
luminosity_dist
distmod
comoving_volume
comoving_volume_element
sound_horizon
```

### Times

```@docs
age
lookback_time
matter_radiation_equality
```

### Densities

```@docs
ρ_c
ρ_m
ρ_b
ρ_dm
ρ_Λ
ρ_γ
ρ_ν
ρ_r
Ω_m
Ω_b
Ω_dm
Ω_k
Ω_γ
Ω_ν
Ω_r
Ω_Λ
```

### Equation Solving

```@docs
z_at_value
```

### Miscellaneous

```@docs
scale_factor
∇scale_factor
```

# Private Methods
The following methods are used internally but not exported. They are not guaranteed to have stable APIs or be properly documented.

```@docs
Cosmology.nu_relative_density
Cosmology.a2E
Cosmology.hubble_dist0
Cosmology.hubble_time0
Cosmology.w
Cosmology.de_density_scale
```