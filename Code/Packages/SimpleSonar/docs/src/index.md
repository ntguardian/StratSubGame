# SimpleSonar.jl Documentation

```@contents
```

## Mathematics Functions

```@docs
freq_to_wavelength
wavelength_to_freq
```

## Raytracing Functions

```@docs
raytrace
raytrace_step
snell
raytrace_angle_df
ray_position_above_df
get_containing_rays_df
```

## Sonar Equation Tools

```@docs
Sonar
sonar_noise
sonar_passive
piston_di
line_di
attenuation_coef_thorp
spherical_tl
sonar_threshold
sonar_threshold
Base.rand(::Sonar)
detection_prob
raytrace_tl
ray_df_to_tl
```

## Sound Velocity Profile Tools

```@docs
svp
depth_slice_idx
svp_to_interpolation
svp_refine
DataFrame(::svp)
```

