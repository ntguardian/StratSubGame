#!/bin/julia
# SonarEquationTools.jl
# 2023-07-24
# curtis
# Sonar equation tools: structs and functions that work with the sonar equations

# STRUCTS ----------------------------------------------------------------------

"""
    Sonar

Sonar equation descriptors

Since multiple types of sonar equations require different treatments, the
abstract typing system provides the potential of common methods for all sonar
types.

...
# Subtypes
- [`sonar_noise`](#SimpleSonar.sonar_nosie): Sonar equations for noise
                                             environments
- [`sonar_passive`](#SimpleSonar.sonar_passive): Sonar equations for passive
                                                 listening
...
"""
abstract type Sonar end

@doc raw"""
    sonar_noise

Sonar equation parameters for noise conditions

The sonar equation for noise conditions is
```math
\text{SL} - 2\text{TL} + \text{TS} = \text{NL} - \text{DI} + \text{DT}
```
with SL representing the source level, TL the transmission loss, TS the target
strength, NL the noise level, and DT the detection threshold (all in decibels).

...
# Fields
- `sl :: Real`: Source level
- `tl :: Real`: Transmission loss
- `ts :: Real`: Transmission strength
- `nl :: UnivariateDistribution`: Noise level
- `di :: Real`: Directivity index
- `dt :: Real`: Detection threshold
...

See also [`Sonar`](#SimpleSonar.Sonar)
"""
struct sonar_noise <: Sonar
    sl :: Real
    tl :: Real
    ts :: Real
    nl :: UnivariateDistribution
    di :: Real
    dt :: Real
end

@doc raw"""
    sonar_passive

Sonar equation parameters for passive sonar conditions

The sonar equation for passive sonar detection are
```math
\text{SL} - \text{TL} = \text{NL} - \text{DI} + \text{DT}
```
with SL representing the source level, TL the transmission loss, NL the noise
level, and DT the detection threshold (all in decibels).

...
# Fields
- `sl :: Real`: Source level
- `tl :: Real`: Transmission loss
- `nl :: UnivariateDistribution`: Noise level
- `di :: Real`: Directivity index
- `dt :: Real`: Detection threshold
...

See also [`Sonar`](#SimpleSonar.Sonar)
"""
struct sonar_passive <: Sonar
    sl :: Real
    tl :: Real
    nl :: UnivariateDistribution
    di :: Real
    dt :: Real
end

@doc raw"""
    sonar_reverb

Sonar equation parameters for passive sonar conditions

The sonar equation for active sonar with reverberation is
```math
\text{SL} - 2\text{TL} + \text{TS} = \text{RL} + \text{DT}
```
with SL representing the source level, TL the transmission loss, RL the
reverberation level, TS the target strength, and DT the detection threshold (all
in decibels).

...
# Fields
- `sl :: Real`: Source level
- `tl :: Real`: Transmission loss
- `ts :: Real`: Target strength
- `rl :: UnivariateDistribution`: Reverberation level
- `dt :: Real`: Detection threshold
...

See also [`Sonar`](#SimpleSonar.Sonar)
"""
struct sonar_reverb <: Sonar
    sl :: Real
    tl :: Real
    ts :: Real
    rl :: UnivariateDistribution
    dt :: Real
end

# FUNCTIONS --------------------------------------------------------------------

@doc raw"""
    piston_di(diameter :: Real, wavelength :: Real) :: Real

Directivity index of a piston sonar

This is a theoretically derived directivity index of a piston sonar, which has
been found to be ``\log_{10}\left(\left(\frac{\pi D}{\lambda}\right)^2\right)``, with ``D``
being the diameter of the sonar and ``\lambda`` the wavelength of the sound to be
detected.

...
# Arguments
- `diameter :: Real`: Diameter of the sonar
- `wavelength :: Real`: Wavelength of the sound to be detected
...

See also [`line_di`](#SimpleSonar.line_di)

# Examples
```jldoctest
julia> piston_di(200.0, 49.3)
[...]
```
"""
function piston_di(diameter :: Real, wavelength :: Real) :: Real
    20 * (log10(π) + log10(diameter) - log10(wavelength))
end

@doc raw"""
    line_di(elements :: Unsigned, spacing :: Real,
            wavelength :: Real) :: Real

Directivity index of a line sonar

For a line transducer directivity index, theoretical calculations suggest the
directivity index is
```math
\log_{10}\left(\frac{n}{1 + \frac{2}{n}
\sum_{\rho=1}^{n-1}\frac{(n - \rho)\sin(2\rho\pid/\lambda)}{2\rho\pid/\lambda}}\right)
``` 
where ``n`` is the
number of elements in the array, ``d`` the spacing of the elements, and ``\lambda`` the
wavelength of the sound wave to be detected.

...
# Arguments
- `elements :: Unsigned`: Number of elements in the sonar array
- `spacing :: Real`: The spacing of the elements in the array
- `wavelength :: Real`: The wavelength of the sound to be detected
...

See also [`piston_di`](#SimpleSonar.piston_di)

# Examples
```jldoctest
julia> line_di(20, 10.0, 49.3)
[...]
```
"""
function line_di(elements   :: Unsigned,
                 spacing    :: Real,
                 wavelength :: Real) :: Real
    log10(elements / (1 + 2 / elements *
                      sum([(elements - idx) *
                           sin(2 * idx * π * spacing / wavelength) /
                           (2 * idx * π * spacing / wavelength)
                           for idx in 1:(elements - 1)])))
end

"""
    attenuation_coef_thorp(freq) :: Real

Thorp's attenuation coefficient

Thorp's computation is for 39°F/4°C at 3000 ft deep

...
# Arguments
- `freq :: Real`: Frequency (kHz)
...

See also [`spherical_tl`](#SimpleSonar.spherical_tl)

# Examples
```jldoctest
julia> attenuation_coef_thorp(0.1)
[...]
```
"""
function attenuation_coef_thorp(freq :: Real) :: Real
    0.1 * freq^2 / (1 + freq^2) + 40 * freq^2 / (4100 + freq^2) + 2.75 *
        10^(-4) * freq^2 + 0.003
end

"""
    spherical_tl(attenuation :: Real, range :: Real) :: Real

Spherical spreading transmission loss

Transmission loss, as described by spherical spreading

...
# Arguments
- `attenuation :: Real`: Attenuation, in decibels per kiloyard
- `range :: Real`: Range, in yards
...

See also [`attenuation_coef_thorp`](#SimpleSonar.attenuation_coef_thorp)

# Examples
```jldoctest
julia> spherical_tl(attenuation_coef_thorp(0.1), 2 * 300.0)
[...]
```
"""
function spherical_tl(attenuation :: Real, range :: Real) :: Real
    20 * log10(range) + attenuation * range / 1000.0
end

@doc raw"""
    sonar_threshold(se :: sonar_noise) :: Real

Compute ``\text{SL} - 2\text{TL} + \text{TS} + \text{DI}``

If this quantity exceeds ``\text{DT} - \text{NL}``, a detection occured.

...
# Arguments
- `se :: sonar_noise`: Sonar equation object
...

See also [`sonar_noise`](#SimpleSonar.sonar_noise), [`Base.rand`](#Base.rand)

# Examples
```jldoctest
julia> se_ex = sonar_noise(1.0, 1.0, 1.0, Normal(1.0, 1.0), 15.0, 15.0)
julia> sonar_threshold(se_ex)
[...]
```
"""
function sonar_threshold(se :: sonar_noise) :: Real
    se.sl - 2 * se.tl + se.ts + se.di
end

@doc raw"""
    sonar_threshold(se :: sonar_passive) :: Real

Compute ``\text{SL} - \text{TL} + \text{DI}``

If this quantity exceeds ``\text{DT} - \text{NL}``, a detection occured.

...
# Arguments
- `se :: sonar_passive`: Sonar equation object
...

See also [`sonar_noise`](#SimpleSonar.sonar_noise), [`Base.rand`](#SimpleSonar.Base.rand)

# Examples
```jldoctest
julia> se_ex = sonar_passive(1.0, 1.0, Normal(1.0, 1.0), 15.0, 15.0)
julia> sonar_threshold(se_ex)
[...]
```
"""
function sonar_threshold(se :: sonar_passive) :: Real
    se.sl -  se.tl + se.di
end

"""
    Base.rand(rng :: AbstractRNG, se :: Sonar) :: Bool

Random detection for `Sonar` objects

Randomly determine a detection event of a `Sonar` object by determining if the
detection threshold was exceeded. Note that only the noise level is treated as
random in this very simple model, and there is no possibility of signal
fluctuation.

...
# Arguments
- `se :: Sonar`: Sonar object to randomize
...

See also [`Sonar`](#SimpleSonar.Sonar)

# Examples
```julia-repl
julia> se_ex = sonar_noise(1.0, 1.0, 1.0, Normal(1.0, 1.0), 15.0, 15.0)
julia> Base.rand(se_ex)
false
```
"""
function Base.rand(se :: Sonar) :: Bool
    sonar_threshold(se) > se.dt + rand(se.nl)
end

@doc raw"""
    detection_prob(se :: Sonar) :: Float64

Sonar equation to detection probability

Determines the probability that the sonar signal will exceed the detection
threshold plus noise, given random noise. Note that this is not necessarily
figuring out the set detection probability as determined by the ROC curve used
for determining the detection threshold.

...
# Arguments
- `se :: Sonar`: Sonar object to randomize
...

See also [`Sonar`](#SimpleSonar.Sonar)

# Examples
```jldoctest
julia> se_ex = sonar_noise(1.0, 1.0, 1.0, Normal(1.0, 1.0), 15.0, 15.0)
julia> detection_prob(se_ex)
[...]
```
"""
function detection_prob(se :: Sonar) :: Float64
    cdf(se.nl, sonar_threshold(se) - se.dt)
end

@doc raw"""
    raytrace_tl(lower_angle :: Real, upper_angle :: Real, depth1 :: Real,
                depth2 :: Real,  range :: Real) :: Real

Computes transmission loss from a raytracing diagram

The formula for target loss given a raytracing diagram is $\text{TL} = 10 \log
\left(r\Deltah/\Delta\theta\right)$, where $r$ is range, $\Deltah$ is the difference in depth
between two rays separated by $\Delta\theta$ radians in angle and adjacent at the origin
of the noise.

...
# Arguments
- `lower_angle :: Real`: Lower angle, in radians
- `upper_angle :: Real`: Upper angle, in radians
- `depth1 :: Real`: Depth in yards
- `depth2 :: Real`: Depth in yards
- `range :: Real`: Range to target, in yards
...

See also [`get_containing_rays_df`](#SimpleSonar.get_containing_rays_df),
[`raytrace_angle_df`](#SimpleSonar.raytrace_angle_df)

# Examples
```jldoctest
julia> raytrace_tl(0.01, 0.02, 100, 102, 20000)
[...]
```
"""
function raytrace_tl(lower_angle :: Real, upper_angle :: Real, depth1 :: Real,
                     depth2 :: Real,  range :: Real) :: Real
    if range ≤ 0
        error("Only positive range allowed")
    end
    if lower_angle ≥ upper_angle
        error("Must have lower_angle < upper_angle")
    end
    
    10 * log10(range * abs(depth1 - depth2) / (upper_angle - lower_angle))
end

"""
    ray_df_to_tl(ray_df :: DataFrame, range :: Real, depth :: Real) :: Real

Obtain transmission loss estimate from `DataFrame` of traced sound rays

When multiple rays pass above or below a point, which ones to use for computing
transmission loss is not clear. I opt here to use the minimum transmission loss.

...
# Arguments
- `ray_df :: DataFrame`: `DataFrame` of rays traced with a function such as
                         [`raytrace_angle_df`](#SimpleSonar.raytrace_angle_df)
- `range :: Real`: The range of the sensor
- `depth :: Real`: The depth of the sensor
...

See also [`ray_position_above_df`](#SimpleSonar.ray_position_above_df),
[`get_containing_rays_df`](#SimpleSonar.get_containing_rays_df),
[`raytrace_angle_df`](#SimpleSonar.raytrace_angle_df)

# Examples
```jldoctest
julia> svp_mat = [# depth    velocity
                        0      1540.4
                       10      1540.5
                       20      1540.7
                       30      1534.4
                       50      1523.3
                       75      1519.6
                      100      1518.5
                      125      1517.9
                 ].* 3.28084
[...]
julia> test_svp = svp(svp_mat[:, 1], svp_mat[:, 2])
[...]
julia> fine_test_svp = svp_refine(test_svp, max_depth = 200)
[...]
julia> rad = raytrace_angle_df(fine_test_svp, 50, [-2 * π/180, 2 * π/180], 20)[:ray]
[...]
julia> ray_df_to_tl(rad, 20, 50)
[...]
```
"""
function ray_df_to_tl(ray_df :: DataFrame, range :: Real, depth :: Real) :: Real
    containing_ray_pos = ray_position_above_df(ray_df, range,
                                               depth) |>
        get_containing_rays_df
    
    if nrow(containing_ray_pos) == 0
        return Inf
    else
        return  min([raytrace_tl(
                    r.lower_angle, r.upper_angle, r.depth_lower_angle/3,
                    r.depth_upper_angle/3, range/3
                   ) for r in eachrow(containing_ray_pos)]...)
    end
end

