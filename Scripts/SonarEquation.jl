#!/bin/julia
# SonarEquation.jl
# 2023-05-09
# curtis
# This is a one-line description of the file.

# ArgParse: A package for handling command line arguments
using ArgParse

# PACKAGES ---------------------------------------------------------------------

using Distributions
using Interpolations
using Plots
using StatsPlots
using ProgressBars
using DataFrames
using Pipe

# STRUCTS ----------------------------------------------------------------------

"""
    Sonar

Sonar equation descriptors

Since multiple types of sonar equations require different treatments, the
abstract typing system provides the potential of common methods for all sonar
types.

...
# Subtypes
- [`sonar_noise`](@sonar_nosie): Sonar equations for noise environments
- [`sonar_passive`](@sonar_passive): Sonar equations for passive listening
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

See also [`Sonar`](@Sonar)
"""
struct sonar_noise <: Sonar
    sl :: Real
    tl :: Real
    ts :: Real
    nl :: UnivariateDistribution
    di :: Real
    dt :: Real
end

"""
    svp

Sound velocity profile object

Tracks the sound velocity profile as a function of depth

...
# Fields
- `depth :: Vector{Real}`: Depth of velocity
- `velocity :: Vector{Real}`: Velocity at some depth
...
"""
struct svp
    depth :: Vector{Real}
    velocity :: Vector{Real}
    function svp(depth, velocity)
        if length(depth) != length(velocity)
            error("depth and velocity of unequal length")
        end
        new(sort(depth), velocity[sortperm(depth)])
    end
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

See also [`Sonar`](@Sonar)
"""
struct sonar_passive <: Sonar
    sl :: Real
    tl :: Real
    nl :: UnivariateDistribution
    di :: Real
    dt :: Real
end

# FUNCTIONS --------------------------------------------------------------------

@doc raw"""
    freq_to_wavelength(velocity :: Real, freq :: Real) :: Real

Convert frequency to wavelength

Returns the result of $λ=v/f$, where $λ$ is the wavelength of a sound wave, $v$
is the velocity of sound in water, and $f$ is the frequency of the sound wave.

...
# Arguments
- `velocity :: Real`: Sound velocity
- `freq :: Real`: Wave frequency
...

See also [`wavelength_to_freq`](@wavelength_to_freq)

# Examples
```jldoctest
julia> freq_to_wavelength(4930.0, 100.0)
49.3
```
"""
function freq_to_wavelength(velocity :: Real, freq :: Real) :: Real
    velocity / freq
end

@doc raw"""
    wavelength_to_freq(velocity :: Real, wavelength :: Real) :: Real

Convert wavelength to frequency

Returns the result of $f=v/λ$, where $λ$ is the wavelength of a sound wave, $v$
is the velocity of sound in water, and $f$ is the frequency of the sound wave.

...
# Arguments
- `velocity :: Real`: Sound velocity
- `wavelength :: Real`: Sound wavelength
...

See also [`freq_to_wavelength`](@freq_to_wavelength)

# Examples
```jldoctest
julia> wavelength_to_freq(4930.0, 49.3)
100.0
```
"""
function wavelength_to_freq(velocity :: Real,
                            wavelength :: Real) :: Real
    freq_to_wavelength(velocity, wavelength)
end

@doc raw"""
    piston_di(diameter :: Real, wavelength :: Real) :: Real

Directivity index of a piston sonar

This is a theoretically derived directivity index of a piston sonar, which has
been found to be $\log_{10}\left(\left(\frac{π D}{λ}\right)^2\right), with $D$
being the diameter of the sonar and $λ$ the wavelength of the sound to be
detected.

...
# Arguments
- `diameter :: Real`: Diameter of the sonar
- `wavelength :: Real`: Wavelength of the sound to be detected
...

See also [`line_di`](@line_di)

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
    line_di(elements :: Int64, spacing :: Real,
            wavelength :: Real) :: Real

Directivity index of a line sonar

For a line transducer directivity index, theoretical calculations suggest the
directivity index is $\log_{10}\left(\frac{n}{1 + \frac{2}{n}
\sum_{ρ=1}^{n-1}\frac{(n - ρ)\sin(2ρπd/λ)}{2ρπd/λ}}\right)$, where $n$ is the
number of elements in the array, $d$ the spacing of the elements, and $λ$ the
wavelength of the sound wave to be detected.

...
# Arguments
- `elements :: Int64`: Number of elements in the sonar array
- `spacing :: Real`: The spacing of the elements in the array
- `wavelength :: Real`: The wavelength of the sound to be detected
...

See also [`piston_di`](@piston_di)

# Examples
```jldoctest
julia> line_di(20, 10.0, 49.3)
[...]
```
"""
function line_di(elements   :: Int64,
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

See also [`spherical_tl`](@spherical_tl)

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

See also [`attenuation_coef_thorp`](@attenuation_coef_thorp)

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

Compute $\text{SL} - 2\text{TL} + \text{TS} + \text{DI}$

If this quantity exceeds $\text{DT} - \text{NL}$, a detection occured.

...
# Arguments
- `se :: sonar_noise`: Sonar equation object
...

See also [`sonar_noise`](@sonar_noise), [`Base.rand`](@Base.rand)

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

Compute $\text{SL} - \text{TL} + \text{DI}$

If this quantity exceeds $\text{DT} - \text{NL}$, a detection occured.

...
# Arguments
- `se :: sonar_passive`: Sonar equation object
...

See also [`sonar_noise`](@sonar_noise), [`Base.rand`](@Base.rand)

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

See also [`Sonar`](@Sonar)

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

See also [`Sonar`](@Sonar)

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
    depth_slice_idx(svp :: svp, depth :: Real) :: 

Get index of depth slices

The returned integer corresponds to the correct depth slice of the [`svp`](@svp)
object. The input depth can exceed the lowest depth in `svp`, but not be too
shallow.

...
# Arguments
- `svp :: svp`: Sound velocity profile
- `depth :: Real`: Input depth
...

See also [`svp`](@svp)

# Examples
```jldoctest
julia> svp_mat = [# depth    velocity
                        0      1540.4
                       10      1540.5
                       20      1540.7
                       30      1534.4
                 ].* 3.28084
[...]
julia> test_svp = svp(svp_mat[:, 1], svp_mat[:, 2])
[...]
julia> depth_slice_idx(svp_mat, 25)
[...]
```
"""
function depth_slice_idx(svp :: svp, depth :: Real) :: UInt64
    if depth < svp.depth[1]
        throw(DomainError(depth, "Depth lower bound violated"))
    end

    sum(svp.depth.< depth)
end

@doc raw"""
    raytrace(svp :: svp, depth :: Real, angle :: Real, max_position :: Real;
    max_iter :: UInt128 = typemax(UInt128)) :: Vector{Tuple{Real, Real}}

Conduct raytracing

This function iteratively conducts raytracing using the depth velocity slices
provided by `svp`. It appropriately applies Snell's law to determine the next
angle of the ray as it travels through the medium. The last entries of `svp` are
presumed to be the ocean bottom, an interface the ray will bounce off; the
bottom and ocean surface (at a depth of 0) are presumed to be perfectly
reflective. (If a depth of 0 is not in `svp`, the raytracing algorithm will
terminate at the most shallow depth if the ray would extend beyond it.)

This function returns a vector of tuples, each generated by
[`raytrace_step`](@raytrace_step).

...
# Arguments
- `svp :: svp`: [`svp`](@svp) `struct` with the sound velocity profile data
- `depth :: Real`: The depth of the ray emitter
- `angle :: Real`: Initial angle of the ray
- `max_position :: Real`: The largest position before the raytrace terminates;
                          note that it is possible for a ray to extend beyond
                          the maximum position
- `max_iter :: UInt128`: Maximum number of iterations
...

See also [`raytrace_step`](@raytrace_step), [`svp`](@svp), [`snell`](@snell)

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
                      150      1517.3
                      200      1516.6
                      250      1516.5
                      300      1516.2
                      400      1516.4
                      500      1517.2
                      600      1518.2
                      700      1519.5
                      800      1521.0
                      900      1522.6
                     1000      1524.1
                     1100      1525.7
                     1200      1527.3
                     1300      1529.0
                     1400      1530.7
                     1500      1532.4
                     1750      1536.7
                     2000      1541.0
                 ].* 3.28084
[...]
julia> test_svp = svp(svp_mat[:, 1], svp_mat[:, 2])
[...]
julia> raytrace(test_svp, 500, 1 * π/180, 1000, max_iter = 1000)
[...]
```
"""
function raytrace(svp :: svp, depth :: Real, angle :: Real,
                  max_position :: Real;
                  max_iter :: UInt128 = typemax(UInt128)) :: Vector{Tuple{Real,
                                                                          Real}}
    if angle > π/2 || angle ≤ -π/2
        throw(DomainError(angle, "Must use radian angle ∈ (-π/2, π/2]"))
    elseif angle == 0
        return [(0, depth), (max_position, depth)]
    elseif angle == π/2
        return [(0, depth), (0, 0)]
    end

    trace = Tuple{Real, Real}[(0, depth)]
    depth_idx = depth_slice_idx(svp, depth)
    for iter in 1:max_iter
        velocity = svp.velocity[depth_idx]
        push!(trace,
              angle > 0 ? raytrace_step(svp.depth[depth_idx], depth, angle,
                                        trace[end][1]) :
                          raytrace_step(depth, svp.depth[depth_idx + 1], angle,
                                        trace[end][1]))
        if trace[end][1] ≥ max_position
            break
        end
        depth = trace[end][2]
        if depth == 0
            depth_idx = 1
        elseif depth == svp.depth[end]
            depth_idx = length(svp.depth) - 1
        else
            depth_idx = depth_slice_idx(svp, depth) + (angle < 0 ? 1 : 0)
            if depth_idx == 0
                break
            end
        end
        if depth ∈ Set([0, svp.depth[end]])
            angle *= -1
        else
            angle = snell(angle, velocity, svp.velocity[depth_idx])
        end
    end

    return trace
end

@doc raw"""
    raytrace_step(depth1 :: Real, depth2 :: Real, angle :: Real,
                  position :: Real) :: Tuple{Real, Real}

Compute the next position of a sound ray

A raytracing step draws a straight line between two layers of ocean depth, given
the angle of the ray. If $d_1$ is the first depth, $d_2$, the second depth with
$d_1 < d_2$, the ray enters at angle $θ$ and at horizontal position $p$, we can
determine the ray's new position and which depth slice it will appear. If $θ >
0$, the ray is travelling from $d_2$ to $d_1$; the ray's starting position is
$(p, d_2)$, and its new position will be

$$\left(\cot(θ)(d_2 - d_1) + p, d_1\right).$$

If the ray has an angle $θ < 0$, it is then starting from position $(p, d_1)$,
and will exit at

$$\left(-\cot(θ)(d_2 - d_1) + p, d_2\right).$$

The tuples above are returned by this function.

Exceptional cases include $θ = 0$ and $\left|θ\right| > π/2$. For horizontal
rays, the ray will continue onward forever. For vertical rays, it will go
directly up to the surface. Rays cannot go directly downward; it is presumed
that the ray will interact with something first and never reach the surface
(basically it will bounce off the emitter).

...
# Arguments
- `depth1 :: Real`: Upper depth
- `depth2 :: Real`: Lower depth
- `angle :: Real`: Angle of sound ray
- `position :: Real`: Current position of the sound ray in terms of distance
                      from the source
...

See also [`raytrace`](@raytrace)

# Examples
```jldoctest
julia> raytrace_step(10, 20, π/4, 100)
(10, 110.0)
```
"""
function raytrace_step(depth1 :: Real, depth2 :: Real, angle :: Real,
                       position :: Real) :: Tuple{Real, Real}
    if angle > π/2 || angle ≤ -π/2
        throw(DomainError(angle, "Must use radian angle ∈ (-π/2, π/2]"))
    elseif angle == 0
        return (max_position, depth)
    elseif angle == π/2
        return (0, 0)
    end
    if depth1 > depth2
        throw(DomainError(depth1, "Must have depth1 ≤ depth2"))
    end

    angle > 0 ? (cot(angle) * (depth2 - depth1) + position, depth1) :
        (-cot(angle) * (depth2 - depth1) + position, depth2)
end

@doc raw"""
    snell(angle :: Real, velocity1 :: Real, velocity2 :: Real) :: Real

Uses Snell's law to compute ray angle in new layer based on velocity comparison

This essentially solves $cos(θ_1)/v_1 = cos(θ_2)/v_2$ for $θ_2$, where $v_1$ is
the sound velocity of the depth layer being left, $θ_1$ the angle at which the
sound ray is leaving, and $v_2$ is the velocity of sound in the new layer.

...
# Arguments
- `angle :: Real`: Ray entering angle; must be between $-π/2$ and $π/2$
- `velocity1 :: Real`: Sound velocity at exiting layer
- `velocity2 :: Real`: Sound velocity at entering layer
...

See also [`raytrace`](@raytrace)

# Examples
```jldoctest
julia> snell(0.2, 5050, 5052)
[...]
```
"""
function snell(angle :: Real, velocity1 :: Real, velocity2 :: Real) :: Real
    if angle ≥ π/2 || angle ≤ -π/2
        throw(DomainError(angle, "Must use radian angle ∈ (-π/2, π/2)"))
    end

    acos_arg = velocity2 / velocity1 * cos(angle)
    if abs(acos_arg) ≤ 1
        return sign(angle) * acos(acos_arg)
    end

    -sign(angle) * acos(1/acos_arg)
end

"""
    svp_to_interpolation(svp_obj :: svp, itp_type :: InterpolationType,
                              ext_type :: BoundaryCondition) :: Any

Turn a [`svp`](@svp) into an interpolation function

This function's utility comes from its ability to take limited sound velocity
profile data and convert it into a function that allows for a better selection
of slices for raytracing. That allows for feeding more fine grained data to a
raytracing algorithm.

...
# Arguments
- `svp_obj :: svp`: [`svp`](@svp) `struct` to be turned into an interpolation
- `itp_type`: Set how interpolation will be done
- `ext_type`: Set how extrapolation will be done
...

See also [`svp`](@svp)

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
                      150      1517.3
                      200      1516.6
                      250      1516.5
                      300      1516.2
                      400      1516.4
                      500      1517.2
                      600      1518.2
                      700      1519.5
                      800      1521.0
                      900      1522.6
                     1000      1524.1
                     1100      1525.7
                     1200      1527.3
                     1300      1529.0
                     1400      1530.7
                     1500      1532.4
                     1750      1536.7
                     2000      1541.0
                 ].* 3.28084
[...]
julia> test_svp = svp(svp_mat[:, 1], svp_mat[:, 2])
[...]
julia> test_svp_itp = svp_to_interpolation(test_svp)
[...]
julia> fine_depths = [(0.0:1.0:14000.0)...]
[...]
julia> fine_test_svp = svp(fine_depths, test_svp_itp(fine_depths))
[...]
```
"""
function svp_to_interpolation(svp_obj :: svp,
                              itp_type = Gridded(Linear()),
                              ext_type = Linear()) :: Any
    extrapolate(interpolate((svp_obj.depth,), svp_obj.velocity, itp_type),
                ext_type)
end

"""
    svp_refine(svp_obj :: svp, Δ :: Real, min_depth :: Real = 0,
               max_depth :: Real = svp_obj.depth[end];
               itp_type = Gridded(Linear()), ext_type = Linear()) :: svp

Increase the resolution of a [`svp`](@svp) via interpolation and extrapolation

The resulting [`svp`](@svp) will have a constant distance between depths and may
extend to otherwise unseen depths.

...
# Arguments
- `svp_obj :: svp`: [`svp`](@svp) object to interpolate
- `Δ :: Real`: Difference between depth levels in new [`svp`](@svp)
- `min_depth :: Real = 0`: Minimum depth of new [`svp`](@svp)
- `max_depth :: Real = svp_obj.depth[end]`: Maximum depth of new [`svp`](@svp)
- `itp_type = Gridded(Linear())`: Set how interpolation will be done
- `ext_type = Linear()`: Set how extrapolation will be done
...

See also [`svp`](@svp), [`svp_to_interpolation`](@svp_to_interpolation)

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
```
"""
function svp_refine(svp_obj :: svp;
                    Δ :: Real = 1,
                    min_depth :: Real = 0,
                    max_depth :: Real = svp_obj.depth[end],
                    itp_type = Gridded(Linear()),
                    ext_type = Linear()) :: svp
    svp_itp = svp_to_interpolation(svp_obj, itp_type, ext_type)
    fine_depths = [(min_depth:Δ:max_depth)...]
    svp(fine_depths, svp_itp(fine_depths))
end

"""
    raytrace_angle_df(svp_obj :: svp, depth :: Real, angles :: Vector{Real},
                      max_position :: Real;
                      max_iter :: UInt128 = typemax(UInt128)) :: DataFrame

Raytrace by angle

Generates multiple ray tracings depending on angle and generates a `DataFrame`
with rays depending on initial angle.

...
# Arguments
- `svp_obj`: A [`svp`](@svp) object
- `depth`: Depth of sound emitter
- `angles`: Angles to trace
- `max_position`: Maximum position of a ray
- `position_grid`: Positions at which ray position should be known; may require
                   interpolating the traced rays, yet another source of
                   inaccuracy, especially when bounce
- `itp_type`: Type of interpolator
- `max_iter = typemax(UInt128)`: Maximum iteration number for raytracing
...

See also [`raytrace`](raytrace)

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
julia> raytrace_angle_df(fine_test_svp, 50, [-2 * π/180, 2 * π/180], 20)
[...]
```
"""
function raytrace_angle_df(svp_obj, depth,
                           angles,
                           max_position;
                           position_grid = [(0.0:1.0:max_position)...],
                           itp_type = Gridded(Linear()),
                           max_iter = typemax(UInt128)) :: DataFrame
    multiray = raytrace.(Ref(svp_obj), depth, angles, max_position,
                         max_iter = max_iter)
    if isnothing(position_grid)
        vcat([DataFrame(angle = ang,
                        range = [p[1] for p in ray],
                        depth = [p[2] for p in ray])
              for (ray, ang) in zip(multiray, angles)]..., cols = :union)
    else
        vcat([DataFrame(angle = ang,
                        range = position_grid,
                        depth = (interpolate((Interpolations.deduplicate_knots!([convert(Float64, p[1]) for p in ray]),),
                                             [p[2] for p in ray],
                                             itp_type))(position_grid))
              for (ray, ang) in zip(multiray, angles)]..., cols = :union)
    end
end

"""
    ray_position_above_df(ray_df :: DataFrame, range :: Real, depth ::) :: DataFrame

Determine whether rays in a raytrace `DataFrame` are above or below a point

The significance of the calculation is determining which rays should be involved
in calculating transmission loss.

...
# Arguments
- `ray_df :: DataFrame`: A `DataFrame` with columns `angle`, `range`, and
             `depth` containing traced ray positions
- `range :: Real`: The range at which to compute transmission loss
- `depth :: Real`: The depth of the point at which to compute transmission loss
...

See also [`raytrace_angle_df`](@raytrace_angle_df)

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
julia> rad = raytrace_angle_df(fine_test_svp, 50, [-2 * π/180, 2 * π/180], 20)
[...]
julia> ray_position_above_df(rad, 20, 50)
[...]
```
"""
function ray_position_above_df(ray_df :: DataFrame, range :: Real, depth :: Real) :: DataFrame
    if range < min((ray_df.range)...)
        error("Input range less than minimum range in ray_df")
    end
    functional_range = max((ray_df.range[ray_df.range .<= range])...)
    ray_df_range = ray_df[ray_df.range .== functional_range, :]
    DataFrame(angle = ray_df_range.angle, above = ray_df_range.depth .> depth,
              depth = ray_df_range.depth) 
end

"""
    get_containing_rays_df(ray_df :: DataFrame) :: DataFrame

Get rays above or below a given position depending on angle and the depth of the
rays

This is essentially the next step after using
[`ray_position_above_df`](@ray_position_above_df).

...
# Arguments
- `ray_df :: DataFrame`: Contains columns `angle` (`Real`), `above` (`Bool`),
                         and `depth` (`Real`), describing ray depth given angle
                         and whether a ray is above or below a certain point
...

See also [`ray_position_above_df`](@ray_position_above_df)

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
julia> rad = raytrace_angle_df(fine_test_svp, 50, [-2 * π/180, 2 * π/180], 20)
[...]
julia> rad_ab = ray_position_above_df(rad, 20, 50)
[...]
julia> get_containing_rays_df(rad_ab)
[...]
```
"""
function get_containing_rays_df(ray_df :: DataFrame) :: DataFrame
    ray_df = deepcopy(ray_df)
    sort!(ray_df, [:angle])
    above_below_vec = ray_df.above[1:end-1] .⊻ ray_df.above[2:end]
    hcat(select(ray_df[[above_below_vec; false], :], :angle => :lower_angle,
                :depth => :depth_lower_angle),
         select(ray_df[[false; above_below_vec], :], :angle => :upper_angle,
                :depth => :depth_upper_angle))
end

@doc raw"""
    raytrace_tl(lower_angle :: Real, upper_angle :: Real, depth1 :: Real,
                depth2 :: Real,  range :: Real) :: Real

Computes transmission loss from a raytracing diagram

The formula for target loss given a raytracing diagram is $\text{TL} = 10 \log
\left(rΔh/Δθ\right)$, where $r$ is range, $Δh$ is the difference in depth
between two rays separated by $Δθ$ radians in angle and adjacent at the origin
of the noise.

...
# Arguments
- `lower_angle :: Real`: Lower angle, in radians
- `upper_angle :: Real`: Upper angle, in radians
- `depth1 :: Real`: Depth in yards
- `depth2 :: Real`: Depth in yards
- `range :: Real`: Range to target, in yards
...

See also [`get_containing_rays_df`](@get_containing_rays_df),
[`raytrace_angle_df`](@raytrace_angle_df)

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

# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION -------------------------

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--opt1"
            help = "an option with an argument"
        "--opt2", "-o"
            help = "another option with an argument"
            arg_type = Int
            default = 0
        "--flag1"
            help = "an option without an argument, i.e. a flag"
            action = :store_true
        "arg1"
            help = "a positional argument"
            required = true
    end

    return parse_args(s)
end

# EXECUTABLE SCRIPT MAIN FUNCTIONALITY -----------------------------------------

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("    $arg => $val")
    end

    svp_mat = [# depth    velocity
                     0      1540.4
                    10      1540.5
                    20      1540.7
                    30      1534.4
                    50      1523.3
                    75      1519.6
                   100      1518.5
                   125      1517.9
                   150      1517.3
                   200      1516.6
                   250      1516.5
                   300      1516.2
                   400      1516.4
                   500      1517.2
                   600      1518.2
                   700      1519.5
                   800      1521.0
                   900      1522.6
                  1000      1524.1
                  1100      1525.7
                  1200      1527.3
                  1300      1529.0
                  1400      1530.7
                  1500      1532.4
                  1750      1536.7
                  2000      1541.0
              ].* 3.28084
    test_svp = svp(svp_mat[:, 1], svp_mat[:, 2])
    fine_test_svp = svp_refine(test_svp, max_depth = 14000)

    angles = [(-5.0:1.0:5.0)...].*π./180
    tgtdepth = 100
    ray_combined_df = raytrace_angle_df(fine_test_svp, 100.0, angles, 30 * 6000)

    @df ray_combined_df plot(:range, :depth, group = :angle, yflip = true)

    detector_depth = 5500
    detector_range = 60000.0
    containing_ray_pos = ray_position_above_df(ray_combined_df, detector_range,
                                               detector_depth) |>
        get_containing_rays_df
    
    mean([raytrace_tl(
              r.lower_angle, r.upper_angle, r.depth_lower_angle/3,
              r.depth_upper_angle/3, detector_range/3
             ) for r in eachrow(containing_ray_pos)])
end

if !isinteractive()
    main()
end

