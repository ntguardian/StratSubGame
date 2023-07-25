#!/bin/julia
# SonarEquation.jl
# 2023-05-09
# curtis
# This is a one-line description of the file.

# ArgParse: A package for handling command line arguments
using ArgParse

# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION -------------------------

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "ddepth"
            arg_type = Float64
            help = "Depth of detector (ft)"
            required = true
        "drange"
            arg_type = Float64
            help = "Range of the detector from the target (ft)"
            required = true
        "--emtdepth", "-E"
            arg_type =Float64
            help = "Depth of sound emitter (ft)"
            default = nothing
        "--tgtrange", "-T"
            arg_type = Float64
            help = "Maximum distance for ray tracing (ft)"
            default = nothing
        "--maxdepth", "-M"
            arg_type = Float64
            default = 36161.0
            help = "The maximum depth (ft) of the ocean, the maximum presuming to be the ocean bottom"
        "--minangle", "-a"
            help = "Minimum angle (degrees) to use for raytracing"
            arg_type = Float64
            default = nothing
        "--maxangle", "-A"
            help = "Maximum angle (degrees) to use for raytracing"
            arg_type = Float64
            default = nothing
        "--stepangle", "-s"
            help = "Difference between angles (degrees) to use for raytracing"
            arg_type = Float64
            default = nothing
        "--svpstep", "-P"
            help = "When dividing SVP, new (constant) distance between slices (ft)"
            arg_type = Float64
            default = nothing
        "--velocity", "-v"
            help = "Sound velocity to use if SVP not in use (ft/sec)"
            arg_type = Float64
            default = nothing
        "--freq", "-f"
            help = "Frequency of sound (Hz)"
            arg_type = Float64
            default = nothing
        "--pdiameter", "-i"
            help = "Diameter of piston hydrophone (ft) (supply if wanting calculations of DI for this type of hydrophone)"
            arg_type = Float64
            default = nothing
        "--lelements", "-e"
            help = "Number of elements in a line hydrophone (supply if wanting calculation of DI for this type of hydrophone)"
            arg_type = UInt128
            default = nothing
        "--lspacing", "-p"
            help = "Spacing of elements in a line hydrophone (ft) (supply if wanting calculation of DI for this type of hydrophone)"
            arg_type = Float64
            default = nothing
        "--di", "-I"
            help = "Directivity index (dB); no calculation of DI will be done if supplied"
            arg_type = Float64
            default = nothing
        "--dionly"
            help = "Compute and report DI only"
            action = :store_true
        "--tl", "-t"
            help = "Transmission loss (dB); no calculation of TL will be done if supplied"
            arg_type = Float64
            default = nothing
        "--tlonly"
            help = "Compute and report TL only"
            action = :store_true
        "--raycsv", "-r"
            help = "Location to save ray tracing results (for use later)"
            arg_type = String
            default = nothing
        "--rayloadcsv", "-R"
            help = "Location of CSV file with precomputed ray tracing paths (will skip ray tracing computation)"
            arg_type = String
            default = nothing
        "--bouncecsv", "-b"
            help = "Location to save ray tracing bounce information (for use later)"
            arg_type = String
            default = nothing
        "--svpcsv", "-S"
            help = "Location of CSV file containing SVP information (depth (ft), velocity (ft/sec))"
            arg_type = String
            default = nothing
        "--sl", "-l"
            help = "Source level (dB)"
            arg_type = Float64
            default = nothing
        "--nlmean", "-N"
            help = "Mean noise level (dB)"
            arg_type = Float64
            default = nothing
        "--nlsd", "-n"
            help = "Standard deviation of noise level (dB)"
            arg_type = Float64
            default = nothing
        "--dt", "-d"
            help = "Detection threshold (dB)"
            arg_type = Float64
            default = nothing
    end

    return Dict([(Symbol(key), val) for (key, val) in parse_args(s)])
end

if !isinteractive()
    parsed_args = parse_commandline()
end

# PACKAGES ---------------------------------------------------------------------

using Distributions
using Interpolations
using Plots
using StatsPlots
using DataFrames
using CSV

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
        if length(depth) ≠ length(velocity)
            error("depth and velocity of unequal length")
        end
        new(sort(depth), velocity[sortperm(depth)])
    end
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
    for i in eachindex(x)
        @inbounds x[i] = d.radius * func_vec[1 + mod(i, 2)](d.θ_min +
            (d.θ_max - d.θ_min) * rand(rng)) + d.center[1 + mod(i, 2)]
    end
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
    line_di(elements :: Unsigned, spacing :: Real,
            wavelength :: Real) :: Real

Directivity index of a line sonar

For a line transducer directivity index, theoretical calculations suggest the
directivity index is $\log_{10}\left(\frac{n}{1 + \frac{2}{n}
\sum_{ρ=1}^{n-1}\frac{(n - ρ)\sin(2ρπd/λ)}{2ρπd/λ}}\right)$, where $n$ is the
number of elements in the array, $d$ the spacing of the elements, and $λ$ the
wavelength of the sound wave to be detected.

...
# Arguments
- `elements :: Unsigned`: Number of elements in the sonar array
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

This function returns a tuple with an element `trace` being a vector of tuples
with the position and depth of the ray, and an element `bounce`, a vector of
tuples with elements `position`, `depth`, and `angle`, giving the grazing angle
and location of either surface or bottom bounces of the ray.

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
julia> raytrace(test_svp, 500, 1 * π/180, 1000, max_iter = UInt128(1000))
[...]
```
"""
function raytrace(svp :: svp, depth :: Real, angle :: Real,
                  max_position :: Real;
                  max_iter :: UInt128 = typemax(UInt128))
    bounce = NamedTuple{(:position, :depth, :angle), Tuple{Real, Real, Real}}[]
    if angle > π/2 || angle ≤ -π/2
        throw(DomainError(angle, "Must use radian angle ∈ (-π/2, π/2]"))
    elseif angle == 0
        return (trace=[(0, depth), (max_position, depth)], bounce=bounce)
    elseif angle == π/2
        return (trace=[(0, depth), (0, 0)],
                bounce=[(position=0, depth=0, angle=angle)])
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
            push!(bounce, (position=trace[end][1], depth=trace[end][2],
                           angle=angle))
            angle *= -1
        else
            angle = snell(angle, velocity, svp.velocity[depth_idx])
        end
    end

    return (trace=trace, bounce=bounce)
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
                      max_iter :: UInt128 = typemax(UInt128)) :: NamedTuple{(:ray, :bounce), Tuple{DataFrame, DataFrame}}

Raytrace by angle

Generates multiple ray tracings depending on angle and generates a `DataFrame`
with rays depending on initial angle. Also collects bounce information.

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
                           max_iter = typemax(UInt128)) :: NamedTuple{(:ray, :bounce), Tuple{DataFrame, DataFrame}}
    multiray = raytrace.(Ref(svp_obj), depth, angles, max_position,
                         max_iter = max_iter)
    if isnothing(position_grid)
        ray_df = vcat([DataFrame(angle = ang,
                                 range = [p[1] for p in ray[:trace]],
                                 depth = [p[2] for p in ray[:trace]])
                       for (ray, ang) in zip(multiray, angles)]..., cols = :union)
    else
        ray_df = vcat([DataFrame(angle = ang,
                                 range = position_grid,
                                 depth = (interpolate((Interpolations.deduplicate_knots!([convert(Float64, p[1]) for p in ray[:trace]]),),
                                                      [p[2] for p in ray[:trace]],
                                                      itp_type))(position_grid))
                       for (ray, ang) in zip(multiray, angles)]..., cols = :union)
    end

    bounce_df = vcat([DataFrame(angle = ang,
                                position = [r[:position] for r in ray[:bounce]],
                                depth = [r[:depth] for r in ray[:bounce]],
                                graze = [r[:angle] for r in ray[:bounce]])
                      for (ray, ang) in zip(multiray, angles)]..., cols = :union)
    
    return (ray=ray_df, bounce=bounce_df)
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
julia> rad = raytrace_angle_df(fine_test_svp, 50, [-2 * π/180, 2 * π/180], 20)[:ray]
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
julia> rad = raytrace_angle_df(fine_test_svp, 50, [-2 * π/180, 2 * π/180], 20)[:ray]
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

"""
    DataFrame(svp_obj :: svp) :: DataFrame

Convert [`svp`](@svp) structs to `DataFrame`s.

...
# Arguments
- `svp_obj :: svp`: The SVP to convert
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
julia> DataFrame(test_svp)
[...]
```
"""
function DataFrame(svp_obj :: svp) :: DataFrame
    DataFrame(depth = svp_obj.depth, velocity = svp_obj.velocity)
end

"""
    ray_df_to_tl(ray_df :: DataFrame, range :: Real, depth :: Real) :: Real

Obtain transmission loss estimate from `DataFrame` of traced sound rays

When multiple rays pass above or below a point, which ones to use for computing
transmission loss is not clear. I opt here to use the minimum transmission loss.

...
# Arguments
- `ray_df :: DataFrame`: `DataFrame` of rays traced with a function such as
                         [`raytrace_angle_df`](@raytrace_angle_df)
- `range :: Real`: The range of the sensor
- `depth :: Real`: The depth of the sensor
...

See also [`ray_position_above_df`](@ray_position_above_df),
[`get_containing_rays_df`](@get_containing_rays_df),
[`raytrace_angle_df`](raytrace_angle_df)

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
                    r.depth_upper_angle/3, drange/3
                   ) for r in eachrow(containing_ray_pos)]...)
    end
end

# EXECUTABLE SCRIPT MAIN FUNCTIONALITY -----------------------------------------

"""
See [`parse_commandline`](@parse_commandline) for argument description
"""
function main(;
              ddepth :: Real,
              drange :: Real,
              emtdepth :: Union{Real, Nothing} = nothing,
              tgtrange :: Union{Nothing, Real} = nothing,
              maxdepth :: Union{Nothing, Real} = nothing,
              minangle :: Union{Nothing, Real} = nothing,
              maxangle :: Union{Nothing, Real} = nothing,
              stepangle :: Union{Nothing, Real} = nothing,
              svpstep :: Union{Nothing, Real} = nothing,
              velocity :: Union{Nothing, Real} = nothing,
              freq :: Union{Nothing, Real} = nothing,
              pdiameter :: Union{Nothing, Real} = nothing,
              lelements :: Union{Nothing, Unsigned} = nothing,
              lspacing :: Union{Nothing, Real} = nothing,
              di :: Union{Nothing, Real} = nothing,
              dionly :: Bool = false,
              tl :: Union{Nothing, Real} = nothing,
              tlonly :: Bool = false,
              raycsv :: Union{Nothing, String} = nothing,
              raydf :: Union{Nothing, DataFrame} = nothing,    # Not used by command line interface
              rayloadcsv :: Union{Nothing, String} = nothing,
              bouncecsv :: Union{Nothing, String} = nothing,
              svpcsv :: Union{Nothing, String} = nothing,
              svpmat :: Union{Nothing, Matrix{<:Real}} = nothing,    # Not used by command line interface
              sl :: Union{Nothing, Real} = nothing,
              nlmean :: Union{Nothing, Real} = nothing,
              nlsd :: Union{Nothing, Real} = nothing,
              dt :: Union{Nothing, Real} = nothing
             )
    svp_obj = nothing
    if !isnothing(svpcsv)
        svp_df = DataFrame(CSV.File(svpcsv))
        svp_obj = svp(svp_df.depth, svp_df.velocity)
        svp_itp = svp_to_interpolation(svp_obj, Gridded(Linear()), Linear())
        velocity = svp_itp(ddepth)
    elseif !isnothing(svpmat)
        svp_obj = svp(svpmat[:,1], svpmat[:,2])
        svp_itp = svp_to_interpolation(svp_obj, Gridded(Linear()), Linear())
        velocity = svp_itp(ddepth)
    end
    
    if isnothing(di)
        if (isnothing(velocity) || isnothing(freq))
            error("Must have velocity and freq if di not specified")
        end
        wavelength = freq_to_wavelength(velocity, freq)
        if !isnothing(pdiameter) && !isnothing(pelements)
            # Assuming piston hydrophone
            di = piston_di(pdiameter, wavelength)
        elseif !isnothing(lelements) && !isnothing(lspacing)
            # Assuming line sonar
            di = line_di(lelements, lspacing, wavelength)
        elseif dionly
            error("Not enough information to determine directivity index")
        end
    end
    if dionly
        println(di)
        return nothing
    end

    if isnothing(tl)
        # Logic here could be tightened; it seems that the current logic implied
        # here will result in duplicate code
        if !isnothing(raydf)
            ray_combined_df = raydf
            tl = ray_df_to_tl(ray_combined_df, drange, ddepth)
        elseif !isnothing(rayloadcsv)
            ray_combined_df = DataFrame(CSV.File(rayloadcsv))
            tl = ray_df_to_tl(ray_combined_df, drange, ddepth)
        elseif !isnothing(svp_obj)
            if (isnothing(maxdepth) || isnothing(minangle) ||
                isnothing(maxangle) || isnothing(stepangle) ||
                isnothing(tgtrange) || isnothing(emtdepth))
                error("Called in raytracing mode; must have all of maxdepth, minangle, stepangle, maxangle, tgtrange, emtdepth")
            end
            fine_svp_obj = svp_refine(svp_obj, max_depth = maxdepth,
                                      Δ = svpstep)

            # plot(fine_svp_obj.velocity, fine_svp_obj.depth, yflip = true,
            #      legend = false)

            if maxangle < minangle
                error("Must have minangle ≤ maxangle")
            end
            angles = [(minangle:stepangle:maxangle)...].*π./180

            ray_combined_result = raytrace_angle_df(fine_svp_obj, emtdepth, angles,
                                                    tgtrange)
            ray_combined_df = ray_combined_result[:ray]
            ray_bounce_bounce = ray_combined_result[:bounce]

            if !isnothing(raycsv)
                CSV.write(raycsv, ray_combined_df)
            end
            if !isnothing(bouncecsv)
                CSV.write(bouncecsv, ray_combined_bounce)
            end

            # This is duplicate code; I don't know how to do it better
            tl = ray_df_to_tl(ray_combined_df, drange, ddepth)
        else
            if isnothing(freq)
                error("Need freq if computing transmission loss without ray tracing")
            end
            tl = spherical_tl(attenuation_coef_thorp(freq), drange / 3)
        end
    end
    if tlonly
        println(tl)
        return nothing
    end

    # @df ray_combined_df plot(:range, :depth, group = :angle, yflip = true,
    #                          legend=false)

    if isnothing(sl) || isnothing(nlmean) || isnothing(nlsd) || isnothing(dt)
        error("sl, nlmean, nlsd, dt not all set when attempting to predict with
              passive sonar equation.")
    end

    son_set = sonar_passive(sl, tl, Normal(nlmean, nlsd), di, dt)
    println(detection_prob(son_set))
end

if !isinteractive()
    main(; parsed_args...)
    
    exit()
end

# Parameters for running main if desired
test_args = Dict(
    :ddepth => 100,
    :drange => 1.2 * 10^5,
    :minangle => -1.0,
    :maxangle => 1.0,
    :stepangle => 1.0,
    :emtdepth => 1000,
    :tgtrange => 30 * 6000,
    :maxdepth => 14000,
    :freq => 100,
    :lelements => Unsigned(20),
    :lspacing => 10,
    :svpmat => [# depth    velocity 
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
               ].* 3.28084,
    :sl => 100,
    :nlmean => 72,
    :nlsd => 2,
    :dt => 15
)

# main(; test_args...)

