#!/bin/julia
# SSBNPatrolASWDES.jl
# 2023-07-11
# curtis
# A discrete event simulation for examining the probability that a SSBN would be
# detected by a listening MPRA and listening SOSUS arrays

# ArgParse: A package for handling command line arguments
using ArgParse

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

if !isinteractive()
    parsed_args = parse_commandline()
end

# PACKAGES ---------------------------------------------------------------------

using Distributions
using Distances
using LinearAlgebra
using ResumableFunctions
using ConcurrentSim
using .Threads
using Random

import Base: length, rand, rand!, +, -
import Distributions: pdf
import Statistics: mean

# STRUCTS ----------------------------------------------------------------------

"""
    Arena

Defines arenas for submarine simulation

While slightly overengineered, an arena is an important part of determining if
objects and positions are in the arena or out of the arena. The arena also
ensures proper setting. It is very simple: it only needs a center and a radius,
so all arenas are circular.

...
# Fields
- `center::Vector{<:Real}`: The center of the arena
- `radius::Real`: The radius of the arena
...
"""
struct Arena
    center::Vector{<:Real}
    radius::Real
    function Arena(center::Vector{<:Real}, radius::Real)
        if length(center) ≠ 2
            error("Arenas must be two-dimensional")
        end
        if radius ≤ 0
            error("Arena radius must be positive")
        end
        new(center, radius)
    end
end

@doc raw"""
    ArcUniform

A multivariate random variable distributed uniformly over an arc

Mathematically the random variable can be constructed as follows: Let $Θ ~
UNIF(θ_{min}, θ_{max})$, then the arc uniform random variable $X$ is equal in
distribution to the random vector
$$\left(\begin{matrix}c_x + r \cos(Θ) \\ c_y + r \sin(Θ)\end{matrix}\right)$$
where $r$ is the radius of the arc's circle, $c_x$ is the $x$ first coordinate
of the circle's center, and $y$ is the second coordinate of the circle's center.

One could argue that calling this random variable a continuous random variable
is false since its support has Lebesgue measure zero. Perhaps a more appropriate
subtyping should be used. This will work for our purposes for now, though.

...
# Fields
- `center::Vector{<:Real}`: The center of the circle over which the random
                            variable is distributed
- `radius::Real`: The radius of the circle over which the random variable is
                  distributed
- `θ_min::Real`: The minimum angle of the arc, in radians
- `θ_max::Real`: The maximum angle of the arc, in radians
...

See also [`MultivariateDistribution`](@MultivariateDistribution)
"""
struct ArcUniform <: ContinuousMultivariateDistribution
    center::Vector{<:Real}
    radius::Real
    θ_min::Real
    θ_max::Real
    function ArcUniform(center::Vector{<:Real}, radius::Real, θ_min::Real,
                        θ_max::Real)
        if θ_min > θ_max
            error("Cannot have θ_min > θ_max")
        end
        if θ_min < 0 || θ_max > 2π
            error("Must have angles bounded: 0 ≤ θ_min ≤ θ_max ≤ 2π")
        end
        if length(center) ≠ 2
            error("Arenas must be two-dimensional")
        end
        if radius ≤ 0
            error("Arena radius must be positive")
        end
        new(center, radius, θ_min, θ_max)
    end
end

# FUNCTIONS --------------------------------------------------------------------

"""
inarena(arena::Arena, position::Vector{<:Real})::Bool

Check that a position is in the arena

If not, return `false`; otherwise, return `true`.

...
# Arguments
- `arena::Arena`: The arena for which position should be checked
- `position::Vecotr{<:Real}`: The position to check
...

See also [`Arena`](@Arena)

# Examples
```jldoctest
julia> inarena(Arena([0;0], 20), [20;1])
false
julia> inarena(Arena([0;0], 20), [10;1])
true
```
"""
function inarena(arena::Arena, position::Vector{<:Real})::Bool
    if length(position) == 2 && euclidean(position, arena.center) ≤ arena.radius
        return true
    else
        return false
    end
end

"""
    ssbn(env::Environment,
         start_pos::Vector{Real},
         pos_dist::ContinuousMultivariateDistribution,
         speed::Real,
         arena::Arena)

Create SSBN simulation instance in a DES

The basic logic of a boomer is:

* Enter the patrol area (the simulation starts)
* If the patrol is not over, find a new position

Hence the SSBN's behavior is fairly simple: move to a new random location upon
reaching the most recent random location.

...
# Arguments
- `env::Environment`: Simulation environment
- `start_pos::Vector{Real}`: Starting position of the submarine
- `pos_dist::ContinuousMultivariateDistribution`: Distribution describing new
                                                  proposed locations to which
                                                  the submarine will travel (may
                                                  be rejected)
- `speed::Real`: Speed of the submarine
- `arena::Arena`: Arena in which the sub travels
...

See also # TODO: SEE ALSO

# Examples
```jldoctest
julia> arn = Arena([0; 0], 10)
[...]
julia> sim = Simulation()
[...]
julia> @process ssbn(sim, [0; 0], MvNormal([0; 0], [1 0; 0 1]), 10/24, arn)
[...]
julia> run(sim, 100)
[...]
```
"""
@resumable function ssbn(env::Environment,
                         start_pos::Vector{<:Real},
                         pos_dist::ContinuousMultivariateDistribution,
                         speed::Real,
                         arena::Arena)
    if !iinarena(start_pos)
        error("start_pos must be in arena")
    end

    pos = start_pos
    destination = pos
    dist_to_travel = 0
    time_to_arrival = 0

    while true
        println("SSBN position: ($(pos[1]), $(pos[2])) at t=$(now(env))")

        destination = sample_pos_in_arena(pos_dist + pos, arena)
        dist_to_travel = euclidean(pos, destination)
        time_to_arrival = dist_to_travel / speed
        @yield timeout(env, time_to_arrival)
        pos = destination
    end
end

"""
    sample_pos_in_arena(pos_dist::ContinuousMultivariateDistribution,
                        arena::Arena;
                        max_iter::Unsigned=typemax(UInt128))::Vector{<:Real}

Sample position until a position in the arena has been found

Ths is simple rejection sampling, drawing from the test distribution until one
in the arena is found. Not very efficient, but should be acceptable for our
purposes.

...
# Arguments
- `pos_dist::ContinuousMultivariateDistribution`: Sampling distribution for
                                                  generating positions
- `arena::Arena`: The arena to find a point in
- `max_iter::Unsigned`: Maximum number of iterations to try before throwing an
                        error (currently the largest number possible, so not
                        really a bound unless set)
...

See also [`Arena`](@Arena)

# Examples
```jldoctest
julia> sample_pos_in_arena(MvNorm([0;0], [1 0; 0 1]), Arena([0;0], 10))
[...]
```
"""
function sample_pos_in_arena(pos_dist::ContinuousMultivariateDistribution,
                             arena::Arena;
                             max_iter::Unsigned=typemax(UInt128))::Vector{<:Real}
    for iter in 1:max_iter
        pos_proposal = rand(pos_dist)
        if inarena(arena, pos_proposal)
            return pos_proposal
        end
    end
    error("Maximum iterations reached; never found a position in arena")
end

@doc raw"""
    sample_arc_uniform(r::Real, θ_min::Real, θ_max::Real)::Vector{<:Real}

Sample uniformly at random the position on a circular arc

The result is a two-dimensional vector.

Sampling is done by first picking Θ such that $Θ ~ \text{UNIF}(θ_{\text{min}},
θ_{\text{max}})$, the $x$ coordinate is $r \cos(Θ)$, and the $y$ coordinate is
$r \sin(Θ)$.

...
# Arguments
- `r::Real`: Radius of the circular arc to sample in
- `θ_min::Real`: The minimum angle, in radians
- `θ_max::Real`: The maximum angle, in radians
...

See also [`ssbn`](@ssbn)

# Examples
```jldoctest
julia> sample_arc_uniform(10, 0, π/4)
[...]
```
"""
function sample_arc_uniform(r::Real, θ_min::Real, θ_max::Real)::Vector{<:Real}
    # TODO: curtis: FUNCTION BODY -- Wed 12 Jul 2023 01:49:34 AM EDT
end

"""
    length(d::ArcUniform)::Int8

The length of each sample

Since this multivariate random variable is always two-dimensional, this always
returns two.

...
# Arguments
- d::ArcUniiform : The random variable for which the length is needed
...

See also [`size`](@size)

# Examples
```jldoctest
julia> X = ArcUniform([0; 0], 1, 0, 2π)
ArcUniform(center=[0, 0], radius=1, θ_min=0, θ_max=6.283185307179586)
julia> length(X)
2
```
"""
function length(d::ArcUniform)::Int8
    length(d.center)
end

"""
    size(d::ArcUniform)::Tuple{Int8}

The sample size of an [`ArcUniform`](@ArcUniform)-distributed random variable

This is just `(length(d),)`.

...
# Arguments
- `d::ArcUniform`: A distribution
...

See also [`length`](@length)

# Examples
```jldoctest
julia> X = ArcUniform([0; 0], 1, 0, 2π)
ArcUniform(center=[0, 0], radius=1, θ_min=0, θ_max=6.283185307179586)
julia> size(X)
(2,)
```
"""
function size(d::ArcUniform)::Tuple{Int8}
    (length(d),)
end

"""
    eltype(d::ArcUniform)::DataType

The default element type of a sample

This is the type of elements generated by the `rand` method.

...
# Arguments
- `d::ArcUniform`: The distribution
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([0; 0], 1, 0, 2π)
ArcUniform(center=[0, 0], radius=1, θ_min=0, θ_max=6.283185307179586)
julia> eltype(X)
Float64
```
"""
function eltype(d::ArcUniform)::DataType
    Float64
end

"""
    insupport(d::ArcUniform, x::AbstractArray)::Bool

TODO: BRIEF DESCRIPTION

TODO: EXTENDED DESCRIPTION (OPTIONAL)

...
# Arguments
- `d::ArcUniform`: The distribution
- `x::AbstractArray`: The vector or matrix with the points to check for support
                      membership
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> insupport(X, [2 ; 2])
true
julia> insupport(X, [1 ; 2])
false
julia> insupport(X, [4 ; 2])
false
julia> insupport(X, [2 2 + 1/2; 2 1 + √3/2])
false
```
"""
function insupport(d::ArcUniform, x::AbstractArray)::Bool
    if x isa Vector{<:Real}
        return insupport_vec(d, x)    # Implemented elsewhere as a method
    end

    all((insupport_vec(d, convert(Vector, u)) for u in eachcol(x)))
end

"""
    insupport_vec(d::ArcUniform, x::Vector{<:Real})::Bool

Check if a vector is in the support of an arc uniform random variable

A point is in the support of the arc uniform random variable if the distance
between the point and the center of the arc's circle is the radius and if the
angle of the point is between the minimum and maximum angles of the arc.

...
# Arguments
- `d::ArcUniform`: The distribution
- `x::Vector{<:Real}`: The vector with the points to check for support
                       membership
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> insupport_vec(X, [2 ; 2])
true
julia> insupport_vec(X, [1 ; 2])
false
julia> insupport_vec(X, [4 ; 2])
false
```
"""
function insupport_vec(d::ArcUniform, x::Vector{<:Real})::Bool
    pre_θ = acos((x[1] - d.center[1])/euclidean(x, d.center))
    θ = x[2] > d.center[2] ? pre_θ : 2π - pre_θ

    euclidean(x, d.center) ≈ d.radius && d.θ_min ≤ θ ≤ d.θ_max
end

@doc raw"""
    pdf(d::ArcUniform, x::AbstractArray)

Probability density function of an arc uniform random variable

This is not really a probability density function because the Lebesgue measure
of the support is zero, but set that aside for a moment. If a point is in the
support of the random variable, the value of the probability density function is
$\frac{1}{θ_{max} - θ_{min}}$, and zero otherwie.

...
# Arguments
- `d::ArcUniform`: The distribution
- `x::AbstractArray`: The vector or matrix with the points to check for support
                      membership
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> pdf(X, [2 ; 2])
1.2732395447351628
julia> pdf(X, [1 ; 2])
0.0
julia> pdf(X, [4 ; 2])
0.0
julia> pdf(X, [2 2 + 1/2; 2 1 + √3/2])
2-element Vector{Float64}:
 1.2732395447351628
 1.2732395447351628
```
"""
function pdf(d::ArcUniform, x::AbstractArray{T, M} where T<:Real) where M
    if x isa Vector{<:Real}
        return pdf_vec(d, x)    # Implemented elsewhere as a method
    end

    Vector([pdf_vec(d, convert(Vector, u)) for u in eachcol(x)])
end

"""
    pdf_vec(d::ArcUniform, x::Vector{<:Real})::Float64

Probability density function of an arc uniform random variable

This is not really a probability density function because the Lebesgue measure
of the support is zero, but set that aside for a moment. If a point is in the
support of the random variable, the value of the probability density function is
$\frac{1}{θ_{max} - θ_{min}}$, and zero otherwie.

...
# Arguments
- `d::ArcUniform`: The distribution
- `x::Vector{<:Real}`: The vector of the point at which to evaluate the PDF
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> pdf_vec(X, [2 ; 2])
1.2732395447351628
julia> pdf_vec(X, [1 ; 2])
0.0
julia> pdf_vec(X, [4 ; 2])
0.0
```
"""
function pdf_vec(d::ArcUniform, x::Vector{<:Real})::Float64
    if !insupport_vec(d, x)
        return 0
    else
        return 1/(d.θ_max - d.θ_min)
    end
end

@doc raw"""
    logpdf(d::ArcUniform, x::AbstractArray)

Log probability density function of an arc uniform random variable

This is simply log applied to the pdf as handled by [`pdf`](@pdf)

...
# Arguments
- `d::ArcUniform`: The distribution
- `x::AbstractArray`: The vector or matrix with the points to check for support
                      membership
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> logpdf(X, [2 ; 2])
0.24156447527049052
julia> logpdf(X, [1 ; 2])
-Inf
julia> logpdf(X, [4 ; 2])
-Inf
julia> logpdf(X, [2 2 + 1/2; 2 1 + √3/2])
2-element Vector{Float64}:
 0.24156447527049052
 0.24156447527049052
```
"""
function logpdf(d::ArcUniform, x::AbstractArray{T, M} where T<:Real) where M
    log.(pdf(d, x))
end

@doc raw"""
    mean(d::ArcUniform)::Vector{Float64}

The mean of an arc uniform random variable.

For an arc uniform random variable with circle center at
$$\left(\begin{matrix}c_1 \\ c_2\end{matrix}\right),$$
with radius $r$ and between angles $θ_{min}≤θ_{max}$, the mean is:
$$\frac{r}{θ_{max}-θ_{min}}\left(\begin{matrix}\sin(θ_{max}) - \sin(θ_{min}) +
c_1 \\ \cos(θ_{min}) - \cos(θ_{max}) +
c_2\end{matrix}\right).$$

...
# Arguments
- `d::ArcUniform`: The distribution
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> mean(X)
2-element Vector{Float64}:
 2.919402318048382
 2.173555860892269
```
"""
function mean(d::ArcUniform)::Vector{Float64}
    d.radius / (d.θ_max - d.θ_min) * [sin(d.θ_max) - sin(d.θ_min) + d.center[1]
                                      ; cos(d.θ_min) - cos(d.θ_max) + d.center[2]]
end

@doc raw"""
    var(d::ArcUniform)::Vector{Float64}

The variance of an arc uniform random variable's components

For an arc uniform random variable with center at the origin, it can be shown
that the second moment of the random variable's first component is
$$\frac{r^2}{2 (θ_{max} - θ_{min})} (\cos(θ_{max}) \sin(θ_{max}) - \cos(θ_{min}) -
\sin(θ_{min}) + θ_{max} - θ_{min}).$$
The second moment for the random variable's second component is
$$\frac{r^2}{2 (θ_{max} - θ_{min})} (\cos(θ_{min}) \sin(θ_{min}) - \cos(θ_{max}) -
\sin(θ_{max}) + θ_{max} - θ_{min}).$$
Combining these with the expression for the mean provides a formula for the
variance.

...
# Arguments
- `d::ArcUniform`: The distribution
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> var(X)
2-element Vector{Float64}:
 0.04261837940312782
 0.0077404170450887655
```
"""
function var(d::ArcUniform)::Vector{Float64}
    cosdiff = cos(d.θ_max) - cos(d.θ_min)
    sindiff = sin(d.θ_max) - sin(d.θ_min)
    denlen = d.θ_max - d.θ_min
    r = d.radius
    [r^2 / (2 * denlen) * (cos(d.θ_max) * sin(d.θ_max) - cos(d.θ_min) *
                           sin(d.θ_min) + denlen) - r^2 / denlen^2 * sindiff^2 ;
     r^2 / (2 * denlen) * (cos(d.θ_min) * sin(d.θ_min) - cos(d.θ_max) *
                           sin(d.θ_max) + denlen) - r^2 / denlen^2 * cosdiff^2]
end

@doc raw"""
    cov(d::ArcUniform)::Matrix{Float64}

The covariance matrix of an arc uniform random variable

The mean and component variance of the arc uniform random variable has been
discussed in the documentation of [`mean`](@mean) and [`var`](@var). For the
covariance, it can be shown that the expectation of the components' product is
$$\frac{r^2(\sin^2(θ_{max}) - \sin^2(θ_{min}))}{2 (θ_{max} - θ_{min})}.$$
The rest follows from the usual calculations.

...
# Arguments
- `d::ArcUniform`: The distribution
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> cov(X)
2×2 Matrix{Float64}:
 0.0426184  0.654059
 0.654059   0.00774042
```
"""
function cov(d::ArcUniform)::Matrix{Float64}
    cosdiff = cos(d.θ_max) - cos(d.θ_min)
    sindiff = sin(d.θ_max) - sin(d.θ_min)
    denlen = d.θ_max - d.θ_min
    r = d.radius
    comp_cov = r^2 * (sin(d.θ_max)^2 - sin(d.θ_min)^2)/(2 * denlen) +
               r^2 / denlen^2 * cosdiff * sindiff
    comp_var = var(d)
    
    [comp_var[1] comp_cov ; comp_cov comp_var[2]]
end

"""
    cor(d::ArcUniform)::Matrix{Float64}

The correlation matrix of an arc uniform random variable

Computed in the usual way one would compute a correlation matrix.

...
# Arguments
- `d::ArcUniform`: The distribution
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> cor(X)
2×2 Matrix{Float64}:
  1.0       -0.960153
 -0.960153   1.0
```
"""
function cor(d::ArcUniform)::Matrix{Float64}
    cov_mat = cov(d)
    denom = sqrt(*(diag(cov_mat)...))
    
    [1 cov_mat[1, 2] / denom ; cov_mat[2, 1] / denom 1]
end

@doc raw"""
    entropy(d::ArcUniform, b::Real)::Float64

The entropy of an arc uniform distribution

The entropy of uniform random variables is easily computed since the density of
a uniform random variable is constant; it will be $-\log_b(θ_{max} - θ_{min})$.

...
# Arguments
- `d::ArcUniform`: The distribution
- `b::Real`: The base of the logarithm
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> entropy(X, 2)
0.3485038705276813
julia> entropy(X)
0.2415644752704905
```
"""
function entropy(d::ArcUniform, b::Real)::Float64
    -log(b, (d.θ_max - d.θ_min))
end

@doc raw"""
    entropy(d::ArcUniform)::Float64

The entropy of an arc uniform distribution

The entropy of uniform random variables is easily computed since the density of
a uniform random variable is constant; it will be $-\log(θ_{max} - θ_{min})$.

...
# Arguments
- `d::ArcUniform`: The distribution
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> entropy(X, 2)
0.3485038705276813
julia> entropy(X)
0.2415644752704905
```
"""
function entropy(d::ArcUniform)::Float64
    -log((d.θ_max - d.θ_min))
end

@doc raw"""
    rand(rng::AbstractRNG, d::ArcUniform)::Vector{Float64}

Generate a random arc uniform vector

Generate a random realization of an arc uniform random variable, via randomly
sampling an angle between $θ_{min}$ and $Θ_{max}$, letting the first coordinate
be the cosine of the random number and the second coordinate be the sine,
multiply by the radius, then add the center vector.

...
# Arguments
- `rng::AbstractRNG`: Random number generator
- `d::ArcUniform`: The distribution
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> rand(X)
[...]
```
"""
function rand(rng::AbstractRNG, d::ArcUniform)::Vector{Float64}
    Θ = d.θ_min + (d.θ_max - d.θ_min) * rand(rng)

    d.radius * [cos(Θ) ; sin(Θ)] + d.center
end

"""
    rand!(rng::AbstractRNG, d::ArcUniform, x::AbstractArray{<:Real})

Randomly sample from arc uniform random variables in place

This function will try to detect if `x` is a matrix with one column of length
two, and applies an appropriate transform to make sure one column or row
corresponds to one of the dimensions of the random variable, and the other the
remaining random variable, thus avoiding intermixing the dimensions of the arc
uniform random variable within a column. If `x` is not a matrix with either two
columns or two rows, it is on the user to ensure that the resulting random
variables make sense. The function iterates through the elements of `x` using
linear indexing (column major).

...
# Arguments
- `rng::AbstractRNG`: Random number generator
- `d::ArcUniform`: The distribution
- `x::AbstractArray{<:Real}`: The array to contain the random numbers
...

See also [`rand`](@!rand)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> mat = [0.0 0.0; 0.0 0.0; 0.0 0.0]
3×2 Matrix{Float64}:
 0.0  0.0
 0.0  0.0
 0.0  0.0
julia> rand!(MersenneTwister(), X, mat)
[...]
julia> altmat = transpose(mat)
[...]
julia> rand!(MersenneTwister(), X, altmat)
[...]
```
"""
function rand!(rng::AbstractRNG,
               d::ArcUniform,
               x::AbstractArray{<:Real})
    func_vec = [sin, cos]
    if ndims(x) == 2 && Base.size(x)[2] == 2
        trfm = transpose
    else
        trfm = identity
    end
    rnum = 0.0
    for i in eachindex(x)
        rnum = (mod(i, 2) == 0) ? rnum : rand(rng)
        @inbounds trfm(x)[i] = d.radius * func_vec[1 + mod(i, 2)](d.θ_min +
            (d.θ_max - d.θ_min) * rnum) + d.center[2 - mod(i, 2)]
    end
end

"""
    loglikelihood(d::ArcUniform, x::Matrix{<:Real})::Float64

Compute the log likelihood of an arc uniform random variable

The function will check whether the number of rows or the number of columns is
two, and use the result to determine whether rows or columns represent samples.
If neither is two, an error is thrown.

...
# Arguments
- `d::ArcUniform`: The distribution
- `x::Matrix{<:Real}`: The data
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> mat = [0.0 0.0; 0.0 0.0; 0.0 0.0]
3×2 Matrix{Float64}:
 0.0  0.0
 0.0  0.0
 0.0  0.0
julia> rand!(MersenneTwister(), X, mat)
[...]
julia> loglikelihood(X, mat)
[...]
```
"""
function loglikelihood(d::ArcUniform, x::Matrix{<:Real})::Float64
    if Base.size(x)[1] == 2
        itermode = eachcol
    elseif Base.size(x)[2] == 2
        itermode = eachrow
    else
        error("Either the rows or the columns must be 2")
    end

    llsum = 0.0
    for vec in itermode(x)
        llsum += logpdf(d, vec)[1]
    end

    llsum
end

"""
    +(d::ArcUniform, x::Vector{<:Real})::ArcUniform

Add an arc uniform random variable to a vector

Adding a constant vector to an arc uniform random variable results in a shift of
the random variable's center parameter

...
# Arguments
- `d::ArcUniform`: The random variable
- `x::Vector{<:Real}`: The vector shift
...

See also [`ArcUniform`](@ArcUniform)

# Examples
```jldoctest
julia> X = ArcUniform([2; 1], 1, π/4, π/2)
ArcUniform(center=[2, 1], radius=1, θ_min=0.7853981633974483, θ_max=1.5707963267948966)
julia> X + [1, 1]
[...]
julia> X - [1, 1]
[...]
julia> [1, 1] + X
[...]
```
"""
function +(d::ArcUniform, x::Vector{<:Real})::ArcUniform
    ArcUniform(d.center + x, d.radius, d.θ_min, d.θ_max)
end
+(x::Vector{<:Real}, d::ArcUniform)::ArcUniform = d + x
-(d::ArcUniform, x::Vector{<:Real})::ArcUniform = d + (-x)
# There is no x - d; I have not implemented a negative version of this random
# variable

# EXECUTABLE SCRIPT MAIN FUNCTIONALITY -----------------------------------------

function main( ; patrollength::Real=100, sims::Unsigned=1000)
    # Experimental code
    sim = Simulation()

    @threads for i in 1:sims
        simenv = Simulation()

        @process ssbn(simenv)

        run(simenv, patrollength)
    end
end

if !isinteractive()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("    $arg => $val")
    end
    main()
end


