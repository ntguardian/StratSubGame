# Raytrace.jl
# 2023-07-27
# curtis
# Tools for raytracing

# FUNCTIONS --------------------------------------------------------------------

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
the angle of the ray. If ``d_1`` is the first depth, ``d_2``, the second depth with
``d_1 < d_2``, the ray enters at angle ``θ`` and at horizontal position ``p``, we can
determine the ray's new position and which depth slice it will appear. If ``θ >
0``, the ray is travelling from ``d_2`` to ``d_1``; the ray's starting position is
``(p, d_2)``, and its new position will be

``\left(\cot(θ)(d_2 - d_1) + p, d_1\right).``

If the ray has an angle ``θ < 0``, it is then starting from position ``(p, d_1)``,
and will exit at

``\left(-\cot(θ)(d_2 - d_1) + p, d_2\right).``

The tuples above are returned by this function.

Exceptional cases include ``θ = 0`` and ``\left|θ\right| > π/2``. For horizontal
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

This essentially solves ``cos(θ_1)/v_1 = cos(θ_2)/v_2`` for ``θ_2``, where ``v_1`` is
the sound velocity of the depth layer being left, ``θ_1`` the angle at which the
sound ray is leaving, and ``v_2`` is the velocity of sound in the new layer.

...
# Arguments
- `angle :: Real`: Ray entering angle; must be between ``-π/2`` and ``π/2``
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

