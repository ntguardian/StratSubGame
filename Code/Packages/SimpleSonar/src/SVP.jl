# SVP.jl
# 2023-07-27
# curtis
# Tools for working with sound velocity profiles (SVPs)

# STRUCTS ----------------------------------------------------------------------

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
    depth_slice_idx(svp :: svp, depth :: Real) :: 

Get index of depth slices

The returned integer corresponds to the correct depth slice of the [`svp`](#SimpleSonar.svp)
object. The input depth can exceed the lowest depth in `svp`, but not be too
shallow.

...
# Arguments
- `svp :: svp`: Sound velocity profile
- `depth :: Real`: Input depth
...

See also [`svp`](#SimpleSonar.svp)

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

"""
    svp_to_interpolation(svp_obj :: svp, itp_type :: InterpolationType,
                              ext_type :: BoundaryCondition) :: Any

Turn a [`svp`](#SimpleSonar.svp) into an interpolation function

This function's utility comes from its ability to take limited sound velocity
profile data and convert it into a function that allows for a better selection
of slices for raytracing. That allows for feeding more fine grained data to a
raytracing algorithm.

...
# Arguments
- `svp_obj :: svp`: [`svp`](#SimpleSonar.svp) `struct` to be turned into an interpolation
- `itp_type`: Set how interpolation will be done
- `ext_type`: Set how extrapolation will be done
...

See also [`svp`](#SimpleSonar.svp)

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

Increase the resolution of a [`svp`](#SimpleSonar.svp) via interpolation and extrapolation

The resulting [`svp`](#SimpleSonar.svp) will have a constant distance between depths and may
extend to otherwise unseen depths.

...
# Arguments
- `svp_obj :: svp`: [`svp`](#SimpleSonar.svp) object to interpolate
- `Δ :: Real`: Difference between depth levels in new [`svp`](#SimpleSonar.svp)
- `min_depth :: Real = 0`: Minimum depth of new [`svp`](#SimpleSonar.svp)
- `max_depth :: Real = svp_obj.depth[end]`: Maximum depth of new [`svp`](#SimpleSonar.svp)
- `itp_type = Gridded(Linear())`: Set how interpolation will be done
- `ext_type = Linear()`: Set how extrapolation will be done
...

See also [`svp`](#SimpleSonar.svp), [`svp_to_interpolation`](#SimpleSonar.svp_to_interpolation)

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
    DataFrame(svp_obj :: svp) :: DataFrame

Convert [`svp`](#SimpleSonar.svp) structs to `DataFrame`s.

...
# Arguments
- `svp_obj :: svp`: The SVP to convert
...

See also [`svp`](#SimpleSonar.svp)

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

