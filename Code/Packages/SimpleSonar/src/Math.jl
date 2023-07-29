# Math.jl
# 2023-07-27
# curtis
# Math functions used when working with sonar equations

# FUNCTIONS --------------------------------------------------------------------

@doc raw"""
    freq_to_wavelength(velocity :: Real, freq :: Real) :: Real

Convert frequency to wavelength

Returns the result of ``λ=v/f``, where ``λ`` is the wavelength of a sound wave, ``v``
is the velocity of sound in water, and ``f`` is the frequency of the sound wave.

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

Returns the result of ``f=v/λ``, where ``λ`` is the wavelength of a sound wave, ``v``
is the velocity of sound in water, and ``f`` is the frequency of the sound wave.

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

