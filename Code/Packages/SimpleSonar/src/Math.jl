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
- `freq :: Real`: Wave frequency
...

See also [`wavelength_to_freq`](#SimpleSonar.wavelength_to_freq)

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

See also [`freq_to_wavelength`](#SimpleSonar.freq_to_wavelength)

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

