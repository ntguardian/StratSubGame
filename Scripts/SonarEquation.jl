#!/bin/julia
# SonarEquation.jl
# 2023-05-09
# curtis
# This is a one-line description of the file.

# ArgParse: A package for handling command line arguments
using ArgParse

# STRUCTS ----------------------------------------------------------------------

abstract type Sonar end

struct sonar_noise <: Sonar
    sl :: Float64
    tl :: Float64
    ts :: Float64
    nl :: Float64
    di :: Float64
    dt :: Float64
end

struct sonar_passive <: Sonar
    sl :: Float64
    tl :: Float64
    nl :: Float64
    di :: Float64
    dt :: Float64
end

# FUNCTIONS --------------------------------------------------------------------

function freq_to_wavelength(velocity :: Float64, freq :: Float64) :: Float64
    velocity / freq
end

function wavelength_to_freq(velocity   :: Float64,
                            wavelength :: Float64) :: Float64
    freq_to_wavelength(velocity = velocity, freq = wavelength)
end

function piston_di(diameter :: Float64, wavelength :: Float64) :: Float64
    20 * (log10(π) + log10(diameter) - log10(wavelength))
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
end

main()


