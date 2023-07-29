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
using SimpleSonar

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

