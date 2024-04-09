module SimpleSonar

using Distributions
using Interpolations
using DataFrames
import Base.rand

export ray_df_to_tl
export piston_di
export line_di
export attenuation_coef_thorp
export spherical_tl
export sonar_threshold
export detection_prob
export Sonar
export sonar_noise
export sonar_passive
include("SonarEquationTools.jl")

export svp
export piston_di
export line_di
export attenuation_coef_thorp
export spherical_tl
export sonar_threshold
export sonar_threshold
export detection_prob
export raytrace_tl
export depth_slice_idx
export svp_to_interpolation
export svp_refine
export DataFrame
include("SVP.jl")

export raytrace
export raytrace_step
export snell
export raytrace_angle_df
export ray_position_above_df
export get_containing_rays_df
export raytrace
export raytrace_step
export snell
export raytrace_angle_df
export ray_position_above_df
export get_containing_rays_df
include("Raytrace.jl")

export freq_to_wavelength
export wavelength_to_freq
include("Math.jl")

end
