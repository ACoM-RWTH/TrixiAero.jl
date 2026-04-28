# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Abstract base type used for dispatch of `analyze` for quantities
# requiring gradients of the velocity field.
abstract type VariableViscous end

include("analysis_surface_integral_2d.jl")
include("analysis_surface_integral_3d.jl")
end # muladd
