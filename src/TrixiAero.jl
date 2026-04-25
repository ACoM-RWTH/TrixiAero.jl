"""
    TrixiAero

High-fidelity aerodynamic simulations with Trixi.jl.
"""
module TrixiAero

using Trixi
# import bunch of non-exported stuff from Trixi to avoid writing `Trixi.` everywhere
using Trixi: @printf, @sprintf, DiscreteCallback, AbstractEquations,
             AbstractSemidiscretization

using MuladdMacro: @muladd
using StaticArrays: SVector, SMatrix, SArray, MVector, MArray

include("callbacks_step/callbacks_step.jl")

export AnalysisSurfacePointwise

end
