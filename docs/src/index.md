# TrixiAero.jl

**TrixiAero.jl** is a package for high-fidelity aerodynamic simulations built on top of [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).
Currently, it extends Trixi's capabilities with a specialized analysis callback for aerodynamic applications.
In the future, we plan to add more features for nonideal and rarefied gases, for instance.

## Features

- **Surface pressure and friction coefficients**: Compute and save pointwise aerodynamic coefficients along boundaries
- **Seamless Trixi.jl integration**: Built directly on Trixi's DG framework and mesh infrastructure

## Installation

TrixiAero.jl is registered in the Julia General Registry. Install it using:

```julia
using Pkg
Pkg.add("TrixiAero")
```

## Quick Start

To get started, it is best to take a look at the examples.
Currently, the only additional functionality by TrixiAero.jl is the extended 
`AnalysisCallback` which can be used to compute pointwise aerodynamic coefficients along boundaries,
such as `SurfacePressureCoefficient` and `SurfaceFrictionCoefficient`.

## Documentation

See the [API Reference](@ref) for documentation on available/extended functions.
