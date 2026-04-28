# TrixiAero.jl

**TrixiAero.jl** is a package for high-fidelity aerodynamic simulations built on top of [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).
It extends Trixi's capabilities with specialized analysis callbacks and surface integral computations for aerodynamic applications.

## Features

- **Surface pressure and friction coefficients**: Compute and save pointwise aerodynamic coefficients along boundaries
- **Extended analysis capabilities**: Additional analysis callbacks tailored for aerodynamic simulations
- **Seamless Trixi.jl integration**: Built directly on Trixi's DG framework and mesh infrastructure

## Installation

TODO: Register
TrixiAero.jl is registered in the Julia General Registry. Install it using:

```julia
using Pkg
Pkg.add("TrixiAero")
```

## Quick Start

Take a look at the examples

## Documentation

See the [API Reference](@ref) for detailed documentation on available callbacks and functions.
