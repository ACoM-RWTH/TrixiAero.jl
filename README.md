# AeroTrixi.jl
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://acom-rwth.github.io/AeroTrixi.jl/stable/)
[![CI](https://github.com/ACoM-RWTH/AeroTrixi.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ACoM-RWTH/AeroTrixi.jl/actions/workflows/ci.yml)
[![Codecov](https://codecov.io/gh/ACoM-RWTH/AeroTrixi.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ACoM-RWTH/AeroTrixi.jl)
[![Coveralls](https://coveralls.io/repos/github/ACoM-RWTH/AeroTrixi.jl/badge.svg?branch=main)](https://coveralls.io/github/ACoM-RWTH/AeroTrixi.jl?branch=main)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

**AeroTrixi.jl** is a package for high-fidelity aerodynamic simulations built on top of [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).
Currently, it extends Trixi's capabilities with a specialized analysis callback for aerodynamic applications.
In the future, we plan to add more features for nonideal and rarefied gases, for instance.

## Features

- **Surface pressure and friction coefficients**: Compute and save pointwise aerodynamic coefficients along boundaries

## Installation

AeroTrixi.jl is registered in the Julia General Registry. Install it using:

```julia
using Pkg
Pkg.add("AeroTrixi")
```

## Quick Start

To get started, it is best to take a look at the examples.
Currently, the only additional functionality by AeroTrixi.jl is the extended 
`AnalysisCallback` which can be used to compute pointwise aerodynamic coefficients along boundaries,
such as `SurfacePressureCoefficient` and `SurfaceFrictionCoefficient`.
