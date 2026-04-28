"""
    examples_dir()

Return the directory where the example files provided with TrixiAero.jl are located. If TrixiAero.jl is
installed as a regular package (with `]add TrixiAero`), these files are read-only and should *not* be
modified. To find out which files are available, use, e.g., `readdir`:

# Examples
```@example
readdir(examples_dir())
```
"""
examples_dir() = pkgdir(TrixiAero, "examples")
