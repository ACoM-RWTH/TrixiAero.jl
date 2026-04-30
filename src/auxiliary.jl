"""
    examples_dir()

Return the directory where the example files provided with AeroTrixi.jl are located. If AeroTrixi.jl is
installed as a regular package (with `]add AeroTrixi`), these files are read-only and should *not* be
modified. To find out which files are available, use, e.g., `readdir`:

# Examples
```@example
readdir(examples_dir())
```
"""
examples_dir() = pkgdir(AeroTrixi, "examples")
