# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

"""
    AnalysisCallback(semi; interval=0,
                           save_analysis=false,
                           output_directory="out",
                           analysis_filename="analysis.dat",
                           extra_analysis_errors=Symbol[],
                           extra_analysis_integrals=(),
                           analysis_pointwise=())

Analyze a numerical solution every `interval` time steps and print the
results to the screen. If `save_analysis`, the results are also saved in
`joinpath(output_directory, analysis_filename)`.

Additional errors can be computed, e.g. by passing
`extra_analysis_errors = (:l2_error_primitive, :linf_error_primitive)`
or `extra_analysis_errors = (:conservation_error,)`.

If you want to omit the computation (to safe compute-time) of the `default_analysis_errors`, specify
`analysis_errors = Symbol[]`.
Note: `default_analysis_errors` are `:l2_error` and `:linf_error` for all equations.
If you want to compute `extra_analysis_errors` such as `:conservation_error` solely, i.e.,
without `:l2_error, :linf_error` you need to specify
`analysis_errors = [:conservation_error]` instead of `extra_analysis_errors = [:conservation_error]`.

Further scalar functions `func` in `extra_analysis_integrals` are applied to the numerical
solution and integrated over the computational domain. Some examples for this are
`entropy`, `energy_kinetic`, `energy_internal`, and `energy_total`.
You can also write your own function with the same signature as the examples listed above and
pass it via `extra_analysis_integrals`.
The default `analysis_integrals` is `(entropy_timederivative,)`.
You can also request `extra_analysis_integrals` such as `LiftCoefficientPressure` or
`DragCoefficientPressure` by constructing an `AnalysisSurfaceIntegral` with one of 
the previously mentioned functions.

Similarly, pointwise, i.e., per quadrature/interpolation point, quantities such at
`SurfacePressureCoefficient` or `SurfaceFrictionCoefficient` can be computed.
Instances of these need to be passed into `AnalysisSurfacePointwise` which is then in turn
passed to `analysis_pointwise`.

See the developer comments about `Trixi.analyze`, `Trixi.pretty_form_utf`, and
`Trixi.pretty_form_ascii` for further information on how to create custom analysis quantities.

In addition, the analysis callback records and outputs a number of quantities that are useful for
evaluating the computational performance, such as the total runtime, the performance index
(time/DOF/rhs!), the time spent in garbage collection (GC), or the current memory usage (alloc'd
memory).
"""
mutable struct AnalysisCallback{Analyzer, AnalysisIntegrals, AnalysisPointwise,
                                InitialStateIntegrals, Cache}
    start_time::Float64
    start_time_last_analysis::Float64
    ncalls_rhs_last_analysis::Int
    start_gc_time::Float64
    const interval::Int
    const save_analysis::Bool
    const output_directory::String
    const analysis_filename::String
    const analyzer::Analyzer
    const analysis_errors::Vector{Symbol}
    const analysis_integrals::AnalysisIntegrals
    const analysis_pointwise::AnalysisPointwise # TrixiAero addition
    initial_state_integrals::InitialStateIntegrals
    const cache::Cache
end

function Base.show(io::IO, ::MIME"text/plain",
                   cb::DiscreteCallback{<:Any, <:AnalysisCallback})
    @nospecialize cb # reduce precompilation time

    if get(io, :compact, false)
        show(io, cb)
    else
        analysis_callback = cb.affect!

        setup = Pair{String, Any}["interval" => analysis_callback.interval,
                                  "analyzer" => analysis_callback.analyzer]
        for (idx, error) in enumerate(analysis_callback.analysis_errors)
            push!(setup, "│ error " * string(idx) => error)
        end
        for (idx, integral) in enumerate(analysis_callback.analysis_integrals)
            push!(setup, "│ integral " * string(idx) => integral)
        end
        # TrixiAero addition
        for (idx, quantity) in enumerate(analysis_callback.analysis_pointwise)
            push!(setup, "│ pointwise " * string(idx) => quantity)
        end
        push!(setup,
              "save analysis to file" => analysis_callback.save_analysis ? "yes" : "no")
        if analysis_callback.save_analysis
            push!(setup, "│ filename" => analysis_callback.analysis_filename)
            push!(setup,
                  "│ output directory" => abspath(normpath(analysis_callback.output_directory)))
        end
        summary_box(io, "AnalysisCallback", setup)
    end
end

# The functions in "analysis_trixi.jl" are copy-paste from the upstream Trixi.jl repo.
# They add no new functionality, but are required to make the herein "overloaded" `AnalysisCallback` work with Trixi.jl.
include("analysis_trixi.jl")

# This is the actual constructor
function AnalysisCallback(mesh, equations::AbstractEquations, solver, cache;
                          interval = 0,
                          save_analysis = false,
                          output_directory = "out",
                          analysis_filename = "analysis.dat",
                          extra_analysis_errors = Symbol[],
                          analysis_errors = union(default_analysis_errors(equations),
                                                  extra_analysis_errors),
                          extra_analysis_integrals = (),
                          analysis_integrals = union(default_analysis_integrals(equations),
                                                     extra_analysis_integrals),
                          analysis_pointwise = (),
                          RealT = real(solver),
                          uEltype = eltype(cache.elements),
                          kwargs...)
    # Decide when the callback is activated.
    # With error-based step size control, some steps can be rejected. Thus,
    #   `integrator.iter >= integrator.stats.naccept`
    #    (total #steps)       (#accepted steps)
    # We need to check the number of accepted steps since callbacks are not
    # activated after a rejected step.
    condition = (u, t, integrator) -> interval > 0 &&
        (integrator.stats.naccept % interval == 0 || isfinished(integrator))

    analyzer = SolutionAnalyzer(solver; kwargs...)
    cache_analysis = create_cache_analysis(analyzer, mesh, equations, solver, cache,
                                           RealT, uEltype)

    analysis_callback = AnalysisCallback(0.0, 0.0, 0, 0.0,
                                         interval, save_analysis, output_directory,
                                         analysis_filename,
                                         analyzer,
                                         analysis_errors, Tuple(analysis_integrals),
                                         Tuple(analysis_pointwise), # TrixiAero addition
                                         SVector(ntuple(_ -> zero(uEltype),
                                                        Val(nvariables(equations)))),
                                         cache_analysis)

    return DiscreteCallback(condition, analysis_callback,
                            save_positions = (false, false),
                            initialize = initialize!)
end

# This method is just called internally from `(analysis_callback::AnalysisCallback)(integrator)`
# and serves as a function barrier. Additionally, it makes the code easier to profile and optimize.
function (analysis_callback::AnalysisCallback)(io, du, u, u_ode, t, semi, iter)
    @unpack analyzer, analysis_errors, analysis_integrals, analysis_pointwise = analysis_callback
    cache_analysis = analysis_callback.cache
    _, equations, _, _ = mesh_equations_solver_cache(semi)

    # Calculate and print derived quantities (error norms, entropy etc.)
    # Variable names required for L2 error, Linf error, and conservation error
    if any(q in analysis_errors
           for q in (:l2_error, :linf_error, :conservation_error, :residual)) &&
       mpi_isroot()
        print(" Variable:     ")
        for v in eachvariable(equations)
            @printf("  %-15s", varnames(cons2cons, equations)[v])
        end
        println()
    end

    if :l2_error in analysis_errors || :linf_error in analysis_errors ||
       :l1_error in analysis_errors
        # Calculate L2/Linf errors
        l2_error, linf_error, l1_error = calc_error_norms(u_ode, t, analyzer, semi,
                                                          cache_analysis)

        if mpi_isroot()
            # L2 error
            if :l2_error in analysis_errors
                print(" L2 error:    ")
                for v in eachvariable(equations)
                    @printf("  % 10.8e", l2_error[v])
                    print(io, " ", l2_error[v])
                end
                println()
            end

            # Linf error
            if :linf_error in analysis_errors
                print(" Linf error:  ")
                for v in eachvariable(equations)
                    @printf("  % 10.8e", linf_error[v])
                    print(io, " ", linf_error[v])
                end
                println()
            end

            # L1 error
            if :l1_error in analysis_errors
                print(" L1 error:    ")
                for v in eachvariable(equations)
                    @printf("  % 10.8e", l1_error[v])
                    print(io, " ", l1_error[v])
                end
                println()
            end
        end
    end

    # Conservation error
    if :conservation_error in analysis_errors
        @unpack initial_state_integrals = analysis_callback
        state_integrals = integrate(u_ode, semi)

        if mpi_isroot()
            print(" |∑U - ∑U₀|:  ")
            for v in eachvariable(equations)
                err = abs(state_integrals[v] - initial_state_integrals[v])
                @printf("  % 10.8e", err)
                print(io, " ", err)
            end
            println()
        end
    end

    # Residual (defined here as the vector maximum of the absolute values of the time derivatives)
    if :residual in analysis_errors
        mpi_print(" max(|Uₜ|):   ")
        for v in eachvariable(equations)
            # Calculate maximum absolute value of Uₜ
            res = maximum(abs, view(du, v, ..))
            if mpi_isparallel()
                # TODO: Debugging, here is a type instability
                # Base.max instead of max needed, see comment in src/auxiliary/math.jl
                global_res = MPI.Reduce!(Ref(res), Base.max, mpi_root(), mpi_comm())
                if mpi_isroot()
                    res::eltype(du) = global_res[]
                end
            end
            if mpi_isroot()
                @printf("  % 10.8e", res)
                print(io, " ", res)
            end
        end
        mpi_println()
    end

    # L2/L∞ errors of the primitive variables
    if :l2_error_primitive in analysis_errors ||
       :linf_error_primitive in analysis_errors ||
       :l1_error_primitive in analysis_errors
        l2_error_prim, linf_error_prim, l1_error_prim = calc_error_norms(cons2prim,
                                                                         u_ode, t,
                                                                         analyzer,
                                                                         semi,
                                                                         cache_analysis)

        if mpi_isroot()
            print(" Variable:     ")
            for v in eachvariable(equations)
                @printf("  %-15s", varnames(cons2prim, equations)[v])
            end
            println()

            # L2 error
            if :l2_error_primitive in analysis_errors
                print(" L2 error prim.: ")
                for v in eachvariable(equations)
                    @printf("%10.8e   ", l2_error_prim[v])
                    print(io, " ", l2_error_prim[v])
                end
                println()
            end

            # L∞ error
            if :linf_error_primitive in analysis_errors
                print(" Linf error pri.:")
                for v in eachvariable(equations)
                    @printf("%10.8e   ", linf_error_prim[v])
                    print(io, " ", linf_error_prim[v])
                end
                println()
            end

            # L1 error
            if :l1_error_primitive in analysis_errors
                print(" L1 error prim.: ")
                for v in eachvariable(equations)
                    @printf("%10.8e   ", l1_error_prim[v])
                    print(io, " ", l1_error_prim[v])
                end
                println()
            end
        end
    end

    # additional integrals
    analyze_integrals(analysis_integrals, io, du, u, t, semi)

    # additional pointwise quantities
    analyze_pointwise(analysis_pointwise, du, u, t, semi, iter)

    return nothing
end

# Iterate over tuples of pointwise analysis quantities in a type-stable way using "lispy tuple programming".
function analyze_pointwise(analysis_quantities::NTuple{N, Any}, du, u, t,
                           semi, iter) where {N}

    # Extract the first pointwise analysis quantity and process it; keep the remaining to be processed later
    quantity = first(analysis_quantities)
    remaining_quantities = Base.tail(analysis_quantities)

    analyze(quantity, du, u, t, semi, iter)

    # Recursively call this method with the unprocessed pointwise analysis quantities
    analyze_pointwise(remaining_quantities, du, u, t, semi, iter)
    return nothing
end

# terminate the type-stable iteration over tuples
function analyze_pointwise(analysis_quantities::Tuple{}, du, u, t, semi, iter)
    nothing
end

# `analyze` function that passes also the iteration number `iter`along.
# Required for callbacks that handle the write-off of results to disk themselves,
# such as `AnalysisSurfacePointwise`.
function analyze(quantity, du, u, t, semi::AbstractSemidiscretization, iter)
    mesh, equations, solver, cache = mesh_equations_solver_cache(semi)
    analyze(quantity, du, u, t, mesh, equations, solver, cache, iter)
end
end # @muladd

# specialized implementations specific to some solvers
include("analysis_surface_integral.jl")
include("analysis_surface_pointwise.jl")

# This version of `analyze` is used for `AnalysisSurfacePointwise` such as `SurfacePressureCoefficient`.
# We need the iteration number `iter` to be passed in here 
# as for `AnalysisSurfacePointwise` the writing to disk is handled by the callback itself.
function analyze(quantity::AnalysisSurfacePointwise{Variable},
                 du, u, t,
                 semi::AbstractSemidiscretization,
                 iter) where {Variable}
    mesh, equations, solver, cache = mesh_equations_solver_cache(semi)
    # Call the `Variable`-specific `analyze` function
    analyze(quantity, du, u, t, mesh, equations, solver, cache, semi, iter)
end
# Special analyze for `SemidiscretizationHyperbolicParabolic` such that
# precomputed gradients are available. Required for `AnalysisSurfacePointwise` equipped 
# with `VariableViscous` such as `SurfaceFrictionCoefficient`.
# As for the inviscid version, we need to pass in the iteration number `iter` as 
# for `AnalysisSurfacePointwise` the writing to disk is handled by the callback itself.
function analyze(quantity::AnalysisSurfacePointwise{Variable},
                 du, u, t,
                 semi::SemidiscretizationHyperbolicParabolic,
                 iter) where {
                              Variable <:
                              VariableViscous}
    mesh, equations, solver, cache = mesh_equations_solver_cache(semi)
    equations_parabolic = semi.equations_parabolic
    cache_parabolic = semi.cache_parabolic
    # Call the `Variable`-specific `analyze` function
    analyze(quantity, du, u, t, mesh, equations, equations_parabolic, solver, cache, semi,
            cache_parabolic, iter)
end
