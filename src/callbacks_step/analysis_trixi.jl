# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# This is the convenience constructor that gets called from the elixirs
function AnalysisCallback(semi::AbstractSemidiscretization; kwargs...)
    mesh, equations, solver, cache = mesh_equations_solver_cache(semi)
    return AnalysisCallback(mesh, equations, solver, cache; kwargs...)
end

# This method gets called from OrdinaryDiffEq's `solve(...)`
function initialize!(cb::DiscreteCallback{Condition, Affect!}, u_ode, t,
                     integrator) where {Condition, Affect! <: AnalysisCallback}
    semi = integrator.p
    du_ode = first(get_tmp_cache(integrator))
    return initialize!(cb, u_ode, du_ode, t, integrator, semi)
end

# This is the actual initialization method
# Note: we have this indirection to allow initializing a callback from the AnalysisCallbackCoupled
function initialize!(cb::DiscreteCallback{Condition, Affect!}, u_ode, du_ode, t,
                     integrator, semi) where {Condition, Affect! <: AnalysisCallback}
    initial_state_integrals = integrate(u_ode, semi)
    _, equations, _, _ = mesh_equations_solver_cache(semi)

    analysis_callback = cb.affect!
    analysis_callback.initial_state_integrals = initial_state_integrals
    @unpack save_analysis, output_directory, analysis_filename, analysis_errors, analysis_integrals = analysis_callback

    if save_analysis && mpi_isroot()
        mkpath(output_directory)

        # write header of output file
        open(joinpath(output_directory, analysis_filename), "w") do io
            print(io, "#timestep ")
            print(io, "time ")
            print(io, "dt ")
            if :l2_error in analysis_errors
                for v in varnames(cons2cons, equations)
                    print(io, "l2_" * v * " ")
                end
            end
            if :linf_error in analysis_errors
                for v in varnames(cons2cons, equations)
                    print(io, "linf_" * v * " ")
                end
            end
            if :l1_error in analysis_errors
                for v in varnames(cons2cons, equations)
                    print(io, "l1_" * v * " ")
                end
            end
            if :conservation_error in analysis_errors
                for v in varnames(cons2cons, equations)
                    print(io, "cons_" * v * " ")
                end
            end
            if :residual in analysis_errors
                for v in varnames(cons2cons, equations)
                    print(io, "res_" * v * " ")
                end
            end
            if :l2_error_primitive in analysis_errors
                for v in varnames(cons2prim, equations)
                    print(io, "l2_" * v * " ")
                end
            end
            if :linf_error_primitive in analysis_errors
                for v in varnames(cons2prim, equations)
                    print(io, "linf_" * v * " ")
                end
            end
            if :l1_error_primitive in analysis_errors
                for v in varnames(cons2prim, equations)
                    print(io, "l2_" * v * " ")
                end
            end

            for quantity in analysis_integrals
                print(io, pretty_form_ascii(quantity), " ")
            end
            # Pointwise quantities are not saved in `analysis_filename`,
            # i.e., `analysis.dat` but handle their own output.

            println(io)
            return nothing
        end
    end

    # Record current time using a high-resolution clock
    analysis_callback.start_time = time_ns()

    # Record current time for performance index computation
    analysis_callback.start_time_last_analysis = time_ns()

    # Record current number of `rhs!` calls for performance index computation
    analysis_callback.ncalls_rhs_last_analysis = ncalls(semi.performance_counter)

    # Record total time spent in garbage collection so far using a high-resolution clock
    # Note: For details see the actual callback function below
    analysis_callback.start_gc_time = Base.gc_time_ns()

    analysis_callback(u_ode, du_ode, integrator, semi)
    return nothing
end

# This method gets called from OrdinaryDiffEq's `solve(...)`
function (analysis_callback::AnalysisCallback)(integrator)
    semi = integrator.p
    du_ode = first(get_tmp_cache(integrator))
    u_ode = integrator.u
    return analysis_callback(u_ode, du_ode, integrator, semi)
end
end # @muladd
