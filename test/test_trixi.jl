using Test: @test, @testset
using TrixiTest
using TrixiAero: examples_dir # important for right examples directory in tests

# Stripped down version from Trixi.jl

macro test_trixi_include(expr, args...)
    local add_to_additional_ignore_content = [
        # NOTE: These warnings arose from Julia 1.10 onwards
        r"WARNING: Method definition .* in module .* at .* overwritten .*.\n",
        # Some examples include an elixir with adaptive time stepping setting `tspan = (0.0, 0.0)`
        # to just get the definition of the problem and spatial discretization. In this case,
        # OrdinaryDiffEq.jl throws the following warning, which we can safely ignore in our tests:
        r"┌ Warning: Verbosity toggle: dt_epsilon \n│  Initial timestep too small \(near machine epsilon\), using default: dt = 0.0\n└ @ OrdinaryDiffEqCore ~/.julia/packages/OrdinaryDiffEqCore.*\n"
    ]
    # if `maxiters` is set in tests, it is usually set to a small number to
    # run only a few steps - ignore possible warnings coming from that
    if any(expr.args[1] == (:maxiters) for expr in args)
        push!(add_to_additional_ignore_content,
              r"┌ Warning: Verbosity toggle: max_iters \n│  Interrupted\. Larger maxiters is needed\..*\n└ @ Trixi .+\n",
              r"┌ Warning: Interrupted\. Larger maxiters is needed\..*\n└ @ Trixi .+\n")
    end
    args = append_to_kwargs(args, :additional_ignore_content,
                            add_to_additional_ignore_content)
    ex = quote
        @test_trixi_include_base($expr, $(args...))
    end
    return esc(ex)
end
