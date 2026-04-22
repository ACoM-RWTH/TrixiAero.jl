using Test
using Aqua
using TrixiAero

@testset "TrixiAero" begin
    @test true
end

@testset "Aqua" begin
    Aqua.test_all(TrixiAero)
end
