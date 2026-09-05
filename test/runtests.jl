using MetriplecticRelaxation
using Random
using Test

# The runs here excite fields with random degrees of freedom, which is deliberate: a
# conservation or degeneracy identity tested on a single smooth mode reports round-off and hides
# a structural defect entirely. The seed is fixed so that a failure is reproducible.
Random.seed!(0x5c1e9a3b)

@testset "$(rpad("Scaffolding Tests",80))" begin
    # Placeholder until the first model lands: asserts only that the package LOADS, which is
    # what a suite with no other content can honestly claim.
    @test isdefined(Main, :MetriplecticRelaxation)
end
