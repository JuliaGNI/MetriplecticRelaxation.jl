#
# The package-quality guards. `detect_ambiguities` returned 2 and now returns 0, and nothing but
# this file holds it there: a new dependency version can reintroduce an ambiguity without this
# repository changing at all.
#
using Aqua
using MetriplecticRelaxation
using Test

@testset "$(rpad("Quality Guards", 80))" begin

    # `recursive` controls whether SUBMODULES are descended into, not dependencies. This module
    # has none, so both settings give the same answer and `false` is the cheaper one. Both
    # findings this guards sat on `*(::PoissonMap, ::AbstractVector)`: `PoissonMap <:
    # AbstractMatrix`, so that signature met ArrayLayouts' `*(::AbstractMatrix, ::LayoutVector)`
    # and FillArrays' `*(::AbstractMatrix{T}, ::AbstractZeros{T,1})` with neither side the more
    # specific.
    @testset "$(rpad("the module carries no method ambiguities", 76))" begin
        ambiguities = Test.detect_ambiguities(MetriplecticRelaxation; recursive = false)
        @test isempty(ambiguities)
        isempty(ambiguities) ||
            foreach(a -> println(stderr, "  ambiguous: ", a), ambiguities)
    end

    # All eight checks were measured green before this was wired in, so a failure here is a
    # regression and not a backlog.
    #
    # `persistent_tasks = false`, and it is the one check deliberately off. Aqua 0.8.16 runs
    # it by generating a temporary project and `Pkg.develop`-ing this package into it
    # (`persistent_tasks.jl:93-95`), then precompiling in a subprocess (`:114`). That
    # project has no manifest, so the resolve has to satisfy this package's `rev = "main"`
    # GitHub `[sources]` — which makes every run of the suite need the network and track a
    # moving branch. Against that it can find nothing here: this package has no `__init__`, no
    # `@async`, no `Threads.@spawn` and no `Timer`, so there is no task for it to catch.
    @testset "$(rpad("Aqua", 76))" begin
        Aqua.test_all(MetriplecticRelaxation; persistent_tasks = false)
    end
end
