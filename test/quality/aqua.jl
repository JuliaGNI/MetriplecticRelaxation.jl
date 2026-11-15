#
# The package-quality guards. The module has no method ambiguities, and a new dependency version
# can add one without any change in this repository.
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

    # Every check is green, so a failure here is a regression and not a backlog.
    @testset "$(rpad("Aqua", 76))" begin
        Aqua.test_all(MetriplecticRelaxation)
    end
end
