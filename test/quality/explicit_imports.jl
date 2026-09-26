#
# The module file's 37 `unused-import` findings from `fatou lint` were shown to be false rather
# than acted on, and this file is what stops a real one from coming back silently.
#
using ExplicitImports
using MetriplecticRelaxation
using Test

@testset "$(rpad("Explicit Imports", 80))" begin
    # This is the tool that answers what `fatou lint`'s `unused-import` rule only
    # approximates: the rule reads one file at a time and does not follow `include`, so in a
    # Julia package it flags the module file's load-bearing imports. ExplicitImports loads
    # the module and analyses real bindings.
    #
    # Staleness only. The other checks it offers are not asserted here: `field` is imported
    # from GeometricBrackets, which neither exports it nor declares it public, and that is a
    # known upstream defect recorded at the top of `src/MetriplecticRelaxation.jl`. Asserting
    # it would make this file red for something no change here can fix.
    @testset "$(rpad("no stale explicit import", 76))" begin
        @test ExplicitImports.check_no_stale_explicit_imports(MetriplecticRelaxation) ===
              nothing
    end
end
