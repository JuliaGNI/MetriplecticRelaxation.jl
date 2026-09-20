#
# The quality guards. Two results were settled by hand and neither had anything holding it
# there: `detect_ambiguities` returned 2 and now returns 0, and the module file's 37
# `unused-import` findings were shown to be false rather than acted on. Both can come back
# silently, and this file is what stops that.
#
# The ambiguity guard carries NO dependency, so it runs under the documented local command
#
#     julia --project=. --startup-file=no test/runtests.jl
#
# as well as under CI's `Pkg.test()`. That matters: it is the guard for the result most likely
# to regress, because a new dependency version can reintroduce an ambiguity without this
# repository changing at all.
#
# Aqua and ExplicitImports are test dependencies, declared in the root `Project.toml`'s
# `[extras]` and `[targets]`, and they reach this file by two different routes. On CI,
# `Pkg.test()` resolves them into the test environment. Locally they come from the shared
# `@v#.#` environment, which is on Julia's default load path beside `--project=.` — neither is
# in this repository's `Manifest.toml`, which is hand-seeded and gitignored and which no `Pkg`
# call may touch.
#
# When neither route supplies them the testset marks itself broken and warns. It does not pass.
# A guard that reports success having checked nothing is worse than no guard, because it also
# removes the reason to look.
#

const QUALITY_DEPS = all(p -> Base.find_package(p) !== nothing, ("Aqua", "ExplicitImports"))

if QUALITY_DEPS
    import Aqua
    import ExplicitImports
end

@testset "$(rpad("Quality Guards", 80))" begin

    # `recursive = false` is the whole module and not its dependencies' internals. Both findings
    # this guards sat on `*(::PoissonMap, ::AbstractVector)`: `PoissonMap <: AbstractMatrix`, so
    # that signature met ArrayLayouts' `*(::AbstractMatrix, ::LayoutVector)` and FillArrays'
    # `*(::AbstractMatrix{T}, ::AbstractZeros{T,1})` with neither side the more specific.
    @testset "$(rpad("the module carries no method ambiguities", 76))" begin
        ambiguities = Test.detect_ambiguities(MetriplecticRelaxation; recursive = false)
        @test isempty(ambiguities)
        isempty(ambiguities) ||
            foreach(a -> println(stderr, "  ambiguous: ", a), ambiguities)
    end

    if !QUALITY_DEPS
        @warn """Aqua and ExplicitImports are not in this environment, so the import and \
                 package-quality guards did NOT run. See the header of test/runquality.jl."""
        @testset "$(rpad("the import and package-quality guards ran", 76))" begin
            @test_broken QUALITY_DEPS
        end
    else
        # This is the tool that answers what `fatou lint`'s `unused-import` rule only
        # approximates: the rule reads one file at a time and does not follow `include`, so in a
        # Julia package it flags the module file's load-bearing imports. ExplicitImports loads
        # the module and analyses real bindings.
        #
        # Staleness only. The other checks it offers are not asserted here: `field` is imported
        # from PoissonBrackets, which neither exports it nor declares it public, and that is a
        # known upstream defect recorded at the top of `src/MetriplecticRelaxation.jl`. Asserting
        # it would make this file red for something no change here can fix.
        @testset "$(rpad("no stale explicit import", 76))" begin
            @test ExplicitImports.check_no_stale_explicit_imports(MetriplecticRelaxation) ===
                  nothing
        end

        # All eight checks were measured green before this was wired in, so a failure here is a
        # regression and not a backlog.
        @testset "$(rpad("Aqua", 76))" begin
            Aqua.test_all(MetriplecticRelaxation)
        end
    end
end
