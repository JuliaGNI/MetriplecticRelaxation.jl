#
# The PASS/FAIL harness the verification scripts report through.
#
# The output format is copied deliberately from `common.py` in the manuscripts' `Scripts/`
# folders, down to the two leading spaces and the three spaces before the detail, so that
# a converted script can be checked against the Python it replaces with a plain `diff`.
# Once the Python is retired that constraint is free to go, but until then it is what makes
# the comparison mechanical rather than a matter of reading two tables side by side.
#
#     include(joinpath(@__DIR__, "check.jl"))
#     using .Checks: header, check, summary
#
# Import the three names explicitly rather than with a bare `using`: `summary` also exists in
# `Base`, and a bare `using` leaves the two ambiguous rather than preferring ours.
#
#     header("1. the bracket is antisymmetric")
#     check("J + J' = 0", norm(J + J') < 1e-14, @sprintf("%.2e", norm(J + J')))
#     summary("verify_something.jl")
#
# `summary` exits 1 if any check failed, so `run_all.jl` and CI see a nonzero status.
#
# `check_exact` and `check_refined` are the same thing for a claim about a continuum identity
# rather than a matrix: one settles it against an absolute threshold, the other against the
# factor by which the residual drops when the resolution is doubled.  They came in with the
# four-bracket manuscript, which has no Python original, so the format note above does not
# constrain them -- but they print through `check`, so one line per claim still holds.
#
# The functions here know nothing about grids.  Which resolutions a refinement uses is the
# calling manuscript's business and lives in `torustools.jl`.

module Checks

using Printf

export header, check, check_exact, check_refined, summary, fmt, relerr, normerr,
       failures, reset_failures!

const _failures = String[]

"""
    header(text)

A blank line, `text`, and a rule of matching length.  Sections of a verification script are
numbered in the Python; keep the numbering so the two outputs line up.
"""
function header(text)
    println()
    println(text)
    println("-"^length(text))
    return nothing
end

"""
    check(label, condition, detail = "")

Record and print one claim.  Returns the boolean, so a check can gate what follows it.

`detail` is the measured quantity behind the verdict -- a residual, a rate, a rank.  A bare
PASS says only that someone chose the tolerance well; the number says what was actually
seen, and it is the number that gets compared against the Python.
"""
function check(label, condition, detail = "")
    ok = Bool(condition)
    ok || push!(_failures, label)
    mark = ok ? "PASS" : "FAIL"
    println("  [$mark] $label" * (isempty(detail) ? "" : "   $detail"))
    return ok
end

"""
    check_exact(label, residual; atol = 1e-11)

One claim whose identity is exact algebra plus integration by parts, so that on a spectral
grid the residual must be at roundoff.  Passes if `abs(residual) <= atol`.
"""
function check_exact(label, residual; atol = 1e-11)
    check(label, abs(residual) <= atol, @sprintf("%.2e   (atol %.0e)", abs(residual), atol))
end

"""
    check_refined(label, coarse, fine; atol = 1e-11, minrate = 40)

One claim whose identity involves fields that are not band-limited -- a quotient, a power, a
logarithm -- so that the residual is discretisation error and no absolute threshold applies.
Passes if the residual on the finer grid is already at roundoff, or if it fell by at least
`minrate` when the resolution was doubled.

`minrate` defaults to 40 against the 256 an 8th-order scheme predicts.  The margin absorbs
the constant in the leading error term, which is not close to one for the quotient
identities, without admitting anything that converges at 5th order or worse.

The detail records both residuals and the factor between them, because the factor *is* the
evidence here: a bare PASS would say only that someone chose `minrate` well, whereas 251x
against a predicted 256 says the identity holds and the scheme is behaving.
"""
function check_refined(label, coarse, fine; atol = 1e-11, minrate = 40.0)
    c, f = abs(coarse), abs(fine)
    rate = c / max(f, 1e-300)
    return check(label, f <= atol || rate >= minrate,
        @sprintf("%.2e -> %.2e   %.0fx", c, f, rate))
end

"""
    summary(title)

Print the tally and exit 1 if anything failed.  Call it once, at the end of a script.
"""
function summary(title)
    println()
    if !isempty(_failures)
        println("$title: $(length(_failures)) FAILURE(S): " * join(_failures, ", "))
        exit(1)
    end
    println("$title: all checks passed.")
    return nothing
end

"""
    fmt(x)

Render a value for the `detail` field of [`check`](@ref).

Rationals print as `p/q`, and as `p` when the denominator is one, rather than in Julia's
`p//q` form. That is not a cosmetic whim: it is the form the Python prototypes print, and
keeping it means a `diff` of the two outputs shows only the numbers that genuinely differ.
"""
function fmt(x::Rational)
    isone(denominator(x)) ? string(numerator(x)) :
    string(numerator(x), "/", denominator(x))
end
fmt(x::AbstractFloat) = string(x)
fmt(x) = string(x)

"""
    relerr(a, b)

The relative error of the claim `a == b`, measured against the larger of the two.
"""
relerr(a, b) = abs(a - b) / max(abs(a), abs(b), 1e-300)

"""
    normerr(a, b, scale)

The error of `a == b` measured against an externally supplied `scale`.

For some test fields both sides of an identity vanish, and then [`relerr`](@ref) divides one
small number by another and reports noise.  Pass the integral of the absolute value of the
integrand as `scale`: it is the size of the terms that were actually formed and cancelled,
which is the quantity the residual has to be small compared with.
"""
normerr(a, b, scale) = abs(a - b) / max(abs(scale), 1e-300)

"The labels of the checks that have failed so far."
failures() = copy(_failures)

"Forget the recorded failures.  Only `run_all.jl` needs this, between scripts."
reset_failures!() = (empty!(_failures); nothing)

end # module
