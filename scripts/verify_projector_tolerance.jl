#!/usr/bin/env julia
#
# The control behind `projector_run.jl`'s 𝔠_η membership bound.
#
#     julia --project=scripts scripts/verify_projector_tolerance.jl
#
# A tolerance is only a bound if something plausible fails it. Section 6 of `projector_run.jl`
# asserts that the relaxed state is a member of `eq:u-eta_Euler_periodic` at a relative residual
# below `3e-4`, and that constant is a stated multiple -- 9.0x -- of A3's measured 3.33e-05. This
# script is the evidence that the multiple is not decorative: it takes A3's own run to half its
# final time and requires the under-relaxed state to be REJECTED.
#
# WHY A COARSER MESH IS NOT THE CONTROL, which is the finding this script exists to record.
# The obvious degradation is to re-run on a coarser space and see whether the row notices. It does
# not: measured, the residual is 3.33168e-05 at 24 cells against 3.33164e-05 at 64 -- identical to
# five digits. The residual is not set by the spatial resolution. It is set by how far the
# relaxation has got at T, so the degradation that tests it is a SHORTER RUN. A mesh degradation
# here returns a null result, which reads as a tolerance that cannot be validated when what is
# actually wrong is the choice of knob.
#
# Measured on A3's own space, fitting at each fraction of T:
#
#     t/T      0.02       0.05       0.10       0.25       0.50       1.00
#     rel      8.37e-01   6.55e-01   3.72e-01   6.45e-02   4.96e-03   3.33e-05
#
# The run stops at T/2 because that is where the discriminating statement is: five of those six
# points lie at or below it, and the sixth is what `run_a3.jl` already asserts. Half the steps,
# all of the content -- about two minutes.
#
# `Knowledge/Metriplectic Relaxation/A control that cannot fail proves nothing.md`, mechanism 8,
# is the catalogue entry, and this is the case that added the limit on its own remedy.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION4_RUNS, SplineTorus, Diagnostics, Trace,
                              spline_state, spline_rhs, spline_step, record!,
                              best_fit_euler, l2norm
using Printf

include(joinpath(@__DIR__, "check.jl"))
using .Checks: header, check, summary

println("verify_projector_tolerance.jl  --  the 𝔠_η membership bound rejects an under-relaxed state")

# Must match the bound asserted in `projector_run.jl`'s section 6. Two literals rather than one
# shared constant because that file is a module the drivers include and this script does not; if
# one moves and the other does not, this row is what notices — the margins below are wide enough
# that only a real change in the bound reaches them.
const MEMBERSHIP_TOL = 3e-4

# A3's own space and A3's own step, at the resolution `run_a3.jl` uses.
const spec = SECTION4_RUNS["a3"]
const CELLS = 64
const DEGREE = 3

const t = SplineTorus(CELLS, DEGREE)
const d = Diagnostics(t, spec)
const rhs = spline_rhs(t, spec)

ω̂, _ = spline_state(t, spec)
const tr = Trace(ω̂)
record!(tr, d, 0.0, ω̂)
const H₀ = tr.H[1]

const nsteps = round(Int, spec.T / spec.Δt)
const FRACTIONS = (0.02, 0.05, 0.10, 0.25, 0.50)

@printf("    A3, %d cubic cells, Δt = %.0e, T = %.1f — stepping to T/2 (%d of %d steps)\n",
    CELLS, spec.Δt, spec.T, nsteps ÷ 2, nsteps)
flush(stdout)

const measured = Float64[]
let done = 0
    for f in FRACTIONS
        target = round(Int, f * nsteps)
        for _ in (done + 1):target
            global ω̂ = spline_step(rhs, ω̂, spec.Δt)
        end
        done = target
        (_, resid, _) = best_fit_euler(d, ω̂, H₀)
        push!(measured, resid / l2norm(t, ω̂))
    end
end

# =============================================================================================
header("1. the residual falls monotonically as the relaxation proceeds")

# Without this the section below would pass for a residual that was simply noisy, and the claim
# is that the quantity tracks the relaxation rather than that it happens to be large early on.
for k in 2:length(FRACTIONS)
    check(@sprintf("t/T = %.2f is below t/T = %.2f", FRACTIONS[k], FRACTIONS[k - 1]),
        measured[k] < measured[k - 1],
        @sprintf("%.5e  ->  %.5e", measured[k - 1], measured[k]))
end

# =============================================================================================
header("2. THE CONTROL: an under-relaxed state is REJECTED by the bound")

let rel = measured[end]
    check("at half its relaxation, A3 FAILS the 𝔠_η membership bound",
        rel > MEMBERSHIP_TOL,
        @sprintf("rel = %.5e at t/T = 0.50   tol %.0e   exceeds it by %.0fx",
            rel, MEMBERSHIP_TOL, rel / MEMBERSHIP_TOL))
end

# And by a margin that is not marginal: a bound this state only just failed would be one an
# ordinary fluctuation could carry either way.
let rel = measured[end]
    check("and by more than a factor of ten, so the rejection is not marginal",
        rel > 10 * MEMBERSHIP_TOL,
        @sprintf("rel / tol = %.1f", rel / MEMBERSHIP_TOL))
end

# =============================================================================================
header("3. and the state at t/T = 0.50 is genuinely relaxing, not stalled")

# The complement of section 2. If the run had not moved at all, "an under-relaxed state fails the
# bound" would be true of the initial condition and would say nothing about the trajectory.
let first = measured[1], last = measured[end]
    check("the residual has fallen by more than two orders between t/T = 0.02 and 0.50",
        first / last > 100,
        @sprintf("%.5e  ->  %.5e   a factor %.0f", first, last, first / last))
end

summary("verify_projector_tolerance.jl")
