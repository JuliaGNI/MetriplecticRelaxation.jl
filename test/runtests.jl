using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Explicit imports" include("quality/explicit_imports.jl")
    @safetestset "Torus geometry" include("torus.jl")
    @safetestset "Spectral solver" include("spectral.jl")
    @safetestset "Spline solver" include("spline.jl")
    @safetestset "Diagnostics" include("diagnostics.jl")
    @safetestset "Short runs" include("integration/short_runs.jl")
    @safetestset "Projector rates" include("integration/projector_rates.jl")
    @safetestset "Section 5 problems" include("integration/euler_problems.jl")
    @safetestset "Euler square" include("euler.jl")
    @safetestset "Free state space" include("integration/free_state_space.jl")
    @safetestset "Collision bracket" include("integration/collision_bracket.jl")
    @safetestset "Gibbs entropy" include("integration/gibbs_entropy.jl")
    @safetestset "Section 5 references" include("integration/euler_references.jl")
    @safetestset "Section 5 short runs" include("integration/euler_short_runs.jl")
    @safetestset "Grad-Shafranov problems" include("gradshafranov.jl")
end
