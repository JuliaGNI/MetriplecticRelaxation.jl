using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Explicit imports" include("quality/explicit_imports.jl")
    @safetestset "Torus geometry" include("torus.jl")
    @safetestset "Spectral solver" include("spectral.jl")
    @safetestset "Spline solver" include("spline.jl")
    @safetestset "Diagnostics" include("diagnostics.jl")
    @safetestset "Short runs" include("short_runs.jl")
    @safetestset "Projector rates" include("projector_rates.jl")
    @safetestset "Section 5 problems" include("euler_problems.jl")
    @safetestset "Euler square" include("euler_square.jl")
    @safetestset "Free state space" include("free_state_space.jl")
    @safetestset "Collision bracket" include("collision_bracket.jl")
    @safetestset "Gibbs entropy" include("gibbs_entropy.jl")
    @safetestset "Section 5 references" include("euler_references.jl")
    @safetestset "Section 5 short runs" include("euler_short_runs.jl")
    @safetestset "Grad-Shafranov problems" include("gradshafranov.jl")
end
