using IMAS
using Test

const SAMPLES = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample")

# Plots is a weak dependency (ext/IMASPlotsExt.jl). This test process may already have Plots loaded,
# so the "IMAS alone" behaviour is checked in a fresh Julia process using the same environment.
@testset "plotting is optional" begin
    code = """
        using IMAS
        loaded(name) = any(k.name == name for k in keys(Base.loaded_modules))
        println("plots=", loaded("Plots"), " ext=", Base.get_extension(IMAS, :IMASPlotsExt) !== nothing)
        dd = IMAS.json2imas(joinpath($(repr(SAMPLES)), "D3D_eq_ods.json"); show_warnings=false)
        try
            IMAS.plot(dd.equilibrium)
        catch e
            println("error=", sprint(showerror, e))
        end
        """
    out = read(`$(Base.julia_cmd()) --project=$(Base.active_project()) --startup-file=no -e $code`, String)
    @test occursin("plots=false ext=false", out)
    @test occursin("error=MethodError", out)
    @test occursin("run `using Plots` first", out)
end

ENV["GKSwstype"] = "nul" # headless GR
using Plots
using IMAS.HelpPlots

@testset "plot_recipes" begin
    @test Base.get_extension(IMAS, :IMASPlotsExt) !== nothing
    @test Base.get_extension(IMAS, :IMASGraphRecipesExt) === nothing # circuit graphs also need `using GraphRecipes`

    dd = IMAS.json2imas(joinpath(SAMPLES, "D3D_eq_ods.json"); show_warnings=false)
    @test plot(dd.equilibrium) isa Plots.Plot
    @test plot(dd.equilibrium.time_slice[]; cx=true) isa Plots.Plot

    # (not core_transport: this sample's transport ions have no label, which that recipe needs)
    dd = IMAS.json2imas(joinpath(SAMPLES, "omas_sample.json"); show_warnings=false)
    @testset "$name" for (name, ids) in (("core_profiles", dd.core_profiles), ("core_sources", dd.core_sources),
        ("pf_active", dd.pf_active), ("wall", dd.wall))
        @test plot(ids) isa Plots.Plot
    end

    help_plot(dd.equilibrium)
    help_plot!(dd.equilibrium)
end
