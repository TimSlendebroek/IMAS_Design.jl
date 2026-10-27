module IMASGraphRecipesExt

# Connection-graph plot of pf_active circuits (`plot(dd.pf_active.circuit[1])`): needs `using Plots, GraphRecipes`

using Plots
import GraphRecipes
import IMAS
import IMAS: IMASdd
# mirror the IMAS namespace (what the included code saw inside `module IMAS`): its own and imported names first,
# then the modules it brings in with `using` (names(...; usings=true) would need Julia 1.12)
for n in names(IMAS; all=true, imported=true)
    if Base.isidentifier(n) && n ∉ (:IMAS, :eval, :include) && isdefined(IMAS, n) && !isdefined(@__MODULE__, n)
        @eval import IMAS: $n
    end
end
using AdaptiveArrayPools, DSP, LinearAlgebra, Printf, RecipesBase
using Compat: @compat
using RecipesBase
import PlotUtils

include(joinpath("..", "src", "plot_circuit_graph.jl"))

end
