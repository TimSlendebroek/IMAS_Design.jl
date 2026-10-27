module IMASGraphRecipesExt

# Connection-graph plot of pf_active circuits (`plot(dd.pf_active.circuit[1])`): needs `using Plots, GraphRecipes`

using Plots
import GraphRecipes
import IMAS
import IMAS: IMASdd
# mirror the IMAS namespace: everything visible inside `module IMAS` (its own names, plus what it imports or uses)
for n in names(IMAS; all=true, imported=true, usings=true)
    if Base.isidentifier(n) && n ∉ (:IMAS, :eval, :include) && isdefined(IMAS, n) && !isdefined(@__MODULE__, n)
        @eval import IMAS: $n
    end
end
using RecipesBase
import PlotUtils

include(joinpath("..", "src", "plot_circuit_graph.jl"))

end
