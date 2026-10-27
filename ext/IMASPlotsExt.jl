module IMASPlotsExt

# All IMAS plot recipes and plot helpers. Loaded automatically after `using Plots`.
# Worker processes that never load Plots skip it (and Plots with its dependencies) entirely.

using Plots
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
import IMASdd: @ddtime, @findall
using Printf

include(joinpath("..", "src", "plot.jl"))

end
