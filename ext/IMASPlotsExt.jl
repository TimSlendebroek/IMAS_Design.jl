module IMASPlotsExt

# All IMAS plot recipes and plot helpers. Loaded automatically after `using Plots`.
# Worker processes that never load Plots skip it (and Plots with its dependencies) entirely.

using Plots
import IMAS
import IMAS: IMASdd
# mirror the IMAS namespace: everything visible inside `module IMAS` (its own names, plus what it imports or uses)
for n in names(IMAS; all=true, imported=true, usings=true)
    if Base.isidentifier(n) && n ∉ (:IMAS, :eval, :include) && isdefined(IMAS, n) && !isdefined(@__MODULE__, n)
        @eval import IMAS: $n
    end
end
import IMASdd: @ddtime, @findall
using Printf

include(joinpath("..", "src", "plot.jl"))

end
