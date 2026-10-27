module IMAS

using Printf
using Compat:@compat
import OrderedCollections
import FastInterpolations as FI
using AdaptiveArrayPools
const document = OrderedCollections.OrderedDict()

macro import_all(mod)
    syms = names(eval(mod); all=true)
    syms = filter(n -> Base.isidentifier(n) && n ∉ (mod, :eval, :include, :document), syms)
    imports = [:(import $mod: $(s)) for s in syms]
    return Expr(:toplevel, imports...)
end

function __init__()
    # register IMAS expressions with IMASdd
    IMASdd.set_dynamic_expressions(dynamic_expressions)
    IMASdd.set_onetime_expressions(onetime_expressions)
    # plot(...) without Plots loaded is a MethodError on the RecipesBase stub: say how to fix it
    Base.Experimental.register_error_hint(_plots_missing_hint, MethodError)
end

function _plots_missing_hint(io::IO, exc::MethodError, argtypes, kwargs)
    f = exc.f === Core.kwcall && length(exc.args) >= 2 ? exc.args[2] : exc.f
    if (f === plot || f === plot!) && Base.get_extension(@__MODULE__, :IMASPlotsExt) === nothing
        print(io, "\nIMAS plot recipes load with Plots.jl, which IMAS does not load by default: run `using Plots` first.")
    end
end

#= ====== =#
#= IMASdd =#
#= ====== =#
import IMASdd
@import_all IMASdd
import IMASdd: @ddtime, @findall
import AbstractTrees: print_tree

#= ===== =#
#= UTILS =#
#= ===== =#
include("real.jl")
include("get_from.jl")
include("geometry.jl")
include("math.jl")
include("benchmark.jl")
include("signal.jl")

#= ======= =#
#= EXTRACT =#
#= ======= =#
include(joinpath("extract", "constraints.jl"))
include(joinpath("extract", "objectives.jl"))
include(joinpath("extract", "extract.jl"))

#= ======= =#
#= CONTROL =#
#= ======= =#
include(joinpath("control", "control.jl"))
include(joinpath("control", "fxp.jl"))

#= ======= =#
#= PHYSICS =#
#= ======= =#
include("physics.jl")
include("experiments.jl")

#= =========== =#
#= EXPRESSIONS =#
#= =========== =#
include(joinpath("expressions", "onetime.jl"))
include(joinpath("expressions", "dynamic.jl"))

#= ======== =#
#= PLOTTING =#
#= ======== =#
# Recipes and plot helpers (plot.jl) live in ext/IMASPlotsExt.jl and load with `using Plots`.
# Debug and do_plot code in src/ calls the RecipesBase plot/plot! stubs, which Plots implements.
import RecipesBase: plot, plot!

#= ====== =#
#= EXPORT =#
#= ====== =#
export @ddtime, @findall, help

import HelpPlots # keeps IMAS.HelpPlots available now that plot.jl lives in ext/IMASPlotsExt.jl
import HelpPlots: help_plot, help_plot!
export help_plot, help_plot!
export print_tree

end # module
