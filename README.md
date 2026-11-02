> [!NOTE]
> **About this fork.** IMAS_Design is part of [FUSE_Design](https://github.com/TimSlendebroek/FUSE_Design.jl), a public, stable version of the FUSE
> stack focused on **fusion power-plant design** rather than the time-dependent, digital-twin direction of upstream.
> Versions are `<upstream version at the split>+design.<n>`: the current release is **7.1.8+design.2**, split from
> upstream v7.1.8 ([CHANGELOG.md](CHANGELOG.md)). Plot recipes live in a Plots extension: `using Plots` loads them.
> To install the whole stack, see [Installing the design stack](https://github.com/TimSlendebroek/FUSE_Design.jl#installing-the-design-stack).

# IMAS.jl

## Physics and math routines
IMAS.jl incorporates a comprehensive set of mathematical, physics, and engineering routines. This eliminates the need for individual packages to re-implement these routines.

## Dynamic expressions
IMAS.jl implements the dynamic expressions that ensure consistency and provides an elegant solution to the mismatching-interfaces problem, where preceding models might not furnish all the derived data needed by subsequent models. See under `IMAS/src/expressions/`.

## Plotting
IMAS.jl makes use of Julia's Plots.jl and uses multiple dispatching mechanism to provide contextual (and composable) plotting capabilities throughout the data structure.

## Online documentation
For more details, see the [online documentation](https://projecttorreypines.github.io/IMAS.jl/dev).

![Docs](https://github.com/ProjectTorreyPines/IMAS.jl/actions/workflows/make_docs.yml/badge.svg)
