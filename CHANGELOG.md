# IMAS_Design changelog

IMAS_Design is the [IMAS.jl](https://github.com/ProjectTorreyPines/IMAS.jl) of the FUSE design stack (see FUSE_Design).
Versions are `<upstream IMAS version at the split>+design.<n>`, e.g. `7.1.8+design.2`; they satisfy the same compat
bounds as the upstream version and never clash with upstream release numbers. Releases are git tags `v<version>`.

## 7.1.8+design.2 (2026-11-02, branch `lean-deps`)

`using IMAS` takes ~3 s and ~0.7 GB instead of ~6 s and ~0.9 GB.

- Plot recipes and helpers (`src/plot.jl`) live in `ext/IMASPlotsExt.jl` and load with `using Plots`; the circuit graph
  recipe in `ext/IMASGraphRecipesExt.jl`. Plots and GraphRecipes are weak dependencies, Graphs is dropped.
- `plot(...)` without Plots is a MethodError on the RecipesBase stub; an error hint says to run `using Plots`.
  `IMAS.HelpPlots` stays available.
- Diagnostic plots in error paths (negative ion density in `sivukhin_fraction`, wall intersections in flux surfaces)
  only run when Plots is loaded, so without Plots the real error is raised. `IMAS.plots_loaded()` tells whether
  Plots is loaded.
- Tests: fresh-process check that IMAS works without Plots, recipes after `using Plots`.
- CI: FuseRegistry from the TimSlendebroek fork with ProjectTorreyPines made unreachable (no read token needed).

## 7.1.8+design.1 (2026-10-19, tag on 2580e8ad)

The IMAS of the October 2026 NERSC CMOOP studies: upstream v7.1.8 (2026-07-08) plus

- `area()` of open polygons fixed: PF/OH coil cost was 8-15x too high.
- L-H power threshold without the factor ~2 ∇B-drift difference.
- `alpha_t` in extract, physics and MHD, and as a constraint.
- `pedestal_poloidal_beta(eqt, pped)` for the EPED analytic width law.
- Ohmic power density clamped to non-negative values.
- IMASdd 8.4 compatibility (`li_1`, `q_min` sub-IDS); `ion.rotation_frequency_tor` no longer a dynamic expression.
