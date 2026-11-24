# GREMLINS dev

* Add stop if convergence not reached and invitation to add more iterations

# GREMLINS 0.2.0 - first CRAN release

* Support for Multipartite SBM inference 
* Added a `NEWS.md` file to track changes to the package.

# GREMLINS 0.2.1 - fixing bugs

* Removed the plot function and encourage users to use the `sbm` package which depends on GREMLINS but is much more complete. 

# GREMLINS 0.2.2 - changing parallelization and initialisation

* Now uses future for the parallelization
* Vignettes converted to Quarto. They (and the README) now state that GREMLINS should be used through the package `sbm`,
  with the correspondence between the functions of the two packages.
* Reference to the published paper: Bar-Hen, Barbillon and Donnet (2022), Statistical Modelling <doi:10.1177/1471082X20963254> (also in `citation("GREMLINS")`).
* Initialisation: by default (`initBM = TRUE`), the starting classifications of `multipartiteBM()` and `multipartiteBMFixedModel()`
  are obtained by fitting a SBM / LBM on each network with `blockmodels` (all combinations are tried), instead of a hierarchical clustering.
  `initBM = FALSE` gives back the hierarchical clustering initialisation.
* Missing values (NA) are now handled in the initialisation (blockmodels and NA-safe distances in the hierarchical clustering).
* New argument `classifInit` in `multipartiteBM()` to provide an initial classification (same name as in `multipartiteBMFixedModel()`).
  In both functions, `classifInit` can be a list named with the functional groups, or an unnamed list in the order of `namesFG`; it is checked.
* Order of the functional groups (bugs appearing only when the internal order of the FGs is not a simple swap of `namesFG`, i.e. 3 FGs or more):
  `v_Kinit` in `multipartiteBM()` and `v_K` / `classifInit` in `multipartiteBMFixedModel()` were wrongly reordered.
* The results (`v_K`, `Z`, `list_pi`, `tau`) are given in the internal order of the FGs; `v_K` and `Z` are now always named accordingly
  (`v_K` was previously labelled with `namesFG`, which could be wrong).
* `multipartiteBMFixedModel()`: fixed a crash when some functional groups have only one block, and explicit error when no model converged.
