
<!-- README.md is generated from README.Rmd. Please edit that file -->
GREMLIN
=======

<!-- badges: start -->
<!-- badges: end -->
The goal of GREMLIN is to perform statistical analysis of multipartite networks through a block model approach.

Multipartite networks consist in the joint observation of several networks implying some common individuals. The individuals (or entities represented by nodes) at stake are partitioned into groups defined by their nature. In what follows, these groups will be referred to as .

Installation
------------

You can install the released version of GREMLIN [GitHub](https://github.com/) with:

``` r
#devtools::install_github("Demiperimetre/GREMLIN")
library(GREMLIN)
```

Mathematical Background
-----------------------

### A collection of networks

Assume that *Q* functional groups of individuals are at stake; Let *n*<sub>*q*</sub> be the number of individuals in the *q*-th functional group.

A multipartite network is a collection of networks: each network may be simple (relations inside a functional group) or bipartite (relations between individuals of two functional groups). We index the collection of networks by pairs of functional groups (*q*, *q*′).

The set *E* denotes the list of pairs of functional groups for which we observe an interaction network.

For any pair (*q*, *q*′) ∈ *E*, the interaction network is encoded in a matrix *X*<sup>*q**q*′</sup> : $X^{qq'}\_{ii'} 0 $ if there is an edge from unit *i* of functional group *q* to unit *i*′ of functional group *q*′, 0 otherwise. - If *q* ≠ *q*′, *X*<sup>*q**q*′</sup> is said to be an . - *X*<sup>*q**q*</sup> is an : it is symmetric if the relation inside the functional group *q* is non-oriented, non-symmetric otherwise.

### A block model

Assume that, each functional group *q* is divided into *K*<sub>*q*</sub> blocks (or equivalently clusters). ∀*q* ∈ {1, …, *Q*} and $ i {1,,n\_q}$, let *Z*<sub>*i*</sub><sup>*q*</sup> be the latent random variable such that *Z*<sub>*i*</sub><sup>*q*</sup> = *k* if individual *i* of functional group *q* belongs to cluster *k*. The random variables *Z*<sub>*i*</sub><sup>*q*</sup>'s are assumed to be independent and such that: ∀*k* ∈ {1, …, *K*<sub>*q*</sub>},∀*q* ∈ {1, …, *Q*},∀*i* ∈ {1, …, *n*<sub>*q*</sub>}:

with $\\sum\_{k=1}^{K\_q}\\pi^{q}\_k=1$, ∀*q* ∈ {1, …, *Q*}. Let $\\bZ = \\left(Z^{q}\_i\\right)\_{i\\in \\{1,\\ldots,n\_q\\}, q \\in \\{1,\\ldots,Q\\}}$ denote the set of latent variables.


