
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

### Example of a mutualistic ecological network

``` r
data = read.csv(file = 'data/ecologicalMPnetwork.csv',header = TRUE,skip = 1)
plantnames = data[,2]
data2 = read.csv(file = 'data/ecologicalMPnetwork.csv',header = FALSE)
type_animals = data2[1,]
type_animals = as.factor(as.character(type_animals[-c(1,2)]))
indflovis = which(data2[1,] == "Floral visitor")
indant = which(data2[1,] == "Ant")
indseed = which(data[1,] == "Seed dispersal")

Adj_flovis = data[,indflovis]
rownames(Adj_flovis) = plantnames
Adj_flovis = as.matrix(Adj_flovis)
Adj_ant = data[,indant]
rownames(Adj_ant) = plantnames
Adj_ant = as.matrix(Adj_ant)
Adj_seed = data[,indseed]
rownames(Adj_seed) = plantnames
Adj_seed = as.matrix(Adj_seed)
```

In this example, the functional groups are plants, floral visitors, ants and birds.

GREMLIN requires the global network be encoded in separate matrices for each network. 3 incidence matrices. We then format the data to be able to use our R package GREMLIN i.e. we transform the matrices into an list containing the matrix, its type (incidence, adjacency symmetric, adjacency non symmetric, the functional group in row, the functional group in col).

``` r
PlantFlovis = defineNetwork(Adj_flovis,"inc","plants","flovis")
PlantAnt = defineNetwork(Adj_ant,"inc","plants","ants")
PlantBird = defineNetwork(Adj_seed,"inc","plants","birds")
```

If you want to keep the names of the species, they should be used as rownames and colnames in the matrices.

``` r
PlantFlovis$mat[1:2,1:2]
```

``` r
Q <- 2
v_NQ = c(30,50)
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(GREMLIN)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub!
