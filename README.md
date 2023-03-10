
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GREMLINS


The goal of GREMLINS is to perform statistical analysis of multipartite
networks through a block model approach.

Multipartite networks consist in the joint observation of several
networks implying some common individuals. The individuals (or entities
represented by nodes) at stake are partitioned into groups defined by
their nature. In what follows, these groups will be referred to as .

## Installation

You can install the released version of GREMLINS
[GitHub](https://github.com/) with:

``` r
#devtools::install_github("GrossSBM/GREMLINS")
library(GREMLINS)
```

## Mathematical Background

### A collection of networks

Assume that \(Q\) functional groups of individuals are at stake; Let
\(n_q\) be the number of individuals in the \(q\)-th functional group.

A multipartite network is a collection of networks: each network may be
simple (relations inside a functional group) or bipartite (relations
between individuals of two functional groups). We index the collection
of networks by pairs of functional groups \((q,q')\).

The set \(E\) denotes the list of pairs of functional groups for which
we observe an interaction network.

For any pair \((q,q') \in E\), the interaction network is encoded in a
matrix \(X^{qq'}\) : \(X^{qq’}_{ii’}= 0 \) if there is an edge from unit
\(i\) of functional group \(q\) to unit \(i'\) of functional group
\(q'\), \(0\) otherwise.

- If \(q \neq q'\), \(X^{qq'}\) is said to be an incidence matrix.   

- \(X^{qq}\) is an adjacency matrix: it is symmetric if the relation inside the
functional group \(q\) is non-oriented, non-symmetric otherwise.

### A block model

Assume that, each functional group \(q\) is divided into \(K_q\) blocks
(or equivalently clusters). \(\forall q\in \{1,\ldots,Q\}\) and \(i= 
{1,,n\_q}\), let \(Z^{q}_i\) be the latent random variable such that
\(Z^ q_i =k\) if individual \(i\) of functional group \(q\) belongs to
cluster \(k\). The random variables \(Z^{q}_i\)’s are assumed to be
independent and such that:
\(\forall k \in \{1,\ldots,K_q\}, \forall q \in \{1,\ldots,Q\}, \forall i \in \{1,\ldots,n_q\}\):

$$P(Z^q_i=k) = \pi^{q}_{k} $$

with \(\sum_{k=1}^{K_q}\pi^{q}_{k}\), \(\forall q \in \{1,\ldots,Q\}\).

Then the connections  are distributed conditionally to the blocks as follows: 

$$ X^{qq'}_{ii'} | Z_i^q = k , Z_{i'}^{q'} = k' \sim \mathcal{F}_{qq'}(\theta^{qq'}_{kk'})$$


