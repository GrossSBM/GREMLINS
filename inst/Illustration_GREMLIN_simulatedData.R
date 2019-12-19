## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)


## ----packages, eval=TRUE-------------------------------------------------
#devtools::install_github("Demiperimetre/GREMLIN")
library(GREMLIN)


## ----param FG------------------------------------------------------------
set.seed(302718)
n_FG <- 3 #number of functional groups
v_NQ = c(100,50,60)
v_K <- c(4,3,2) #number of clusters in each functional group
list_pi = lapply(1:n_FG,function(q){v = rgamma(v_K[q],1,1); return(v/sum(v))})
namesFG <- c('A','B','C')


## ----param net-----------------------------------------------------------
E  <-  rbind(c(1,2),c(2,3),c(1,1),c(3,3))
typeInter <- c( "inc",'inc',"diradj", "adj")
v_distrib <- c('gaussian','bernoulli','poisson','bernoulli')
list_theta <- list()
list_theta[[1]] <- list()
list_theta[[1]]$mean  <- matrix(round(rnorm(v_K[E[1,1]] * v_K[E[1,2]],7.5,4 ),1),nrow = v_K[E[1,1]], ncol = v_K[E[1,2]] )
list_theta[[1]]$var  = matrix(round(rgamma(v_K[E[1,1]] * v_K[E[1,2]],7.5,4 ),1),nrow = v_K[E[1,1]], ncol = v_K[E[1,2]] )
for (q in 2:4) {
  list_theta[[q]] <- switch(v_distrib[q],
      bernoulli = matrix(rbeta(v_K[E[q,1]] * v_K[E[q,2]],2,2 ),nrow = v_K[E[q,1]], ncol = v_K[E[q,2]]),
      poisson = matrix(rgamma(v_K[E[q,1]] * v_K[E[q,2]],6,2 ),nrow = v_K[E[q,1]], ncol = v_K[E[q,2]])
  )
}
list_theta[[4]] =0.5 * (list_theta[[4]] + t(list_theta[[4]]))



## ----simul, eval = TRUE--------------------------------------------------
dataSim <- rMBM(v_NQ ,E , typeInter, v_distrib, list_pi,
                list_theta, seed = NULL,
                namesFG =namesFG,keepClassif  = TRUE)
list_Net <- dataSim$list_Net
length(list_Net)

save(dataSim,file='vignettes/dataSim.Rda')

## ----MBM, echo = FALSE, eval = TRUE--------------------------------------
if (!file.exists('vignettes/resMBM_Simu.Rda')) {
  res_MBMsimu <- multipartiteBM(list_Net,
                      v_distrib = v_distrib,
                      v_Kmin = 1,
                      v_Kmax = 10,
                      v_Kinit = NULL,
                      verbose = TRUE,
                      save=FALSE, initBM = TRUE)
  save(RES_MBM,file="vignettes/resMBM_Simu.Rda")
} else {load("vignettes/resMBM_Simu.Rda")}


## ----MBM eval false, echo = TRUE, eval = FALSE---------------------------
#>   res_MBMsimu <- multipartiteBM(list_Net,
#>                       v_distrib = v_distrib,
#>                       v_Kmin = 1,
#>                       v_Kmax = 10,
#>                       v_Kinit = NULL,
#>                       verbose = TRUE,
#>                       save=FALSE, initBM = TRUE)

