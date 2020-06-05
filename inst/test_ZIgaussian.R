
library(GREMLIN)

## ----param FG--------------------------------------------------------------------------------------------------------------------------------

n_FG <- 1 #number of functional groups (FG)
v_NQ <-  c(60) #size og each FG
v_K  <- c(3) #number of clusters in each functional group
list_pi = lapply(1:n_FG,function(q){v = rgamma(v_K[q],1,1); return(v/sum(v))}) #proportion of each block in each eahc FG
namesFG <- c('A')


## ----param 1 net-----------------------------------------------------------------------------------------------------------------------------
E  <-  rbind(c(1,1))
typeInter <- c( "adj")
v_distrib <- c('ZIgaussian')


## ----param net-------------------------------------------------------------------------------------------------------------------------------
list_theta <- list()
list_theta[[1]] <- list()
list_theta[[1]]$mean  <- matrix(round(rnorm(v_K[E[1,1]] * v_K[E[1,2]],7.5,4 ),1),nrow = v_K[E[1,1]], ncol = v_K[E[1,2]] )
list_theta[[1]]$var  = matrix(round(rgamma(v_K[E[1,1]] * v_K[E[1,2]],7.5,4 ),1),nrow = v_K[E[1,1]], ncol = v_K[E[1,2]] )
list_theta[[1]]$p0  = matrix(round(rbeta(v_K[E[1,1]] * v_K[E[1,2]],2,2 ),1),nrow = v_K[E[1,1]], ncol = v_K[E[1,2]] )
list_theta[[1]] = lapply(list_theta[[1]] , function(u){0.5 * (u + t(u))}) # for symetrisation


# ----simul 2, eval = FALSE, echo = TRUE------------------------------------------------------------------------------------------------------
 dataSim <- rMBM(v_NQ ,E , typeInter, v_distrib, list_pi,
                 list_theta, seed = NULL,
                 namesFG = namesFG,keepClassif  = TRUE)
 list_Net <- dataSim$list_Net
 diag(list_Net[[1]]$mat) = NA
 length(list_Net)
 names(list_Net[[1]])
 list_Net[[1]]$typeInter
 list_Net[[1]]$rowFG
 list_Net[[1]]$colFG


## ---- verif---------------------------------------------------------------------------------------------------------
# Y <- dataSim$list_Net[[1]]$mat
# Z = dataSim$classif[[1]]
# pi_true <- table(Z)/length(Z)
# tau_true <- matrix(0,v_NQ[1],v_K[1])
# for (i in 1:v_NQ[1]){tau_true[i,Z[i]] = 1}
# Unit <- matrix(1,)
# Denom <- crossprod(crossprod(Unit, tau[[gr]]), tau[[gc]])
# mu <- crossprod(crossprod(list_Mat[[e]], tau[[gr]]), tau[[gc]])
# Zeros_e  <- Y == 0
# mean_true <- mu / crossprod(crossprod(1 - Zeros_e, tau[[gr]]), tau[[gc]])
# A <- crossprod(crossprod(list_Mat[[e]]^2, tau[[gr]]), tau[[gc]]) /  crossprod(crossprod(1-Zeros_e, tau[[gr]]), tau[[gc]])
# list_theta_e$var <-  A - list_theta_e$mean^2
# list_theta_e$p0 <- crossprod(crossprod(Zeros_e, tau[[gr]]), tau[[gc]]) / Denom


## ----MBM simul eval false, echo = TRUE, eval = FALSE-----------------------------------------------------------------------------------------
  res_MBMsimu <- multipartiteBM(list_Net,
                      v_distrib = v_distrib,
                      v_Kmin = 1,
                      v_Kmax = 10,
                      v_Kinit = NULL,
                      verbose = TRUE,
                      maxiterVE =  100,maxiterVEM =  100,
                      save=FALSE, initBM = TRUE)


## ----estim param, eval=FALSE-----------------------------------------------------------------------------------------------------------------
res_MBMsimu$fittedModel[[1]]$paramEstim$list_theta$AA$mean
list_theta[[1]]$mean

res_MBMsimu$fittedModel[[1]]$paramEstim$list_theta$AA$var
list_theta[[1]]$var

res_MBMsimu$fittedModel[[1]]$paramEstim$list_theta$AA$p0
list_theta[[1]]$p0

res_MBMsimu$fittedModel[[1]]$paramEstim$list_pi
pi_true

table(res_MBMsimu$fittedModel[[1]]$paramEstim$Z[[1]],Z)

##

