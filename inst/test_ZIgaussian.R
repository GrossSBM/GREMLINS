
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
list_theta_true <- list_theta

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
Y <- list(mat = list(list_Net[[1]]$mat))
Z = dataSim$classif[[1]]
pi_true <- table(Z)/length(Z)
tau_true <- matrix(0,v_NQ[1],v_K[1])
for (i in 1:v_NQ[1]){tau_true[i,Z[i]] = 1}
list_MaskNA <- lapply(Y$mat,function(m){1 * (1 - is.na(m))})
list_Mat <- lapply(Y$mat,function(m){m[is.na(m)] = 0; m})

Denom <- crossprod(crossprod(list_MaskNA[[1]], tau_true), tau_true)
mu <- crossprod(crossprod(list_Mat[[1]], tau_true), tau_true)
NonZeros_e  <- list_Mat[[1]] != 0
list_theta_sim = list()
list_theta_sim$mean <- mu / crossprod(crossprod(NonZeros_e, tau_true), tau_true)
A <- crossprod(crossprod(list_Mat[[1]]^2, tau_true), tau_true) /  crossprod(crossprod(NonZeros_e, tau_true), tau_true)
list_theta_sim$var <-  A - list_theta_sim$mean^2
list_theta_sim$p0 <- 1 - crossprod(crossprod(NonZeros_e, tau_true), tau_true) / Denom



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
list_theta_sim$mean


res_MBMsimu$fittedModel[[1]]$paramEstim$list_theta$AA$var
list_theta_true[[1]]$var
list_theta_sim$var

res_MBMsimu$fittedModel[[1]]$paramEstim$list_theta$AA$p0
list_theta_true[[1]]$p0
list_theta_sim$p0


res_MBMsimu$fittedModel[[1]]$paramEstim$list_pi
pi_true

table(res_MBMsimu$fittedModel[[1]]$paramEstim$Z[[1]],Z)

##

