


## ----param FG--------------------------------------------------------------------------------------------------------------------------------
set.seed(302718)
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
 length(list_Net)
 names(list_Net[[1]])
 list_Net[[1]]$typeInter
 list_Net[[1]]$rowFG
 list_Net[[1]]$colFG


## ----simul, eval = TRUE, echo =FALSE---------------------------------------------------------------------------------------------------------
Y <- list_Net[[1]]$mat
Z <- dataSim$classif[[1]]
theta_estim  = list(mean =  matrix(0,v_K[1],v_K[1]), var =  matrix(0,v_K[1],v_K[1]), p0  =  matrix(0,v_K[1],v_K[1]))

for (k in v_K[1]){
  for (l in v_K[1]){

    theta_estim$mean[k,l] <-

  }
}


## ----MBM simul, echo = FALSE, eval = TRUE----------------------------------------------------------------------------------------------------
 # if (!file.exists('resMBM_Simu.Rda')) {
 #   res_MBMsimu <- multipartiteBM(list_Net,
 #                       v_distrib = v_distrib,
 #                       v_Kmin = 1,
 #                       v_Kmax = 10,
 #                       v_Kinit = NULL,
 #                       verbose = TRUE,
 #                       save=FALSE, initBM = TRUE)
 #   save(RES_MBM,file = "resMBM_Simu.Rda")
 # } else {load("resMBM_Simu.Rda")}


## ----MBM simul eval false, echo = TRUE, eval = FALSE-----------------------------------------------------------------------------------------
  res_MBMsimu <- multipartiteBM(list_Net,
                      v_distrib = v_distrib,
                      v_Kmin = 1,
                      v_Kmax = 10,
                      v_Kinit = NULL,
                      verbose = TRUE,
                      save=FALSE, initBM = TRUE, nbCores = 1)


## ----estim param, eval=FALSE-----------------------------------------------------------------------------------------------------------------
##
##
## res_MBMsimu$fittedModel[[1]]$paramEstim$list_theta$AB$mean
## list_theta[[1]]$mean
##

