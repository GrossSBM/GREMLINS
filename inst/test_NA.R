
library(GREMLIN)

## ----param FG--------------------------------------------------------------------------------------------------------------------------------
n_FG <- 3  #number of functional groups
v_K <- c(3,2,2) #number of clusters in each functional group
list_pi <- vector("list", 3); # parameters of clustering
list_pi[[1]] <- c(0.4,0.3,0.3); list_pi[[2]] <- c(0.6,0.4); list_pi[[3]]  <- c(0.6,0.4)
E  = rbind(c(1,2),c(2,3),c(2,2),c(1,3))
v_distrib <- c('bernoulli','poisson','ZIgaussian','gaussian')
typeInter <- c( "inc", "inc"  ,  "adj" ,"inc")

# network 1 bernoulli
list_theta <- list()
theta1 <- rbeta(v_K[E[1,1]] * v_K[E[1,2]],0.9,0.9 )
list_theta[[1]] <- matrix(theta1,nrow = v_K[E[1,1]], ncol = v_K[E[1,2]] )

# network 2 poisson
theta2 <- rgamma(v_K[E[2,1]] * v_K[E[2,2]],7.5,1 )
list_theta[[2]] <- matrix(theta2,nrow = v_K[E[2,1]], ncol = v_K[E[2,2]] )

# network 3  ZI gaussian non dirigé
list_theta[[3]] <- list()
list_theta[[3]]$mean  <- matrix(round(rnorm(v_K[E[3,1]] * v_K[E[3,1]],7.5,4 ),1),nrow = v_K[E[3,1]], ncol = v_K[E[3,2]] )
list_theta[[3]]$var  = matrix(round(rgamma(v_K[E[3,1]] * v_K[E[3,2]],7.5,4 ),1),nrow = v_K[E[3,1]], ncol = v_K[E[3,2]] )
list_theta[[3]]$p0  = matrix(round(rbeta(v_K[E[3,1]] * v_K[E[3,2]],2,2 ),1),nrow = v_K[E[3,1]], ncol = v_K[E[3,2]] )
list_theta[[3]] = lapply(list_theta[[3]] , function(u){0.5 * (u + t(u))}) # for symetrisation

# network 3  Gaussian
list_theta[[4]] <- list()
theta4_mean <- rnorm(v_K[E[4,1]] * v_K[E[4,2]],7.5,1 )
theta4_var <- rgamma(v_K[E[4,1]] * v_K[E[4,2]],7.5,1 )
list_theta[[4]]$mean <- matrix(theta4_mean,nrow = v_K[E[4,1]], ncol = v_K[E[4,2]] )
list_theta[[4]]$var <- matrix(theta4_var,nrow = v_K[E[4,1]], ncol = v_K[E[4,2]] )



v_NQ <-  c(100,50,40)
namesFG <-  c('A','B','C')
dataSim <- rMBM(v_NQ = v_NQ , E = E , typeInter = typeInter,
                v_distrib = v_distrib, list_pi = list_pi,
                list_theta = list_theta, namesFG,keepClassif = TRUE)
list_Net <- dataSim$list_Net
diag(list_Net[[3]]$mat) = NA
classifTrue <- dataSim$classif


## ----MBM simul eval false, echo = TRUE, eval = FALSE-----------------------------------------------------------------------------------------
  res_MBMsimu <- multipartiteBM(list_Net,
                      v_distrib = v_distrib,
                      v_Kmin = 1,
                      v_Kmax = 10,
                      v_Kinit = NULL,
                      verbose = TRUE,
                      maxiterVE =  100,maxiterVEM =  100,
                      save=FALSE, initBM = TRUE)

for (i in 1:n_FG){print(table(res_MBMsimu$fittedModel[[1]]$paramEstim$Z[[i]],dataSim$classif[[i]]))}

table(res_MBMsimu$fittedModel[[1]]$paramEstim$Z[[1]])



############# NA data at random in any matrix
epsilon =  10/100
list_Net_NA <- list_Net
for (m in 1:nrow(E)){
   U <-  sample(c(1,0),v_NQ[E[m,1]]*v_NQ[E[m,2]],replace=TRUE,prob  = c(epsilon, 1-epsilon))
   matNA <- matrix(U,v_NQ[E[m,1]],v_NQ[E[m,2]])
   list_Net_NA[[m]]$mat[matNA== 1] = NA
   if (list_Net_NA[[m]]$typeInter == 'adj') {
     M <- list_Net_NA[[m]]$mat
     diag(M) <- NA
     M[lower.tri(M)] = t(M)[lower.tri(M)]
     list_Net_NA[[m]]$mat <- M
     }
}

res_MBMsimu_NA <- multipartiteBM(list_Net_NA,
                              v_distrib = v_distrib,
                              v_Kmin = 1,
                              v_Kmax = 10,
                              v_Kinit = NULL,
                              verbose = TRUE,
                              maxiterVE =  100,maxiterVEM =  100,
                              save=FALSE, initBM = TRUE)

for (i in 1:n_FG){print(table(res_MBMsimu_NA$fittedModel[[1]]$paramEstim$Z[[i]],dataSim$classif[[i]]))}
