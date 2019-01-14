


#
# nFG <- 4;
namesFG <- c('famers','plants')
v_K <- c(3,2)
v_NQ <- c(30,37)
E = matrix(1,2,2)
E[1,] = c(1,1)
E[2,] = c(1,2)
nnet <- 2
typeInter = c('diradj','inc')
# #
#
# #---------------------- parameters
# list_pi = list(c(0.31,0.42,0.27),c(0.65,0.35))
# list_theta   <- lapply(1:nnet,function(i){matrix(0,v_K[E[i,1]],v_K[E[i,2]])})
# list_theta[[1]][1,] = c(0.025,0.123,0.053)
# list_theta[[1]][2,]  = c(0.159, 0.300,0.070)
# list_theta[[1]][3,] =  c(0.374, 0.585,  0.357)
# list_theta[[2]][1,] = c(0.186 ,0.653)
# list_theta[[2]][2,] = c(0.559,0.905 )
# list_theta[[2]][3,] = c( 0.390,0.696)
# #
# save(list_theta,list_pi,file = '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_MIRES/param1/paramSim1.Rdata')
load(file = '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_MIRES/param1/paramSim1.Rdata')
# #------------------------SIMULATION -----------------------------
nSimu = 100 ;
v_distrib = rep('bernoulli',2)
dirSaveSimuData <- '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_MIRES/param1/data'
# #
# for (i in 1:nSimu) {
#   datasim <- rMBM(v_NQ ,E , typeInter = typeInter, v_distrib = v_distrib, list_pi  = list_pi, list_theta = list_theta, seed = NULL, namesFG = namesfg,keepClassif = TRUE)
#   namedata  <- paste('datasim',i,'.Rdata',sep = "")
#   save(datasim,file = paste(dirSaveSimuData,namedata,sep = '/' ))
# }


#---------------- ESTIMATION ------------------------------------
dirSaveSimuRes <- '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_MIRES/param1/res_MBM'
for (i in 1:nSimu)
{
  print(i)
  #load data
  namedata  <- paste('datasim',i,'.Rdata',sep = "")
  load(file = paste(dirSaveSimuData,namedata,sep = '/' ))
  listNet <- datasim$list_Net
  # estim
  res_estim <- multipartiteBM(listNet, namesFG = namesFG, v_distrib = v_distrib , v_Kmin = 1 , v_Kmax = 10 , v_Kinit = c(1,1), initBM = TRUE, save = FALSE , verbose = FALSE,nbCores = 10)
  nameres  <- paste('resMBM_',i,'.Rdata',sep = "")
  save(res_estim, file = paste(dirSaveSimuRes,nameres,sep = '/' ))
}
#

#-----------------------------------------------------
#------------------- EXPLOITATION -------------------
#------------------------------------------------------
#
# truevK <- vK_estim <- resCompar <- matrix(0,nSimu,2)
# for (i in 1:nSimu)
# {
#   print(i)
#   #load data
#   namedata  <- paste('datasim',i,'.Rdata',sep = "")
#   load(file = paste(dirSaveSimuData,namedata,sep = '/' ))
#   #listNet <- datasim$list_Net
#   trueZ <- datasim$classif
#   truevK[i,] <- vapply(1:2,function(q){length(unique(trueZ[[q]]))},1)
#   nameres  <- paste('resMBM_',i,'.Rdata',sep = "")
#   load(file = paste(dirSaveSimuRes,nameres,sep = '/' ))
#   estimZ <- res_estim$fittedModel[[1]]$paramEstim$Z
#   resCompar[i,] <-  comparClassif(trueZ,estimZ)
#   vK_estim[i,] <- res_estim$fittedModel[[1]]$paramEstim$v_K
# }

# vK_estim - truevK
# summary(resCompar)
# par(mfrow=c(2,2))
# for (q in 1:2){plot(density(resCompar[,q])); abline(v=1,col='red')}

