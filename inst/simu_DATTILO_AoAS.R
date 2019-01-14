#rm(list = ls())

namesFG <- c('plants','pollinators','birds','ants')
nFG <- 4;
v_K <- c(7,2,2,1)
v_NQ <- c(141,173,46,30)
E = matrix(1,nFG - 1,2)
E[,2] = c(2,3,4)
nnet =  3
typeInter = c('inc','inc','inc')
v_distrib = rep('bernoulli',3)


#---------------------- parameters : param1

# ### param1
# list_pi = list(c(0.4675,0.1606, 0.1351,0.0784,0.1061, 0.0142,0.0381),c(0.06,0.94),c(0.1,0.9),c(1))
#
# # ## param 2
# # list_pi[[1]] = c(0.46,0.16, 0.15,0.184,0.1061, 0.1,0.1)
# # list_pi[[1]] = list_pi[[1]]/sum(list_pi[[1]])
#
#
# list_theta   <- lapply(1:nnet,function(i){matrix(0,v_K[E[i,1]],v_K[E[i,2]])})
# list_theta[[1]][1,] = c(0.0957,0.0075)
# list_theta[[1]][2,1]  = 0.01
# list_theta[[1]][3,2] = 0.0003
# list_theta[[1]][4,] = c(0.1652, 0.0343)
# list_theta[[1]][5,] = c(0.2018,0.138)
#
# list_theta[[2]][1,2] = 0.0006
# list_theta[[2]][2,] = c(0.5431, 0.0589)
# list_theta[[2]][4,] = c(0.6620,0.1542)
# list_theta[[2]][7,] = c(0.8492,0.3565)
#
# list_theta[[3]][1,] = c(0.0013)
# list_theta[[3]][3,] = 0.0753
# list_theta[[3]][5,] = 0.0163
# list_theta[[3]][6,] = 0.5108
#
# save(list_theta,list_pi,file = '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_DATTILO/param1/paramSim1.Rdata')


#---------------------- parameters : param2
lpi = list(c(0.4675,0.1606, 0.1351,0.0784,0.1061, 0.0142,0.0381),c(0.06,0.94),c(0.1,0.9),c(1))
lpi[[2]] = c(0.1,0.9)
lpi[[1]] = c(0.46,0.16, 0.15,0.184,0.1061, 0.1,0.1)
lpi[[1]] = lpi[[1]]/sum(lpi[[1]])
ltheta   <- lapply(1:nnet,function(i){matrix(0,v_K[E[i,1]],v_K[E[i,2]])})
ltheta[[1]][1,] = c(0.0957,0.0075)
ltheta[[1]][2,1]  = 0.01
ltheta[[1]][3,2] = 0.0003
ltheta[[1]][4,] = c(0.1652, 0.0343)
ltheta[[1]][5,] = c(0.2018,0.138)

ltheta[[2]][1,2] = 0.0006
ltheta[[2]][2,] = c(0.5431, 0.0589)
ltheta[[2]][4,] = c(0.6620,0.1542)
ltheta[[2]][7,] = c(0.8492,0.3565)

ltheta[[3]][1,] = c(0.0013)
ltheta[[3]][3,] = 0.0753
ltheta[[3]][5,] = 0.0163
ltheta[[3]][6,] = 0.5108

list_pi <- lpi
list_theta <- ltheta
save(ltheta,lpi,file = '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_DATTILO/param2/paramSim2.Rdata')

##########  ------------------SIMULATION -----------------------------
nSimu = 100 ;
#
 dirSaveSimuData <- '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_DATTILO/param2/data'
#
for (i in 1:nSimu){
    datasim <- rMBM(v_NQ = v_NQ ,E  = E , typeInter =  typeInter, v_distrib, list_pi = list_pi, list_theta = list_theta, seed = NULL, namesFG = namesFG,keepClassif = TRUE)
    namedata  <- paste('datasim',i,'.Rdata',sep = "")
    save(datasim,file = paste(dirSaveSimuData,namedata,sep = '/' ))
}


#
# #---------------- ESTIMATION ------------------------------------
dirSaveSimuRes <- '/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Avner-Pierre/Ecologie/Code/generalized_multi_BM/res_simu_AoAS/res_simu_AoAS_DATTILO/param2/res_MBM'
for (i in 1:nSimu)
 {
 print(i)
 #load data
 namedata  <- paste('datasim',i,'.Rdata',sep = "")
 load(file = paste(dirSaveSimuData,namedata,sep = '/' ))
 list_Net <- datasim$list_Net
 # estim
 res_estim <- multipartiteBM(list_Net, namesFG = namesFG, v_distrib = v_distrib , v_Kmin = 1 , v_Kmax = 10 , v_Kinit = c(1,1,1,1), initBM = TRUE, save = FALSE , verbose = FALSE,nbCores = 10)
 nameres  <- paste('resMBM_',i,'.Rdata',sep = "")
 save(res_estim, file = paste(dirSaveSimuRes,nameres,sep = '/' ))
}


# #-----------------------------------------------------
# #------------------- EXPLOITATION -------------------
# #------------------------------------------------------
# #
# truev_K <- v_K_estim <- resCompar <- matrix(0,nSimu,4)
# for (i in 1:nSimu)
# {
#   print(i)
#   #load data
#   namedata  <- paste('datasim',i,'.Rdata',sep = "")
#   load(file = paste(dirSaveSimuData,namedata,sep = '/' ))
#   trueZ <- datasim$classif
#   truev_K[i,] <- vapply(1:4,function(q){length(unique(trueZ[[q]]))},1)
#   nameres  <- paste('resMBM_',i,'.Rdata',sep = "")
#   load(file = paste(dirSaveSimuRes,nameres,sep = '/' ))
#   estimZ <- res_estim$fitted.model[[1]]$param_estim$Z
#   resCompar[i,] <-  comparClassif(trueZ,estimZ)
#   v_K_estim[i,] <- res_estim$fitted.model[[1]]$param_estim$v_K
# }
#
# v_K_estim - truev_K
# summary(resCompar)
# par(mfrow=c(2,2))
# for (q in 1:4){plot(density(resCompar[,q])); abline(v=1,col='red')}

