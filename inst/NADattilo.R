
# package and data --------------------------------------------------------

library(GREMLIN)
library(alluvial)
library(pROC)
library(parallel)


setwd("~/Dropbox/Multiplex/Ecologie/Code/generalized_multi_BM/realdata/")

#data
dattilo = read.csv(file = 'dattilo.csv',header = TRUE,skip = 1)
plantnames = dattilo[,2]

dattilo2 = read.csv(file = 'dattilo.csv',header = FALSE)
type_animals = dattilo2[1,]
type_animals = as.factor(as.character(type_animals[-c(1,2)]))

indflovis = which(dattilo2[1,] == "Floral visitor")
indant = which(dattilo2[1,] == "Ant")
indseed = which(dattilo2[1,] == "Seed dispersal")

Adj_flovis = dattilo[,indflovis]
rownames(Adj_flovis) = plantnames
Adj_flovis = as.matrix(Adj_flovis)
Adj_ant = dattilo[,indant]
rownames(Adj_ant) = plantnames
Adj_ant = as.matrix(Adj_ant)
Adj_seed = dattilo[,indseed]
rownames(Adj_seed) = plantnames
Adj_seed = as.matrix(Adj_seed)


# formatting data for GREMLIN ---------------------------------------------

PlantFlovis = defineNetwork(Adj_flovis,"inc","plants","flovis")
PlantAnt = defineNetwork(Adj_ant,"inc","plants","ants")
PlantBird = defineNetwork(Adj_seed,"inc","plants","birds")


# NB rep and options

NREP = 200
Nmiss = 50




# Case 1 NA on plant Flovis -----------------------------------------------------

REPAUC = mclapply(1:NREP,function(ite){

  # simulation NAs in relation matrix
  Adj_flovisNA = Adj_flovis
  wNA = sample(length(Adj_flovisNA),Nmiss)
  Adj_flovisNA[wNA] = NA
  PlantFlovisNA = defineNetwork(Adj_flovisNA,"inc","plants","flovis")


  # multipartite
  RESMBM = multipartiteBM(list(PlantFlovisNA,PlantAnt,PlantBird),namesFG = c("plants","flovis","ants","birds"),v_Kmin = 1,v_Kmax = 10,save=F,initBM=TRUE)
  PRED = predictMBM(RESMBM)
  wNA = PRED[[1]]$wNA
  prob = PRED[[1]]$mat[wNA]
  true = Adj_flovis[wNA]
  RESROC = pROC::roc(response=true,predictor=prob)
  aucMBM = RESROC$auc
  #plot(RESROC)

  ## unipartite sans inventaire
  RESBM = multipartiteBM(list(PlantFlovisNA),namesFG = c("plants","flovis"),v_Kmin = 1,v_Kmax = 10,save=F,initBM=TRUE)
  PRED2 = predictMBM(RESBM)
  wNA2 = PRED2[[1]]$wNA
  prob2 = PRED2[[1]]$mat[wNA2]
  RESROC2 = pROC::roc(response=true,predictor=prob2)
  aucBM = RESROC2$auc

  return(c(aucMBM,aucBM))
},mc.cores = 4)

AUCmissFlovis = Reduce(rbind,REPAUC)

save(AUCmissFlovis,file="AUCnaDattilo.Rdata")



# Case 2 NA on plant ant --------------------------------------------------


REPAUC = mclapply(1:NREP,function(ite){

  # simulation NAs in relation matrix
  Adj_antNA = Adj_ant
  wNA = sample(length(Adj_antNA),Nmiss)
  Adj_antNA[wNA] = NA
  PlantAntNA = defineNetwork(Adj_antNA,"inc","plants","ants")


  # multipartite
  RESMBM = multipartiteBM(list(PlantFlovis,PlantAntNA,PlantBird),namesFG = c("plants","flovis","ants","birds"),v_Kmin = 1,v_Kmax = 10,save=F,initBM=TRUE)
  PRED = predictMBM(RESMBM)
  wNA = PRED[[2]]$wNA
  prob = PRED[[2]]$mat[wNA]
  true = Adj_ant[wNA]
  RESROC = pROC::roc(response=true,predictor=prob)
  aucMBM = RESROC$auc
  #plot(RESROC)

  ## unipartite sans inventaire
  RESBM = multipartiteBM(list(PlantAntNA),namesFG = c("plants","ants"),v_Kmin = 1,v_Kmax = 10,save=F,initBM=TRUE)
  PRED2 = predictMBM(RESBM)
  wNA2 = PRED2[[1]]$wNA
  prob2 = PRED2[[1]]$mat[wNA2]
  RESROC2 = pROC::roc(response=true,predictor=prob2)
  aucBM = RESROC2$auc

  return(c(aucMBM,aucBM))
},mc.cores = 4)

AUCmissAnt = Reduce(rbind,REPAUC)

save(AUCmissFlovis,AUCmissAnt,file="AUCnaDattilo.Rdata")




# Case 3 NA on plant bird --------------------------------------------------


REPAUC = mclapply(1:NREP,function(ite){

  # simulation NAs in relation matrix
  Adj_seedNA = Adj_seed
  wNA = sample(length(Adj_seedNA),Nmiss)
  Adj_seedNA[wNA] = NA
  PlantBirdNA = defineNetwork(Adj_seedNA,"inc","plants","birds")


  # multipartite
  RESMBM = multipartiteBM(list(PlantFlovis,PlantAnt,PlantBirdNA),namesFG = c("plants","flovis","ants","birds"),v_Kmin = 1,v_Kmax = 10,save=F,initBM=TRUE)
  PRED = predictMBM(RESMBM)
  wNA = PRED[[2]]$wNA
  prob = PRED[[2]]$mat[wNA]
  true = Adj_seed[wNA]
  RESROC = pROC::roc(response=true,predictor=prob)
  aucMBM = RESROC$auc
  #plot(RESROC)

  ## unipartite sans inventaire
  RESBM = multipartiteBM(list(PlantBirdNA),namesFG = c("plants","birds"),v_Kmin = 1,v_Kmax = 10,save=F,initBM=TRUE)
  PRED2 = predictMBM(RESBM)
  wNA2 = PRED2[[1]]$wNA
  prob2 = PRED2[[1]]$mat[wNA2]
  RESROC2 = pROC::roc(response=true,predictor=prob2)
  aucBM = RESROC2$auc

  return(c(aucMBM,aucBM))
},mc.cores = 4)

AUCmissBird = Reduce(rbind,REPAUC)

save(AUCmissFlovis,AUCmissAnt,AUCmissBird,file="AUCnaDattilo.Rdata")
