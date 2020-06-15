
# package and data --------------------------------------------------------

library(GREMLIN)
library(alluvial)
library(pROC)
library(parallel)
#data
setwd("~/Dropbox/Multiplex/Ecologie/Code/generalized_multi_BM/realdata/MIRES/")
relations_MIRES = read.csv(file ="ssna_adj_30n_all_ind_ssloop.csv",header=TRUE,skip=0)
row.names(relations_MIRES) <- relations_MIRES[,1]
relations_MIRES <- as.matrix(relations_MIRES[,-1])
n = dim(relations_MIRES)[1]
diag(relations_MIRES) <- 0
adj_relations<- relations_MIRES
for (i in 1:n){for (j in 1:n){adj_relations[i,j] <- as.numeric(relations_MIRES[i,j]>0)}}
#

#############" Qui cultive quoi
species <- read.csv(file ="inventory_species.csv",header=TRUE)
rownames(species) <- species[,1]
species <- as.matrix(species[,-1])

for (i in 1:dim(species)[1]){for (j in 1:dim(species)[2]){species[i,j] <- as.numeric(species[i,j]>0)}}
# Verification de l'adequation de sdeux matrices (memes lignes)
rownames(adj_relations)==rownames(species)



# formatting data for GREMLIN ---------------------------------------------


Relation = defineNetwork(adj_relations,"diradj","ind","ind")
Inventaire = defineNetwork(species,"inc","ind","spe")



# NB rep and options

NREP = 200
Nmiss = 50



# Case 1 missing data on farmers' relations -------------------------------

REPAUC = mclapply(1:NREP,function(ite){

  # simulation NAs in relation matrix
  adj_relationsNA = adj_relations
  diag(adj_relationsNA) =NA # NA on the diagonal since there are no loop
  wNA = sample(length(adj_relationsNA),Nmiss)
  adj_relationsNA[wNA] = NA
  RelationNA = defineNetwork(adj_relationsNA,"diradj","ind","ind")


  # multipartite
  RESMBM = multipartiteBM(list(RelationNA,Inventaire),namesFG = c("ind","spe"),v_Kmin = 1,v_Kmax = 5,save=F,initBM=TRUE)
  PRED = predictMBM(RESMBM)
  wNA = PRED[[1]]$wNA
  wNA = wNA[!(abs(wNA[,1]-wNA[,2])==0),]
  prob = PRED[[1]]$mat[wNA]
  true = adj_relations[wNA]
  RESROC = pROC::roc(response=true,predictor=prob)
  aucMBM = RESROC$auc
  #plot(RESROC)

  ## unipartite sans inventaire
  RESBM = multipartiteBM(list(RelationNA),namesFG = c("ind"),v_Kmin = 1,v_Kmax = 5,save=F,initBM=TRUE)
  PRED2 = predictMBM(RESBM)
  wNA2 = PRED2[[1]]$wNA
  wNA2 = wNA2[!(abs(wNA2[,1]-wNA2[,2])==0),]
  prob2 = PRED2[[1]]$mat[wNA2]
  RESROC2 = pROC::roc(response=true,predictor=prob2)
  aucREL = RESROC2$auc

  return(c(aucMBM,aucREL))
},mc.cores = 4)

AUCmissRel = Reduce(rbind,REPAUC)




# Case 2 missing data on inventory ------------------------------------------------

REPAUC = mclapply(1:NREP,function(ite){

  # simulation NAs in inventory matrix
  speciesNA = species
  wNA = sample(length(species),Nmiss)
  speciesNA[wNA] = NA
  InventaireNA = defineNetwork(speciesNA,"inc","ind","spe")


  # multipartite
  RESMBM = multipartiteBM(list(Relation,InventaireNA),namesFG = c("ind","spe"),v_Kmin = 1,v_Kmax = 5,save=F,initBM=TRUE)
  PRED = predictMBM(RESMBM)
  wNA = PRED[[2]]$wNA
  prob = PRED[[2]]$mat[wNA]
  true = species[wNA]
  RESROC = pROC::roc(response=true,predictor=prob)
  aucMBM = RESROC$auc
  #plot(RESROC)

  ## unipartite sans inventaire
  RESBM = multipartiteBM(list(InventaireNA),namesFG = c("ind","spe"),v_Kmin = 1,v_Kmax = 5,save=F,initBM=TRUE)
  PRED2 = predictMBM(RESBM)
  wNA2 = PRED2[[1]]$wNA
  prob2 = PRED2[[1]]$mat[wNA2]
  RESROC2 = pROC::roc(response=true,predictor=prob2)
  aucINV = RESROC2$auc

  return(c(aucMBM,aucINV))
},mc.cores = 4)

AUCmissInv = Reduce(rbind,REPAUC)


# Case 3 missing data on both relation and inventory ----------------------

REPAUC = mclapply(1:NREP,function(ite){

  adj_relationsNA = adj_relations
  diag(adj_relationsNA) =NA
  wNA = sample(length(adj_relationsNA),Nmiss)
  adj_relationsNA[wNA] = NA
  RelationNA = defineNetwork(adj_relationsNA,"diradj","ind","ind")

  speciesNA = species
  wNA = sample(length(species),Nmiss)
  speciesNA[wNA] = NA
  InventaireNA = defineNetwork(speciesNA,"inc","ind","spe")



  RESMBM = multipartiteBM(list(RelationNA,InventaireNA),namesFG = c("ind","spe"),v_Kmin = 1,v_Kmax = 5,save=F,initBM=TRUE)

  PRED = predictMBM(RESMBM)
  wNA = PRED[[1]]$wNA
  wNA = wNA[!(abs(wNA[,1]-wNA[,2])==0),]
  prob = PRED[[1]]$mat[wNA]
  wNA2 = PRED[[2]]$wNA
  prob2 = PRED[[2]]$mat[wNA2]

  true = adj_relations[wNA]
  true2 = speciesNA[wNA2]

  RESROC = pROC::roc(response=c(true,true2),predictor=c(prob,prob2))
  aucMBM = RESROC$auc


  RESRel = multipartiteBM(list(RelationNA),namesFG = c("ind"),v_Kmin = 1,v_Kmax = 5,save=T,initBM=TRUE)
  PREDRel = predictMBM(RESRel)
  wNA = PREDRel[[1]]$wNA
  wNA = wNA[!(abs(wNA[,1]-wNA[,2])==0),]
  prob = PREDRel[[1]]$mat[wNA]
  RESInv = multipartiteBM(list(InventaireNA),namesFG = c("ind","spe"),v_Kmin = 1,v_Kmax = 5,save=T,initBM=TRUE)
  PREDInv = predictMBM(RESInv)
  wNA = PREDInv[[1]]$wNA
  prob2 = PREDInv[[1]]$mat[wNA]
  RESROC = pROC::roc(response=c(true,true2),predictor=c(prob,prob2))
  aucSEP = RESROC$auc

  return(c(aucMBM,aucSEP))
},mc.cores=4)

AUCmissBoth = Reduce(rbind,REPAUC)



#save(AUCmissRel,AUCmissInv,AUCmissBoth,file="AUCnaMIRES.Rdata")


# plotting results --------------------------------------------------------


library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
dfAUCmissRel = data.frame(AUCmissRel)
names(dfAUCmissRel) = c("All","OnlyRel")
dfAUCmissRel = melt(dfAUCmissRel,value.name = "AUC",variable.name = "type")
pRel = ggplot(dfAUCmissRel,aes(x=type,y=AUC)) + geom_boxplot() + ggtitle(paste("MBM better than single BM :",100*mean(AUCmissRel[,1]>AUCmissRel[,2]),"%"))

dfAUCmissInv = data.frame(AUCmissInv)
names(dfAUCmissInv) = c("All","OnlyInv")
dfAUCmissInv = melt(dfAUCmissInv,value.name = "AUC",variable.name = "type")
pInv = ggplot(dfAUCmissInv,aes(x=type,y=AUC)) + geom_boxplot() + ggtitle(paste("MBM better than single BM :",100*mean(AUCmissInv[,1]>AUCmissInv[,2]),"%"))


dfAUCmissBoth = data.frame(AUCmissBoth)
names(dfAUCmissBoth) = c("All","SingleBoth")
dfAUCmissBoth = melt(dfAUCmissBoth,value.name = "AUC",variable.name = "type")
pBoth = ggplot(dfAUCmissBoth,aes(x=type,y=AUC)) + geom_boxplot() + ggtitle(paste("MBM better than single BM :",100*mean(AUCmissBoth[,1]>AUCmissBoth[,2]),"%"))


grid.arrange(pRel,pInv,pBoth,nrow=1)

