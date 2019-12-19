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


## ----loading dataset, eval=TRUE------------------------------------------
load(file = 'data/MultiPartiteMutalistic.Rda')


## ----transform dataset,  eval=TRUE---------------------------------------
PlantFlovis = defineNetwork(Inc_plant_flovis,"inc","plants","flovis")
PlantAnt = defineNetwork(Inc_plant_ant,"inc","plants","ants")
PlantBird = defineNetwork(Inc_plant_bird,"inc","plants","birds")
str(PlantBird)
list_Net <- list(PlantFlovis,PlantAnt,PlantBird)



## ----example of dataset, eval=TRUE---------------------------------------
PlantFlovis$mat[1:2,1:2]



## ----MBM-----------------------------------------------------------------
if (!file.exists('vignettes/resMBM_Mutualistic.Rda')) {
  RES_MBM = multipartiteBM(
    list_Net,
    namesFG = c('plants','flovis','ants','birds'),
    v_distrib  = c('bernoulli','bernoulli','bernoulli'),
    initBM = TRUE,
    save = TRUE)
  save(RES_MBM,file="vignettes/resMBM_Mutualistic.Rda")
} else {load("vignettes/resMBMinitBM.Rdata")}


## ----MBM param-----------------------------------------------------------
RES_MBM[[1]]


## ----example-------------------------------------------------------------
library(GREMLIN)
## basic example code


## ----cars----------------------------------------------------------------
summary(cars)


## ----pressure, echo = FALSE----------------------------------------------
plot(pressure)

