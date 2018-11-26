rm(list = ls())
library(gtools)
library(GREMLIN)

v_K <- c(3,2,2)
n_FG <- length(v_K)
lpi <- lapply(1:n_FG,function(i) {rdirichlet(1,rep(2,v_K[i]))} )

n_net = 3
vdistrib <- sample(c('bernoulli','poisson'),n_net,replace = TRUE)
#vdistrib <- rep('poisson',n_net)
#vdistrib <- rep('bernoulli',n_net)
type_inter <- rep('NA',n_net)


### chose n_net graphes involving all the functional groups.
stop = FALSE
while (!stop) {
  U <- combinations(n_FG, 2, v = 1:n_FG, set = TRUE, repeats.allowed = TRUE)
  list_net <- sample(1:nrow(U),n_net,replace = FALSE)
  E = U[list_net,]
  ord <- order(E[,1])
  E <- E[ord,]
  if (length(unique(c(E))) == n_FG) {stop = TRUE}
}


for (i in 1:n_net) {
  if (E[i,1] == E[i,2]) {type_inter[i] = sample(c('diradj','adj'),1)}  else{type_inter[i] = 'inc' }
}

ltheta <- list()
for (i in 1:n_net) {
  if (vdistrib[i] == 'bernoulli') {
    ltheta[[i]] <- matrix(rbeta(v_K[E[i,1]] * v_K[E[i,2]],1.5,1.5 ),nrow = v_K[E[i,1]], ncol = v_K[E[i,2]] )
  }else{
    ltheta[[i]] <- matrix(rgamma(v_K[E[i,1]] * v_K[E[i,2]],7.5,1 ),nrow = v_K[E[i,1]], ncol = v_K[E[i,2]] )
  }
  if (type_inter[i] == 'adj') { ltheta[[i]] <- 0.5*(ltheta[[i]] + t(ltheta[[i]]))} # symetrisation
}

# nb of individuals

v_NQ = c(100,50,40)
data_sim <- rMBM(v_NQ ,E , type_inter, vdistrib, lpi, ltheta, seed=NULL, namesfg= c('A','B','D'))

listNet <- data_sim
#vdistrib[1] = 'poisson'


res <- MultipartiteBM(listNet, namesFG = c('A','B','D'), vdistrib = vdistrib , vKmin = 1 , vKmax = c(6,6,6) , vKinit = c(3,3,3), init.BM = FALSE , save=FALSE , verbose = TRUE)


