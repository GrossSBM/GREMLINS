rm(list = ls())
library(gtools)
library(GREMLIN)

myseed  = sample(1:10000,1)
set.seed(myseed)


n_FG <- sample(3:5,1)
n_FG <- 2
v_K <- sample(2:3,n_FG,replace = TRUE)
names(v_K) <- LETTERS[1:n_FG]
list_pi <- lapply(1:n_FG,function(i) {rdirichlet(1,rep(2,v_K[i]))} )

n_net <- sample(3:5,1)
n_net = 2
v_distrib <- sample(c('bernoulli','poisson'),n_net,replace = TRUE)
#vdistrib <- rep('poisson',n_net)
#vdistrib <- rep('bernoulli',n_net)
typeInter <- rep('NA',n_net)


### chose n_net graphes involving all the functional groups.
stop = FALSE
while (!stop) {
  U <- combinations(n_FG, 2, v = 1:n_FG, set = TRUE, repeats.allowed = TRUE)
  list_net <- sample(1:nrow(U),n_net,replace = FALSE)
  E = U[list_net,]
  ord <- order(E[,1],E[,2])
  E <- E[ord,]
  if (length(unique(c(E))) == n_FG) {stop = TRUE}
}


for (i in 1:n_net) {
  if (E[i,1] == E[i,2]) {typeInter[i] = sample(c('diradj','adj'),1)}  else{typeInter[i] = 'inc' }
}



list_theta <- list()
for (i in 1:n_net) {
  if (v_distrib[i] == 'bernoulli') {
    u <- c(1:(2 * v_K[E[i,1]] * v_K[E[i,2]]))/(2 * v_K[E[i,1]]*v_K[E[i,2]])
    u <- round(u[-c(1,length(u))],4)
    mu <- sample(u,v_K[E[i,1]] * v_K[E[i,2]],replace = FALSE)
    list_theta[[i]] <- matrix(mu,nrow = v_K[E[i,1]], ncol = v_K[E[i,2]] )
    #matrix(rbeta(v_K[E[i,1]] * v_K[E[i,2]],alpha,beta ),nrow = v_K[E[i,1]], ncol = v_K[E[i,2]] )
  }else{
    mean_theta_i = sample(c(1:(2 * v_K[E[i,1]] * v_K[E[i,2]])),v_K[E[i,1]] * v_K[E[i,2]],replace = FALSE)
    list_theta[[i]] <- matrix(rgamma(v_K[E[i,1]] * v_K[E[i,2]],mean_theta_i*10,10 ),nrow = v_K[E[i,1]], ncol = v_K[E[i,2]] )
  }
  if (typeInter[i] == 'adj') {list_theta[[i]] <- 0.5*(list_theta[[i]] + t(list_theta[[i]]))} # symetrisation
}

# nb of individuals

v_NQ = sample(c(50:100),n_FG,replace = TRUE)
data_sim <- rMBM(v_NQ ,E , typeInter, v_distrib, list_pi, list_theta, seed = NULL, namesFG = LETTERS[1:length(v_NQ)])


list_Net <- data_sim$list_net
#vdistrib[1] = 'poisson'
length(list_Net)


res <- multipartiteBM(list_Net, namesFG = LETTERS[1:length(v_NQ)], v_distrib = v_distrib , v_Kmin = 1 , v_Kmax = rep(6,n_FG) , v_Kinit = NULL, initBM = TRUE, save = FALSE , verbose = TRUE,nbCores = 10)


