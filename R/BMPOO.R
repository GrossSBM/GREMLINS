CollInteraction = R6Class("CollInteraction", ### classe objet pour décrire les données
      public = list(
            mats = NULL,    #list of mats of interaction
            E = NULL, # 2-column table giving which fgs (functional groups) interact in corresponding mat
            namesFG = NULL, #vector of names of fgs
            namesInd = NULL,
            typeInter = NULL, #vector of type of matrices (diradj=directed adjacency, adj, inc=incidence)
            v_distrib = NULL, #vector of emission distribution (same length as number mats) (poisson, bernoulli...)
            Q = NULL,     #number of functional groups
            cardE = NULL, #number of mats
            v_NQ = NULL, #number of individuals in fgs
            where = NULL, #useful
            Ecode = NULL,
            initialize = function(mats, E_FG , typeInter,v_distrib = NULL, namesInd = NULL)

              {
                self$cardE = length(mats)
                self$mats = mats
                #if (is.null(namesFG)){self$namesFG=unique(as.vector(E_FG))}else {self$namesFG=namesFG}
                self$namesFG <- unique(as.vector(E_FG))# les fonctional groups names sont toujours ceux dans l'ordre d'apparition dans les matrices

                self$Q <- length(self$namesFG)
                self$typeInter <- typeInter
                self$E <- matrix(sapply(E_FG,function(a){which(self$namesFG == a)}),self$cardE,2,byrow = FALSE)


                #browser()

                #init v_distrib, default Bernoulli
                if (is.null(v_distrib)) {self$v_distrib = rep("bernoulli",self$cardE)}
                else{
                  if (length(v_distrib) != self$cardE) stop("number of distributions not consistent with number of interaction matrices")
                  self$v_distrib = v_distrib
                }


                #chack adequation of v_distrib from data
                v_distrib_guessed <- unlist(lapply(mats,function(Net){
                   support <- sort(unique(as.vector(Net)))
                   if (all(is.poswholenumber(support))) {
                     if (length(support) > 2) {return('poisson')}
                     if ((length(support) == 2) & all(support == c(0,1))) {return('bernoulli')}
                     if ((length(support) == 2) & any(support != c(0,1))) {return('poisson')}
                     if ((length(support) < 2) & (support == 0 | support == 1)) {return('bernoulli')}
                   }else{return('continuous')}
                } ))

                w.continuous <- which(v_distrib_guessed == 'continuous')
                w.noncontinuous <- (1:self$cardE)[-w.continuous]

                if (length(w.continuous) > 0) {
                  check <- 1
                  for (u in w.continuous) {check = check * as.numeric(v_distrib[u] %in% c('gaussian','laplace'))}
                  if ( check == 0) {stop('Check distribution for continuous weighted network')}
                }
                if (length(w.noncontinuous) > 0) {
                  if (any(v_distrib_guessed[w.noncontinuous] != self$v_distrib[w.noncontinuous])) {stop('Non adequate distributions')}
                }



                # creating Ecode
                self$Ecode = transfoE(self$E, self$typeInter)


                prov = private$check()
                self$v_NQ = prov$n_q
                self$where <- prov$where_q
                if (is.null(namesInd)) { self$namesInd <- lapply(1:self$Q, function(q){
                  where_q <- self$where[[q]][1,];
                  if (where_q[2] == 1) {namesInd_q = rownames(self$mats[[where_q[1]]])}else{namesInd_q = colnames(self$mats[[where_q[1]]])}
                  return(namesInd_q)})}
                else{self$namesInd <-  lapply(1:self$Q,function(q){namesInd[[q]]})}

                },
            estime = function(classif,tau=NULL){varEMMBM(self,classif,tau)},
            cleanResults = function(R){cleanEstim(self,R)},
            searchNbClusters = function(classifInit,Kmin,Kmax, nbCores = NULL, verbose = TRUE)
            {
              searchKQ(dataR6 = self,classifInit = classifInit,Kmin = Kmin,Kmax = Kmax,nbCores = nbCores,verbose = verbose)
            }
        ),
      private = list(
                 #version de E interne au code
                  check = function()
                  {checkExtract(self$mats,self$Ecode)
                    }
                )

)

#------------------------------- CLASS Multipartite Block Models
MBMfit = R6Class("MBMfit",
              public = list(
                v_K = NULL, #vector of number of blocks in each functional group
                v_distrib = NULL, #vector of emission distribution (same length as number mats) (bernoulli, poisson, gaussian, laplace...)
                list_pi = NULL, #list of vectors (length given in vK) for mixture distribution of Z
                list_theta = NULL,
                v_NQ = NULL, # number of individuas by functional groups.
                E = NULL,
                typeInter = NULL,
                Q = NULL,
                Z = NULL,
                tau = NULL,
                initialize = function(v_K,v_distrib,list_pi = NULL, list_theta = NULL)
                  {
                    self$v_K <- v_K
                    self$v_distrib <- v_distrib
                    self$list_pi <- list_pi
                    self$list_theta <- list_theta

                  }
                )
                 )

MBMfit$set("public",'sim',
             function(seed = NULL, E, v_NQ, typeInter,keepClassif = FALSE){
               self$E <- E;
               self$v_NQ <- v_NQ;
               self$typeInter <- typeInter;
               self$Q <- length(unique(c(E)))
               if (!is.null(seed)) {set.seed(seed)}
               Z <- lapply(1:self$Q,function(q){Zq <- sample(1:self$v_K[q],self$v_NQ[q],replace = TRUE,prob = self$list_pi[[q]])})
               self$Z = Z
               mats <- lapply(1:nrow(self$E),function(e){
                 fg1 <- self$E[e,1]
                 fg2 <- self$E[e,2]
                 list_theta_e <- list_theta[[e]] ###
                 Z_fg1 <- self$Z[[fg1]]
                 Z_fg2 <- self$Z[[fg2]]
                 switch(self$v_distrib[e],
                        bernoulli = {
                          X_e <- matrix(rbinom(self$v_NQ[fg1] * self$v_NQ[fg2],1,list_theta_e[Z_fg1,Z_fg2]),self$v_NQ[fg1],self$v_NQ[fg2])
                          diag(X_e) <- 0;                              },
                        poisson = {
                          X_e <- matrix(rpois(self$v_NQ[fg1] * self$v_NQ[fg2],list_theta_e[Z_fg1,Z_fg2]),self$v_NQ[fg1],self$v_NQ[fg2])
                          },
                        gaussian  = {
                          X_e <- matrix(rnorm(self$v_NQ[fg1] * self$v_NQ[fg2],mean = list_theta_e$mean[Z_fg1,Z_fg2], sd = list_theta_e$sd[Z_fg1,Z_fg2]),self$v_NQ[fg1],self$v_NQ[fg2])
                          },
                       laplace = {
                          X_e <- matrix(rlaplace(self$v_NQ[fg1] * self$v_NQ[fg2], location = 0, scale = list_theta_e[Z_fg1,Z_fg2]),self$v_NQ[fg1],self$v_NQ[fg2])
                          },
                      stop("Enter a valid distribution (poisson or bernoulli or laplace or gaussian)!"))
                 if (self$typeInter[e] == "adj") {X_e[lower.tri(X_e)] = t(X_e)[lower.tri(X_e)]}
                 return(X_e)}
               )
               dataSim <- CollInteraction$new(mats = mats,self$E,self$typeInter,self$v_distrib)
               res <- list(networks = dataSim)
               if (keepClassif) {res$classif <- Z}
               return(res)}

)



