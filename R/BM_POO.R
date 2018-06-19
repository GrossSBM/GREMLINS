coll_interaction=R6::R6Class("coll_interaction", ### classe objet pour décrire les données
      public=list(
            mats=NULL,    #list of mats of interaction
            E=NULL, # 2-column table giving which fgs (functional groups) interact in corresponding mat
            namesfg=NULL, #vector of names of fgs
            names_ind = NULL,
            type_inter=NULL, #vector of type of matrices (diradj=directed adjacency, adj, inc=incidence)
            vdistrib=NULL, #vector of emission distribution (same length as number mats) (poisson, bernoulli...)
            Q=NULL,     #number of functional groups
            card_E=NULL, #number of mats
            v_NQ=NULL, #number of individuals in fgs
            where=NULL, #useful
            Ecode=NULL,
            initialize=function(MATS,E_FG,type,distrib=NULL,names_ind=NULL)
              {
                self$card_E=length(MATS)
                self$mats=MATS
                #if (is.null(namesFG)){self$namesfg=unique(as.vector(E_FG))}else {self$namesfg=namesFG}
                self$namesfg = unique(as.vector(E_FG))# les fonctional groups names sont toujours ceux dans l'ordre d'apparition dans les matrices

                self$Q=length(self$namesfg)
                self$type_inter=type
                self$E=matrix(sapply(E_FG,function(a){which(self$namesfg==a)}),self$card_E,2,byrow=FALSE)

                #init vdistrib, default Bernoulli
                if (is.null(distrib)) self$vdistrib=rep("bernoulli",self$card_E)
                else
                  {
                    if (length(distrib)!=self$card_E) stop("number of distributions not consistent with number of interaction matrices")
                    self$vdistrib=distrib
                  }

                # creating Ecode
                self$Ecode = transfo_E(self$E, self$type_inter)


                prov = private$check()
                self$v_NQ = prov$n_q
                self$where <- prov$where_q
                if (is.null(names_ind)) { self$names_ind <- lapply(1:self$Q, function(q){
                  where_q <- self$where[[q]][1,];
                  if(where_q[2] == 1) {names_ind_q = rownames(self$mats[[where_q[1]]])}else{names_ind_q = colnames(self$mats[[where_q[1]]])}
                  return(names_ind_q)})}
                else{self$names_ind <-  lapply(1:self$Q,function(q){names_ind[[q]]})}

                },
            estime = function(classif){VEM_gen_BM(self,classif)},
            search_nb_clusters = function(vKinit,Kmin,Kmax,nb_cores=NULL,verbose=TRUE){search_KQ(data = self,vKinit = vKinit,Kmin = Kmin,Kmax = Kmax,nb_cores = nb_cores,verbose = verbose)}),
      private=list(
                 #version de E interne au code
                  check=function()
                  {check_extract(self$mats,self$Ecode)
                    }
                )

)

genBMfit=R6::R6Class("genBMfit",
              public=list(
                vK=NULL, #vector of number of blocks in each functional group
                vdistrib=NULL, #vector of emission distribution (same length as number mats) (poisson, bernoulli...)
                lpi=NULL, #list of vectors (length given in vK) for mixture distribution of Z
                ltheta = NULL,
                v_NQ=NULL, # number of individuas by functional groups.
                E=NULL,
                type_inter=NULL,
                Q=NULL,
                Z=NULL,
                tau = NULL,
                initialize=function(vK,vdistrib,lpi=NULL,ltheta=NULL)
                  {
                    self$vK <- vK
                    self$vdistrib <- vdistrib
                    self$lpi <- lpi
                    self$ltheta <- ltheta

                  }
                )
                 )

genBMfit$set("public",'sim',
             function(seed=NULL,E,v_NQ,type_inter){
               self$E <- E;
               self$v_NQ <- v_NQ;
               self$type_inter <- type_inter;
               self$Q<- length(unique(c(E)))
               if (!is.null(seed)){set.seed(seed)}
               Z <- lapply(1:self$Q,function(q){Zq <- sample(1:self$vK[q],self$v_NQ[q],replace=TRUE,prob=self$lpi[[q]])})
               self$Z = Z
               mats <-lapply(1:nrow(self$E),function(e){
                 fg1 <- self$E[e,1]
                 fg2 <- self$E[e,2]
                 ltheta_e <- ltheta[[e]] ###
                 Z_fg1 <- self$Z[[fg1]]
                 Z_fg2 <- self$Z[[fg2]]
                 switch(self$vdistrib[e],
                        bernoulli = {
                          X_e <- matrix(rbinom(self$v_NQ[fg1]*self$v_NQ[fg2],1,ltheta_e[Z_fg1,Z_fg2]),self$v_NQ[fg1],self$v_NQ[fg2])
                          diag(X_e) <- 0;                              },
                        poisson = {X_e <-matrix(rpois(self$v_NQ[fg1]*self$v_NQ[fg2],ltheta_e[Z_fg1,Z_fg2]),self$v_NQ[fg1],self$v_NQ[fg2])},
                        stop("Enter a valid distribution (poisson or bernoulli)!"))
                 if (self$type_inter[e]=="adj"){X_e [lower.tri(X_e)] = t(X_e)[lower.tri(X_e)]}
                 return(X_e)}
               )
               res <- coll_interaction$new(MATS=mats,self$E,self$type_inter,self$vdistrib)
               return(res)}

)


