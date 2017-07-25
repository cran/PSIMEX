simul.replace.uni.herit <-
function (pedigree, lambda, B,  data, model, prior, nitt,thin,burnin) {
 
  
  
  #where to generate more errors
  known <- pedigree[which(!is.na(pedigree$dam)), ]
  #already missing
  unknown <- pedigree[which(is.na(pedigree$dam)), ]
    
  names(known) <- c("animal", "sire", "dam")
  names(unknown) <- c("animal", "sire", "dam")
  
  
  #setting initial values
  h <- matrix(data=NA, nrow=length(lambda), ncol=B)
  se_h <- matrix(data=NA, nrow=length(lambda), ncol=B)
  VA <- matrix(data=NA, nrow=length(lambda), ncol=B)
  VE <- matrix(data=NA, nrow=length(lambda), ncol=B)
  
  for (k in 1:length(lambda)){
    for (j in 1:B) {
      animal <- known$animal
      dam <- known$dam
      sire <- known$sire  
      
     
      
      #choosing lambda over the total number of individuals 
      #and setting their parents as random
      m <- sample(1:length(known$animal), lambda[k]*length(known$animal))
      for (i in m) {
        
        generation <- data[which(data$id == known$sire[i]), ]$year
        others <- data[which(data$year == generation & data$sex == 1),  ]$id
        others <- others[-which(others == known[i, 2])]
        sire[i] <- others[sample(1:length(others), 1)]
        
      }
      
      
      
      #recalculate pedigree
      ped <- data.frame(animal = animal, dam = dam, sire = sire)
      ped <- rbind(ped, unknown[c('animal', 'dam', 'sire')])
      ord <- orderPed(ped)
      ped <- ped[order(ord),]
      
      #recalculate heritability
      model <- MCMCglmm(model$Fixed$formula,random=model$Random$formula,
                        pedigree=ped,data=data,prior=prior,
                        nitt=nitt,thin=thin,burnin=burnin,verbose=F)
      
      
      
      
        random <- all.vars(model$Random$formula)
        if (length(random) > 1)     VEs <- rowSums(model$VCV[ ,random])+model$VCV[ , 'units']
        if (length(random) == 1)     VEs <- model$VCV[ ,random]+model$VCV[ , 'units']
      
      
      
      h_expr <- model$VCV[,"animal"]/VEs
      h[k, j] <- mean(h_expr)
      
      se_h[k, j] <- sd(h_expr)
      
      VA[k,j]<-mean(model$VCV[,"animal"])
      
      VE[k,j]<- mean(VEs)
      
      
    }
    
  }
  
  #store results
  res = list(h, se_h, VA,VE)
  names(res) <- c('heritability', 'se_h', 'VA', 'VE')
  return(res)
}
