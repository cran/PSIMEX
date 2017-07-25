simul.na.bottom <-
function (pedigree, pedigree0, lambda, B, data, model, nbottom) {
  
  
  names(pedigree) <- c("animal", "dam", "sire")
  
  
  keep <- pedigree0[which(pedigree0$generation<nbottom, ), 1]
  keep_ped <- pedigree[keep, ]
  pedigree <- pedigree[-keep, ]
  names(keep_ped) <- c("animal", "dam", "sire")
  
  
  
  #where to generate more errors
  known <- pedigree[which(!is.na(pedigree$dam)), ]
  #already missing
  unknown <- pedigree[which(is.na(pedigree$dam)), ]
  names(known) <- c("animal", "dam", "sire")
  names(unknown) <- c("animal", "dam", "sire")
  
  #setting initial values
  inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol  =  B)
  se_inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol  =  B)
  mean_inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol = B)
  median_inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol = B)
  var_inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol = B)
  pval <- matrix(data  =  NA, nrow  =  length(lambda), ncol = B)
  
  for (k in 1:length(lambda)){
    for (j in 1:B) {
      animal <- known$animal
      dam <- known$dam
      sire <- known$sire  
      
      #initializing inbreeding coefficient
      data$f_inb<- c()
      
      #choosing lambda over the total number of individuals 
      #and setting their parents as random
      m <- sample(1:length(known$animal), lambda[k]*length(known$animal))
      for (i in m) {
        dam[i]<-NA
        sire[i]<-NA
        
      }
      
      #recalculate pedigree
      ped <- data.frame(animal = animal, dam = dam, sire = sire)
      ped <- rbind(ped, unknown, keep_ped)
      ord <- orderPed(ped)
      ped <- ped[order(ord),]
      
      #re calculate inbreeding coefficient
      f_inb <- data.frame(ped$animal, calcInbreeding(ped))
      names(f_inb) <- c("id", "f_inb")
      data <- merge(data, f_inb, by = "id")
      data$f_inb <- as.numeric(data$f_inb)
      
      if(sum(data$f_inb)>0){
        
        #linear model 
        #inbreeding depression is the regression coefficient of f_inb
        type <- as.character(model$call[1])
        if (type == 'lm'){
          model1 <- lm(formula = as.character(model$call[2]), data = data)
        }
        
        if (type == 'glm'){
          model1 <- glm(formula = as.character(model$call[2]), family = model$call$family,  data = data)
        }
        
        
        #inbreeding depression
        inb[k, j] <- summary(model1)$coefficients["f_inb", 1]
        se_inb[k,j] <- summary(model1)$coefficients["f_inb", 2]
        pval [k, j] <- summary(model1)$coefficients["f_inb", 4]
        
        #info on inbreeding coefficient
        mean_inb[k, j] <- mean(data$f_inb)
        median_inb[k, j] <- median(data$f_inb)
        var_inb[k, j] <- var(data$f_inb)
        
      }
      #this must be checked
      if(sum(data$f_inb) == 0) {
        inb[k, j] <- NA
        se_inb[k,j] <- NA
        mean_inb[k, j] <- NA
        median_inb[k, j] <- NA
        var_inb[k, j] <- NA
        pval [k, j] <- NA
        
      }
      
    }
    
  }
  
  #store results
  res = list(inb, se_inb, pval, 
           mean_inb, median_inb, var_inb)
  names(res) <- c('inb_dep', 'se_inb_dep', 'pval_inb_dep', 'mean_inb', 'median_inb', 'var_inb')
  return(res)
}
