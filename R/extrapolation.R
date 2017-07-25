extrapolation <-
function (results, lambda, lambda0, estimate0, fitting.method, B, parameter) {
  
  
  if (parameter == 'inbreeding'){
    
    inb <- t(results$inb_dep)
    se_inb <- t(results$se_inb_dep)
    
    #compute means
    inb_ave<-colMeans(inb, na.rm = TRUE)
    se_inb_ave<-colMeans(se_inb, na.rm = TRUE)
    
  
  #merge all values
  lambda<-c(lambda0, lambda)
  inb_ave<-c(estimate0$inb0, inb_ave)
  se_inb<-c(estimate0$se0, se_inb_ave)

  
  #extrapolation for inbreeding depression
  p.names<-c("inb", "se_inb")
  estimates<-data.frame(inb_ave, se_inb)
  colnames(estimates) <- p.names
  
  inb_pred <- c()
  inb_pred_se <- c()
  AIC <- c()
  se_pred <- c()
  var <- c()
  
  for (i in 1: length(fitting.method)) {
    f <- fitting.method[i]
  #linear case 
  if (f == 'line') {

    extrapolation_inb <- lm(estimates[ ,1] ~ lambda)
    inb_pred[i] <- predict(extrapolation_inb, newdata  =  data.frame(lambda  =  0))
    extrapolation_inb_se <- lm(estimates[,2] ~ lambda)
    inb_pred_se[i] <- predict(extrapolation_inb_se, newdata  =  data.frame(lambda  =  0))
    AIC[i] <- AIC(extrapolation_inb)
    }   
  
  
  #quadratic case
  if (f == 'quad') {

    extrapolation_inb1<- lm(estimates[ ,1] ~ lambda+ I(lambda^2))
    inb_pred[i] <- predict(extrapolation_inb1, newdata  =  data.frame(lambda  =  0))
    extrapolation_inb1_se<- lm(estimates[ ,2] ~ lambda+ I(lambda^2))
    inb_pred_se[i] <-predict(extrapolation_inb1_se, newdata  =  data.frame(lambda  =  0))
    AIC[i] <- AIC(extrapolation_inb1)
    } 
    
  #non linear case 
  if (f == 'nonl') {

    extrapolation_inb2<-fit.nls(lambda, p.names, estimates[, ])
    inb_pred[i] <-predict(extrapolation_inb2[[1]], newdata  =  data.frame(lambda  =  0))
    inb_pred_se[i] <-predict(extrapolation_inb2[[2]], newdata  =  data.frame(lambda  =  0))
    AIC[i] <- AIC(extrapolation_inb2)
   }
  
  #cubic case
  if (f == 'cubi') {
    
    extrapolation_inb<- lm(estimates[ ,1] ~ lambda+ I(lambda^2)+ I(lambda^3))
    inb_pred[i] <-predict(extrapolation_inb, newdata  =  data.frame(lambda  =  0))
    extrapolation_inb_se<- lm(estimates[,2] ~ lambda+ I(lambda^2)+ I(lambda^3))
    inb_pred_se[i] <-predict(extrapolation_inb_se, newdata  =  data.frame(lambda  =  0))
    AIC[i] <- AIC(extrapolation_inb)
    
  }   
  
  #sampling error across simulation (JACKKNIFE)
  S <- c()

  for ( ii in 1:length(inb[ 1,])) {
    diff <- c()
    #calculate the differences per each simulation
    for ( j in 1: B ) {
      diff[j] <- inb[j,ii]-inb_ave[ii]
      
      
    }
    diff <- na.omit(diff)
    S[ii] <- 1/(B-1)*sum(diff^2) 
    
  }
  
  #extrapolated value for sampling error
  lambda1 <- lambda[-1]
  #linear case
  if (f == 'line') {
    
  extrapolation_var <- lm(S ~ lambda1)
  se_pred[i] <- predict(extrapolation_var, newdata = data.frame(lambda1 = 0))
  }
  
  #quadratic case
  
  if (f == 'quad') {
    
  extrapolation_var1 <- lm(S ~ lambda1+ I(lambda1^2))
  se_pred[i] <- predict(extrapolation_var1, newdata = data.frame(lambda1 = 0))
  }
  
  # non linear case
  
  if (f == 'nonl') {
  p.names1 <- c("S")
  estimates1 <- data.frame(S)
  colnames(estimates1) <- p.names1
  ex <- fit.nls(lambda1, p.names1, estimates1[ ])
  
  se_pred[i] <- predict(ex[[1]], newdata = data.frame(lambda1 = 0))
  }
  #cubic case
  if (f == 'cubi') {
    
    extrapolation_var1 <- lm(S ~ lambda1+ I(lambda1^2)+ I(lambda1^3))
    se_pred[i] <- predict(extrapolation_var1, newdata = data.frame(lambda1 = 0))
  }
  
  #comprehensive measure of error
  var[i] <- inb_pred_se[i]^2 + se_pred[i]^2
  }
  
  extrapolation <- list()
  extrapolation$fitting.method <- fitting.method
  extrapolation$inb_dep <- inb_pred
  extrapolation$se <- inb_pred_se
  extrapolation$sampling_se <- se_pred
  extrapolation$var <- var
  extrapolation$AIC <- AIC
 
  }
      
  
   if (parameter =='heritability'){
     h <- t(results$heritability)
     se_h <- t(results$se_h)
     VA <- t(results$VA)
     VE <- t(results$VE)
     
     
     #compute means
     h_ave<-colMeans(h, na.rm = TRUE)
     se_h_ave<-colMeans(se_h, na.rm = TRUE)
     VA_ave<-colMeans(VA, na.rm = TRUE)
     VE_ave<-colMeans(VE, na.rm = TRUE)
     
     #merge all values
     lambda <- c(lambda0, lambda)
     h_ave <- c(estimate0$h0, h_ave)
     se_h <- c(estimate0$seh0, se_h_ave)
     VA <- c(estimate0$VA0, VA_ave)
     VE <- c(estimate0$VE0, VE_ave)
     
     #extrapolation for inbreeding depression
     p.names<-c("h", "se_h", "VA", "VE")
     estimates<-data.frame(h_ave, se_h, VA, VE)
     colnames(estimates) <- p.names
     
     
     h_pred <- c()
     h_pred_se <- c()
     VA_pred <- c()
     VE_pred <- c()
     AIC <- c()
     aicc <- c()
     se_pred <- c()
     var <- c()
     
     for (i in 1: length(fitting.method)) {
       f <- fitting.method[i]
     #linear case 
     if (f == 'line') {
       
       extrapolation_h<- lm(estimates[ ,1] ~ lambda)
       h_pred[i] <-predict(extrapolation_h, newdata  =  data.frame(lambda  =  0))
       extrapolation_h_se <- lm(estimates[,2] ~ lambda)
       h_pred_se[i] <- predict(extrapolation_h_se, newdata  =  data.frame(lambda  =  0))
       extrapolation_VA <- lm(estimates[ ,3] ~ lambda)
       VA_pred[i] <-predict(extrapolation_VA, newdata  =  data.frame(lambda  =  0))
       extrapolation_VE<- lm(estimates[ ,4] ~ lambda)
       VE_pred[i] <-predict(extrapolation_VE, newdata  =  data.frame(lambda  =  0))
       AIC[i] <- AIC(extrapolation_h)
     }   
     
     
     #quadratic case
     if (f == 'quad') {
       
       extrapolation_h<- lm(estimates[ ,1] ~ lambda+ I(lambda^2))
       h_pred[i] <-predict(extrapolation_h, newdata  =  data.frame(lambda  =  0))
       extrapolation_h_se<- lm(estimates[ ,2] ~ lambda+ I(lambda^2))
       h_pred_se[i] <-predict(extrapolation_h_se, newdata  =  data.frame(lambda  =  0))
       extrapolation_VA <- lm(estimates[ ,3] ~ lambda+ I(lambda^2))
       VA_pred[i] <-predict(extrapolation_VA, newdata  =  data.frame(lambda  =  0))
       extrapolation_VE<- lm(estimates[ ,4] ~ lambda+ I(lambda^2))
       VE_pred[i] <-predict(extrapolation_VE, newdata  =  data.frame(lambda  =  0))
       AIC[i] <- AIC(extrapolation_h)
     } 
     
     #non linear case 
     if (f == 'nonl') {
       
       extrapolation_h<-fit.nls(lambda, p.names, estimates[, ])
       h_pred[i] <-predict(extrapolation_h[[1]], newdata  =  data.frame(lambda  =  0))
       h_pred_se[i] <-predict(extrapolation_h[[2]], newdata  =  data.frame(lambda  =  0))
       VA_pred[i] <-predict(extrapolation_h[[3]], newdata  =  data.frame(lambda  =  0))
       VE_pred[i] <-predict(extrapolation_h[[4]], newdata  =  data.frame(lambda  =  0))
       
     }
     #cubic case
     if (f == 'cubi') {
       
       extrapolation_h<- lm(estimates[ ,1] ~ lambda+ I(lambda^2)+I(lambda^3))
       h_pred[i] <-predict(extrapolation_h, newdata  =  data.frame(lambda  =  0))
       extrapolation_h_se<- lm(estimates[ ,2] ~ lambda+ I(lambda^2)+I(lambda^3))
       h_pred_se[i] <-predict(extrapolation_h_se, newdata  =  data.frame(lambda  =  0))
       extrapolation_VA <- lm(estimates[ ,3] ~ lambda+ I(lambda^2)+I(lambda^3))
       VA_pred[i] <-predict(extrapolation_VA, newdata  =  data.frame(lambda  =  0))
       extrapolation_VE<- lm(estimates[ ,4] ~ lambda+ I(lambda^2)+I(lambda^3))
       VE_pred[i] <-predict(extrapolation_VE, newdata  =  data.frame(lambda  =  0))
       AIC[i] <- AIC(extrapolation_h)
     } 
     
     
     #sampling error across simulation (JACKKNIFE)
     S <- c()
     
     for ( ii in 1:length(h[ 1, ])) {
       diff <- c()
       #calculate the differences per each simulation
       for ( j in 1: B ) {
         diff[j] <- h[j,ii]-h_ave[ii]
         
         
       }
       diff <- na.omit(diff)
       S[ii] <- 1/(B-1)*sum(diff^2) 
       
     }
     
     #extrapolated value for sampling error
     lambda1 <- lambda[-1]
     #linear case
     if (f == 'line') {
       
       extrapolation_var <- lm(S ~ lambda1)
       se_pred[i] <- predict(extrapolation_var, newdata = data.frame(lambda1 = 0))
     }
     
     #quadratic case
     
     if (f == 'quad') {
       
       extrapolation_var1 <- lm(S ~ lambda1+ I(lambda1^2))
       se_pred[i] <- predict(extrapolation_var1, newdata = data.frame(lambda1 = 0))
     }
     
     # non linear case
     
     if (f == 'nonl') {
       p.names1 <- c("S")
       estimates1 <- data.frame(S)
       colnames(estimates1) <- p.names1
       ex <- fit.nls(lambda1, p.names1, estimates1[ ])
       
       se_pred[i] <- predict(ex[[1]], newdata = data.frame(lambda1 = 0))
     }
     #cubic case
     
     if (f == 'cubi') {
       
       extrapolation_var1 <- lm(S ~ lambda1+ I(lambda1^2)+I(lambda1^3))
       se_pred[i] <- predict(extrapolation_var1, newdata = data.frame(lambda1 = 0))
     }
     #comprehensive measure of error
     var[i] <- h_pred_se[i]^2 + se_pred[i]^2
     }
     
     extrapolation <- list()
     extrapolation$fitting.method <- fitting.method
     extrapolation$h <- h_pred
     extrapolation$se <- h_pred_se
     extrapolation$sampling_se <- se_pred
     extrapolation$var <- var 
     extrapolation$VA <- VA_pred
     extrapolation$VE <- VE_pred
     extrapolation$AIC <- AIC
 
   }
  return(extrapolation)
}
