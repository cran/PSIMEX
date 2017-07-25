plot.Psimex <-
function(results, extrapolation_results, lambda, lambda0, estimate0, parameter, fitting.method){
  #merge all values
  lambda<-c(lambda0, lambda)
  if (parameter == 'inbreeding'){
    
    inb <- t(results$inb_dep)
    se_inb <- t(results$se_inb_dep)
    
    #compute means
    inb_ave<-colMeans(inb, na.rm = TRUE)
    se_inb_ave<-colMeans(se_inb, na.rm = TRUE)
    
    
    
  inb_ave<-c(estimate0$inb0, inb_ave)
  inb_pred <- extrapolation_results$inb_dep
  
  p.names<-c("inb_ave", "lambda")
  estimates<-data.frame(inb_ave, lambda)
  colnames(estimates) <- p.names

  
  plot(inb_ave[-1] ~lambda[-1], xlim = c(0, max(lambda)), ylim = c(min(inb_pred, min(inb_ave)), max(inb_pred, max(inb_ave))), cex = 1.5,  
       xlab = 'Error proportion', ylab = 'Inbreeding depression', main = "Extrapolation phase")
  for ( i in 1: length(inb_pred)){
  points(i/1000, inb_pred[i], pch = 4+i, cex = 2, col = 'lightblue')
  if (fitting.method[i] == 'line') mod <- lm(inb_ave~lambda, data = estimates)
  if (fitting.method[i] == 'quad') mod <- lm(inb_ave~lambda+I(lambda^2), data = estimates)
  if (fitting.method[i] == 'cubi') mod <- lm(inb_ave~lambda+I(lambda^2)+I(lambda^3), data = estimates)
    lambda2 <- c(0, lambda)
    lines(lambda2, predict(mod, data.frame(lambda = lambda2)), 
          col="lightblue", lty=2)  
  }
  points(lambda[1], inb_ave[1], pch = 18, cex = 2)
  
  text <- c()
  for (  i in 1: length(inb_pred)){
    text <- c(text, paste('SIMEX estimate: ', fitting.method[i]))
  }
  legend('bottomright',   inset = 0,  c(text,  "Naive Estimate"), 
         col = c(rep("lightblue", length(inb_pred)),"black"),  pch = c(5:(4+length(inb_pred)), 18))
  }
  
  if (parameter == 'heritability'){
    h <- t(results$heritability)
    se_h <- t(results$se_h)
    VA <- t(results$VA)
    VE <- t(results$VE)
    
    
    #compute means
    h_ave<-colMeans(h, na.rm = TRUE)
    se_h_ave<-colMeans(se_h, na.rm = TRUE)
    VA_ave<-colMeans(VA, na.rm = TRUE)
    VE_ave<-colMeans(VE, na.rm = TRUE)
    
    h_ave<-c(estimate0$h0, h_ave)
    h_pred <- extrapolation_results$h
    
    p.names<-c("h_ave", 'lambda')
    estimates<-data.frame(h_ave, lambda)
    colnames(estimates) <- p.names
    
    
    
    plot(h_ave[-1] ~lambda[-1], xlim = c(0, max(lambda)),
         ylim = c(min(h_pred, min(h_ave)), max(h_pred, max(h_ave))),
         cex = 1.5,  xlab = 'Error proportion', ylab = 'Heritability', main = "Extrapolation phase")
    for ( i in 1: length(h_pred)){
      points(i/1000, h_pred[i], pch = 4+i, cex = 2, col = 'lightblue')
      if (fitting.method[i] == 'line') mod <- lm(h_ave~lambda, data = estimates)
      if (fitting.method[i] == 'quad') mod <- lm(h_ave~lambda+I(lambda^2), data = estimates)
      if (fitting.method[i] == 'cubi') mod <- lm(h_ave~lambda+I(lambda^2)+I(lambda^3), data = estimates)
      lambda2 <- c(0, lambda)
      lines(lambda2, predict(mod, data.frame(lambda = lambda2)), 
            col="lightblue", lty=2)  
      }
    
    points(lambda[1], h_ave[1], pch = 18, cex = 2)
    
    text <- c()
    for (  i in 1: length(h_pred)){
      text <- c(text, paste('SIMEX estimate: ', fitting.method[i]))
    }
    legend('bottomright',   inset = 0,  c(text,  "Naive Estimate"), 
           col = c(rep("lightblue", length(h_pred)),"black"),  pch = c(5:(4+length(h_pred)), 18))
  }
  
}
