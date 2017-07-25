plot.simul.Psimex <-
function(results, lambda, lambda0, estimate0, parameter){
  
  
  
  #merge all values
  lambda<-c(lambda0, lambda)
  if (parameter == 'inbreeding'){
    inb <- t(results$inb_dep)
    se_inb <- t(results$se_inb_dep)
    
    #compute means
    inb_ave<-colMeans(inb, na.rm = TRUE)
    se_inb_ave<-colMeans(se_inb, na.rm = TRUE)
    
    
  par_ave<-c(estimate0$inb0, inb_ave)
  se<-c(estimate0$se0, se_inb_ave)
  }
  
  if (parameter == 'heritability'){
    h <- t(results$h)
    se_h <- t(results$se_h)
    VA <- t(results$VA)
    VE <- t(results$VE)
    
    
    #compute means
    h_ave<-colMeans(h, na.rm = TRUE)
    se_h_ave<-colMeans(se_h, na.rm = TRUE)
    VA_ave<-colMeans(VA, na.rm = TRUE)
    VE_ave<-colMeans(VE, na.rm = TRUE)
    
    
    par_ave<-c(estimate0$h0, h_ave)
    se<-c(estimate0$seh0, se_h_ave)
  }
  
  
  
  low <- par_ave-1.96*se
  up <- par_ave+1.96*se
  
  plotCI(lambda, par_ave, ui=up, li=low, xlab="Error proportion", 
         ylab="Parameter", col="black", cex=1.5, main="Simulation results", cex.lab=1.3, cex.main=1.3)
  
  
}
