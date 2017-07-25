simul.na <-
function(pedigree, pedigree0, lambda, lambda0,  B, data, model, parameter, way, ntop, nbottom, prior, nitt,thin,burnin){
  lambda <- lambda - lambda0
  
  
   if (way =='uniform' & parameter == 'inbreeding')    results <- simul.na.uni(pedigree, lambda, 
                                                                              B, data, model)
   if (way =='uniform' & parameter == 'heritability')  results <- simul.na.uni.herit(pedigree, lambda,  
                                                                                    B, data, model, 
                                                                                    prior,nitt,thin,burnin)

                                                                                  
   if (way == 'top' & parameter == 'inbreeding')       results <- simul.na.top(pedigree, pedigree0, lambda, 
                                                                               B, data, model, ntop)
   if (way == 'top' & parameter == 'heritability')     results <- simul.na.top.herit(pedigree, pedigree0, lambda, 
                                                                                     B,  data, model, ntop, 
                                                                                     prior,nitt,thin,burnin)
        
    if (way == 'bottom' & parameter == 'inbreeding')    results <- simul.na.bottom(pedigree,pedigree0, lambda,
                                                                                   B,  data, model, nbottom)
   if (way == 'bottom' & parameter == 'heritability')   results <- simul.na.bottom.herit(pedigree, pedigree0, lambda, 
                                                                                        B,  data, model, nbottom,  
                                                                                        prior, nitt,thin,burnin)
         
  
  

  
  return(results)
  
  
  
  
  
}
