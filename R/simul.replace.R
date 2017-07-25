simul.replace <-
function(pedigree, lambda, lambda0, B, data, model, parameter,  way, prior, nitt,thin,burnin){
  lambda <- lambda - lambda0 
  if (parameter == 'inbreeding' & way == 'uniform')  results <- simul.replace.uni(pedigree, lambda, 
                                                                                  B,  data, model)
  if (parameter == 'inbreeding' & way == 'similar') results <- simul.replace.similar(pedigree, lambda, 
                                                                                     B,   data, model)
  
  if (parameter == 'heritability' & way == 'uniform')  results <- simul.replace.uni.herit(pedigree, lambda, B,
                                                                                         data, model, prior = prior, 
                                                                                         nitt, thin, burnin)
  
  if (parameter == 'heritability' & way == 'similar')  results <- simul.replace.similar.herit(pedigree, lambda, B,   
                                                             data, model, prior, 
                                                             nitt,thin,burnin)
  
  
  
  return(results)
  
  
  
  
  
}
