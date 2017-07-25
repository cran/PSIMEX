initial <-
function(model, data, pedigree) {
 
  if (!is.element("f_inb", all.vars(model$call))){
  data$f_inb <- calcInbreeding(pedigree)
  
  vars <- all.vars(model$call$formula)[-1]
  vars <- c(vars, "f_inb")
  covars <- paste(vars, collapse = "+")
  
  y <- all.vars(model$call)[1]
  fla <- paste(y,"~", covars)
  
  model$call$formula <- fla
  model <- update(model)
  }
  
  
  
  init <- list()
  init$model <- model
  init$inb0 <- summary(model)$coefficient["f_inb", 1]
  init$se0 <- summary(model)$coefficient["f_inb", 2]
  init$pval0 <- summary(model)$coefficient["f_inb", 4]
  init$mean_inb0 <- mean(data$f_inb)
  init$median_inb0 <- median(data$f_inb)
  init$var_inb0 <- var(data$f_inb)
  
  return(init)
  
}
