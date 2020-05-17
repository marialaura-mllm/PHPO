print.summary.weibull_phpo <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  if(x$approach=="mle"){
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Coefficients:\n")
    printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)


  }else{
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Coefficients:\n")
    print(x$coefficients)
    cat("--- \n")
    cat("Inference for Stan model: ", x$model_name, ';\n', sep = '')
    cat(x$chains, " chains, each with iter=", x$iter,
      "; warmup=", x$warmup, "; thin=", x$thin, "; \n",
      "post-warmup draws per chain=", x$n_kept[1], ", ",
      "total post-warmup draws=", sum(x$n_kept), ".\n\n", sep = '')

  }
}


summary.weibull_phpo <- function(object, ...){

  if(object$approach=="mle"){
    p <- object$p
    labels <- object$labels
    coefficients <- object$fit$par[3:(p+2)]
    vcov <- solve(-object$fit$hessian[3:(p+2), 3:(p+2)])

    resid <- object$residuals
    se <- sqrt(diag(vcov))
    zval <- coefficients / se
    TAB <- cbind(coef = coefficients,
      "se(coef)" = se,
      z = zval,
      p.value = 2*stats::pnorm(-abs(zval)))
    rownames(TAB) <- labels[1:length(coefficients)]

    res <- list(call=object$call,
      coefficients=TAB, residuals=resid, approach=object$approach)
    class(res) <- "summary.weibull_phpo"
    res

    # Bayesiam output:
  }else{
    labels <- object$labels
    s <- rstan::summary(object$fit, pars=c("beta"))
    TAB <- round(s$summary, digits = 3)
    rownames(TAB) <- labels
    n_kept <- object$fit@sim$n_save - object$fit@sim$warmup2

    res <- list(call=object$call, coefficients=TAB, n_kept=n_kept,
      model_name=object$fit@model_name, chains=object$fit@sim$chains,
      warmup=object$fit@sim$warmup, thin=object$fit@sim$thin,
      iter=object$fit@sim$iter, approach=object$approach)

  }

  class(res) <- "summary.weibull_phpo"

  return(res)
}
