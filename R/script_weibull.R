mle_weibull_PH <- function(time, status, x, ...){

  p <- ncol(x)
  n <- nrow(x)
  tau <- max(time)
  time2 <- time/tau

  a_gamma <- 1
  b_gamma <- 1
  a_lambda <- 0.1
  b_lambda <- 0.1
  mu_beta <- 0
  sigma_beta <- 1

  data_mle_PH <- list(x=x, p=p, n=n, time=time2, tau=tau, status=status, mu_beta=mu_beta,
    sigma_beta=sigma_beta, a_lambda=a_lambda, b_lambda=b_lambda, a_gamma=a_gamma,
    b_gamma=b_gamma, approach=0, type=0)

  fit <- rstan::optimizing(stanmodels$weibull, data=data_mle_PH, hessian=TRUE)

  return(fit)
}


bayes_weibull_PH <- function(time,status, x, ...){

  p <- ncol(x)
  n <- nrow(x)
  tau <- max(time)
  time2 <- time/tau

  a_gamma <- 1
  b_gamma <- 1
  a_lambda <- 0.1
  b_lambda <- 0.1
  mu_beta <- 0
  sigma_beta <- 1

  data_bayes_PH <- list(x=x, p=p, n=n, time=time2, tau=tau, status=status, mu_beta=mu_beta,
    sigma_beta=sigma_beta, a_lambda=a_lambda, b_lambda=b_lambda, a_gamma=a_gamma,
    b_gamma=b_gamma, approach=1, type=0)

  fit <- rstan::sampling(stanmodels$weibull, data=data_bayes_PH)

  return(fit)

}

mle_weibull_PO <- function(time, status, x, ...){

  p <- ncol(x)
  n <- nrow(x)
  tau <- max(time)
  time2 <- time/tau

  a_gamma <- 1
  b_gamma <- 1
  a_lambda <- 0.1
  b_lambda <- 0.1
  mu_beta <- 0
  sigma_beta <- 1

  data_mle_PO <- list(x=x, p=p, n=n, time=time2, tau=tau, status=status,
    mu_beta=mu_beta, sigma_beta=sigma_beta, a_lambda=a_lambda, b_lambda=b_lambda,
    a_gamma=a_gamma, b_gamma=b_gamma, approach=0, type=1)

  fit <- rstan::optimizing(stanmodels$weibull, data=data_mle_PO, hessian=TRUE)

  return(fit)
}


bayes_weibull_PO <- function(time, status, x, ...){

  p <- ncol(x)
  n <- nrow(x)
  tau <- max(time)
  time2 <- time/tau

  a_gamma <- 1
  b_gamma <- 1
  a_lambda <- 0.1
  b_lambda <- 0.1
  mu_beta <- 0
  sigma_beta <- 1

  data_bayes_PO <- list(x=x, p=p, n=n, time=time2, tau=tau, status=status,
    mu_beta=mu_beta, sigma_beta=sigma_beta, a_lambda=a_lambda, b_lambda=b_lambda,
    a_gamma=a_gamma, b_gamma=b_gamma, approach=1, type=1)

  fit <- rstan::sampling(stanmodels$weibull, data=data_bayes_PO)

  return(fit)

}




weibull_phpo <- function(formula, data, approach=c("mle","bayes"), type=c("PH","PO"), a_gamma=1, b_gamma=1, a_lambda=0.1, b_lambda=0.1, mu_beta=0, sigma_beta=1, ...){

  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  resp <- stats::model.response(mf)
  time <- resp[,1]
  status <- resp[,2]
  x <- stats::model.matrix(formula, data = mf, rhs = 1)
  n <- length(time)

  tau <- max(time)
  time <- time/tau
  labels <- colnames(x)[-1]
  x <- matrix(x[,-1], ncol=length(labels))
  p <- ncol(x)

  if(type=="PH"){
    if(approach=="mle"){
        fit <- mle_weibull_PH(time=time, x=x, status=status, tau=tau, ...)
    }else{
        fit <- bayes_weibull_PH(time=time, x=x, status=status, tau=tau, ...)}

  }else{
      if(approach=="mle"){
        fit <- mle_weibull_PO(time=time, x=x, status=status, tau=tau, ...)

      }else{
        fit <- bayes_weibull_PO(time=time, x=x, status=status, tau=tau, ...)}
  }

  output <- list(fit=fit)

  output$n <- n
  output$p <- p

  output$call <- match.call()
  output$formula <- formula
  output$terms <- stats::terms.formula(formula)
  output$mf <- mf
  output$labels <- labels
  output$approach <- approach
  output$type <- type

  class(output) <- "weibull_phpo"
  return(output)
}



