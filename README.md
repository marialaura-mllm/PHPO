# PHPO
An R package for ajusting proportional odds and proportional hazard models using the base function modeled by the *Weibull* distribution.

How to install
-------
```r
library(devtools)
install_github("marialaura-mllm/PHPO")
```

How to use
---------
Use the function **weibull_phpo(formula, data, approach=c("mle","bayes"), type=c("PH","PO"), a_gamma=1, b_gamma=1, a_lambda=0.1, 
b_lambda=0.1, mu_beta=0, sigma_beta=1, ...)** to model the survival data. With the *approach* parameter you can choose the estimation 
method, **mle** for maximum likelihood estimation and **bayes** for estimation via bayesian inference, and with the *type* parameter 
you can choose the model type, **PH** for proportional hazard models and **PO** for proportional odds models.


Example:

```r
#load package
libray(PHPO) 
#load data
library(KMsurv)
data("larynx")
#adjusting the model
mod <- weibull_phpo(Surv(time,delta)~age, approach = "mle", type = "PH", data = larynx)
#observing the model adjust result
summary(mod)
```

