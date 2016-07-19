
# Poisson link function written by Aurelien Belot
# for relative survival according to the Dickman approach
# adapted for R
#-------------------------------------------------------------------------


POISS.RS.SPLIT.R.GLM <- function(link="log", MyData.dcatt){
  poissRS <- make.link("log")
  poissRS$linkfun <- function(mu, dcatt=MyData.dcatt) log(mu-dcatt)
  poissRS$linkinv <- function(eta, dcatt=MyData.dcatt) exp(eta)+dcatt
  poissRS$mu.eta <- function(eta) exp(eta)
  variance <- function(mu) mu
  validmu <- function(mu) all(mu > 0)
  dev.resids <- function(y, mu, wt) 2 * wt * (y * log(ifelse(y ==
        0, 1, y/mu)) - (y - mu))
  aic <- function(y, n, mu, wt, dev) -2 * sum(dpois(y, mu,
        log = TRUE) * wt)
  initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the Poisson family")
        n <- rep.int(1, nobs)
        mustart <- y + 0.1
    })
  structure(list(family="poissonRS", link="poissonRS.link", linkfun=poissRS$linkfun,
    linkinv=poissRS$linkinv, variance=variance, dev.resids=dev.resids, aic=aic,
    mu.eta=poissRS$mu.eta, initialize=initialize, validmu=validmu,
    valideta=poissRS$valideta),
  class="family")
}

