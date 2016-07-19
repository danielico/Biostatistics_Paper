########################################################################
# Necessary fonctions to estimate hazard functions, cumulative hazard  #
# functions and conponents of the first and second derivatives of the  #
# loglikelihood                                                        #
########################################################################

Baseline.Haz.A <- function(t,coef,p,cov){ 
    exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2,cov,cov*t,cov*t^2,cov*(t>p[1])*(t-p[1])^2,cov*(t>p[2])*(t-p[2])^2)%*%coef)
}

d1.Baseline.t.A <- function(t,coef,p,cov){
    t*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2,cov,cov*t,cov*t^2,cov*(t>p[1])*(t-p[1])^2,cov*(t>p[2])*(t-p[2])^2)%*%coef)
}

d1.Baseline.t2.A <- function(t,coef,p,cov){
    t^2*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2,cov,cov*t,cov*t^2,cov*(t>p[1])*(t-p[1])^2,cov*(t>p[2])*(t-p[2])^2)%*%coef)
}

d1.Baseline.t2p.A <- function(t,coef,p,wp,cov){ 
    (t>p[wp])*(t-p[wp])^2*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2,cov,cov*t,cov*t^2,cov*(t>p[1])*(t-p[1])^2,cov*(t>p[2])*(t-p[2])^2)%*%coef)
}

d2.Baseline.t.t2.A <- function(t,coef,p,cov){
    t^3*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2,cov,cov*t,cov*t^2,cov*(t>p[1])*(t-p[1])^2,cov*(t>p[2])*(t-p[2])^2)%*%coef)
}

d2.Baseline.t.t2p.A <- function(t,coef,p,wp,cov){
    t*(t>p[wp])*(t-p[wp])^2*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2,cov,cov*t,cov*t^2,cov*(t>p[1])*(t-p[1])^2,cov*(t>p[2])*(t-p[2])^2)%*%coef)
}

d2.Baseline.t2.t2.A <- function(t,coef,p,cov){
    t^4*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2,cov,cov*t,cov*t^2,cov*(t>p[1])*(t-p[1])^2,cov*(t>p[2])*(t-p[2])^2)%*%coef)
}

d2.Baseline.t2.t2p.A <- function(t,coef,p,wp,cov){
    t^2*(t>p[wp])*(t-p[wp])^2*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2,cov,cov*t,cov*t^2,cov*(t>p[1])*(t-p[1])^2,cov*(t>p[2])*(t-p[2])^2)%*%coef)
}

d2.Baseline.t2p.t2p.A <- function(t,coef,p,p1,p2,cov){
    (t>p[p1])*(t-p[p1])^2*(t>p[p2])*(t-p[p2])^2*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2,cov,cov*t,cov*t^2,cov*(t>p[1])*(t-p[1])^2,cov*(t>p[2])*(t-p[2])^2)%*%coef)
}

Time.2.Knot <- function(t,p){
    (t>p)*(t-p)^2
}
