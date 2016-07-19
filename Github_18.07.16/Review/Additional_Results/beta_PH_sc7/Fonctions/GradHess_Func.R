########################################################################
# Necessary fonctions to estimate hazard functions, cumulative hazard  #
# functions and conponents of the first and second derivatives of the  #
# loglikelihood                                                        #
########################################################################

Baseline.Haz <- function(t,coef,p){ 
    exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2)%*%coef)
}

d1.Baseline.t <- function(t,coef,p){
    t*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2)%*%coef)
}

d1.Baseline.t2 <- function(t,coef,p){
    t^2*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2)%*%coef)
}

d1.Baseline.t2p <- function(t,coef,p,wp){ 
    (t>p[wp])*(t-p[wp])^2*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2)%*%coef)
}

d2.Baseline.t.t2 <- function(t,coef,p){
    t^3*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2)%*%coef)
}

d2.Baseline.t.t2p <- function(t,coef,p,wp){
    t*(t>p[wp])*(t-p[wp])^2*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2)%*%coef)
}

d2.Baseline.t2.t2 <- function(t,coef,p){
    t^4*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2)%*%coef)
}

d2.Baseline.t2.t2p <- function(t,coef,p,wp){
    t^2*(t>p[wp])*(t-p[wp])^2*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2)%*%coef)
}

d2.Baseline.t2p.t2p <- function(t,coef,p,p1,p2){
    (t>p[p1])*(t-p[p1])^2*(t>p[p2])*(t-p[p2])^2*exp(cbind(1,t,t^2,(t>p[1])*(t-p[1])^2,(t>p[2])*(t-p[2])^2)%*%coef)
}

Time.2.Knot <- function(t,p){
    (t>p)*(t-p)^2
}
