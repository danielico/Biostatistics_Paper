#-------------------------------------------------------*
#|                                                      |
#|  Vraisemblance, gradient, Hessien, algo Newton-raphson |
#|  General fonction allowing to adjust several types of models with the following variables : intnum, unsurt, logt, agec, adiagc
#|                                                      |
#|                                                      |
#-------------------------------------------------------*

excess.rate.fp.ml.v4.vagec = function( splitdonnees , modele , theta.init )
{

  i    = 1
  ll   = 100
  llold= 1

  p = ncol( model.matrix( as.formula(modele) , data.frame( intnum = 1 , agec = 1 , adiagc = 1 , unsurt = 1 , logt = 1 ) ) ) #nb de cov ds le modele

  if ( missing(theta.init) )
  {
    theta.init = rep(0,p)
  }

  thetaold = theta.init
  theta    = thetaold

	t0     <- splitdonnees$Entry
	t1     <- splitdonnees$Exit
	tm     <- (t0+t1)/2
	delta  <- splitdonnees$Fail
	mua    <- splitdonnees$MUA
	agec   <- splitdonnees$agec
	adiagc <- splitdonnees$adiagc
	unsurt <- splitdonnees$unsurt
	logt   <- splitdonnees$logt

	X.tm = model.matrix( as.formula(modele) , data.frame( intnum = tm , agec = agec , adiagc = adiagc , unsurt = unsurt , logt = logt ) )
	X.t0 = model.matrix( as.formula(modele) , data.frame( intnum = t0 , agec = agec , adiagc = adiagc , unsurt = unsurt , logt = logt ) )
	X.t1 = model.matrix( as.formula(modele) , data.frame( intnum = t1 , agec = agec , adiagc = adiagc , unsurt = unsurt , logt = logt ) )

while( abs(ll-llold) > 0.0001 | any( abs( ( theta-thetaold )/thetaold ) > 0.0001 ) )
{
  if (i > 50)
  {
    cat("message LRE: Ran out of iterations", i, "and did not converge ","\n","\n")
    break
  }

  if( i >= 2 )
  {
    llold     = ll
    thetaold  = theta
  }

   ftmold = as.numeric( exp( X.tm %*% thetaold ) )
   ft0old = as.numeric( exp( X.t0 %*% thetaold ) )
   ft1old = as.numeric( exp( X.t1 %*% thetaold ) )

	 h <-   t(X.t0) %*% ( ft0old * X.t0 * ( -(t1-t0)/6 ) ) + 4 *t(X.tm) %*% ( X.tm * ftmold * ( -(t1-t0)/6 ) ) + t(X.t1) %*% ( X.t1 * ft1old * ( -(t1-t0)/6 ) + delta * X.t1 * ft1old * mua / ( ft1old + mua )^2 )

	grad <- colSums( -(t1-t0)/6 * ( X.t0 * ft0old + 4 * X.tm * ftmold + X.t1 * ft1old ) + delta * ( X.t1 * ft1old )/( ft1old + mua ) )

	pas <- try( chol2inv(chol(-h)) %*% grad )

  if ( class(pas) == "Error" | class(pas) == "try-error" )
  {
    cat("iterations: ", i, "Class(pas)==Error OU try-error.......: ","\n","\n")
    break
  }

  theta = thetaold + pas


	ftm = as.numeric( exp( X.tm%*%theta ) )
	ft0 = as.numeric( exp( X.t0%*%theta ) )
	ft1 = as.numeric( exp( X.t1%*%theta ) )

	ll = sum( -(t1-t0)/6 * ( ft0 + 4*ftm + ft1 ) + delta * log( ft1 + mua ) )

   cat("iter : ",i,"\n",
      "thetaold= ", round(thetaold,4),"\n",
      "theta= ", round(theta,4),"\n",
      "abs((theta-thetaold)/thetaold)= ", round(abs((theta-thetaold)/thetaold),5),"\n",
      "ll-llold= ", round(ll-llold,5),"\n",
      "loglik.excess.rate.fp.ml.v4.vagec(theta)= ",round(ll,4),"\n","\n" )

   i = i + 1
}

if ( class(pas) == "Error" | class(pas) == "try-error" | i > 50 )
{

  fit=F
  sigma = matrix(-1,2,3)

}else{

  fit=T

  sigma = try( chol2inv( chol( -( t(X.t0) %*% ( ft0 * X.t0 * (-(t1-t0)/6)) + 4 * t(X.tm) %*% ( X.tm * ftm * (-(t1-t0)/6)) + t(X.t1) %*% ( X.t1 * ft1 * (-(t1-t0)/6) + delta * X.t1 * ft1 * mua / ( ft1 + mua )^2 ) ) ) ) )

}

if ( fit == T )
{
  theta = as.vector(theta)
  names(theta) = dimnames( model.matrix( as.formula(modele) , data.frame( intnum = 1 , agec = 1 , adiagc = 1 , unsurt = 1 , logt = 1 ) ) )[[2]]

  if ( class(sigma) != "try-error" )
  {
    dimnames(sigma)=list(names(theta),names(theta))
  }
}

return( list( fit = fit , loglik = ll , theta = theta , sigma = sigma ) )

}

