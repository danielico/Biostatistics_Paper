
## This fonction is the second step to generate simulated survival data

cdatasimulation.PourDiffusion <- function ( donnees , cens.admin = 2013 , min.adiag , modele.simule , beta )
{ 

  tab5      <- donnees$tab5
  tab2      <- donnees$tab2

  n = nrow(tab2)

  # =====================#
  # T1 : Time to death due to cancer        |
  # =====================#
  
  Rescale.in <- function(gl,a,b)
  {
    gl$nodes    <- gl$nodes*(b-a)/2+(a+b)/2
    gl$weights  <- gl$weights*(b-a)/2
    return(gl)
  }
  
  survieFunc.in <- function( x , beta , modele , GL , a = 0 )   
  {
    gg    = Rescale.in(GL,a,x["times"])
    myMat = model.matrix( as.formula(modele) , data.frame( intnum = gg$nodes , unsurt = 1/(gg$nodes+1) , logt = log(gg$nodes+1) , agec = as.numeric(x["agec"] ) ) )
    return( exp( -sum( exp(myMat%*%beta) * gg$weights ) ) )
  }
  
  pred.surv.in <- function(  times , agec , beta, modele)
  {
  	GL=gauss.quad(n=100,kind="legendre")
  	temp <- data.frame( times = rep( times , length(agec) ) , unsurt = rep( 1/(times) , length(agec) ) , logt = rep( log(times) , length(agec) ) , agec = rep( rep( agec , each = length(times) ) ) )                     
  	temp$surv <- apply( temp , 1 , FUN = survieFunc.in , beta = beta , modele = modele , GL = GL )
  	return(temp$surv)
  }
  
  F_inverse_BO_v2 <- function( y , m.agec , m.modele , m.beta , lower = 0 , upper = (cens.admin - min.adiag) )
  {
    return( uniroot( f = function(x){ 1 - pred.surv.in( times = x , agec = m.agec , modele = m.modele , beta = m.beta ) - y } , lower = lower, upper = upper)[1] )
  }

  tab2$runif.T1 = runif(nrow(tab2))
  
  T1 <- c()
  for ( i in 1:nrow(tab2) ) 
  { 
    rm(x.T1)
    x.T1  = try( silent = T , F_inverse_BO_v2( y = tab2$runif.T1[i] , m.agec = tab2$agec[i] , m.modele = modele.simule , m.beta = beta )$root )
    T1[i] = ifelse( class(x.T1) == "try-error" , (cens.admin - min.adiag) , x.T1 )
  }

  # =====================#
  # T2 : Time to death due to other causes        |
  # =====================#
 
  tab5$T2bis <- rexp(n * (cens.admin - min.adiag), tab5$MUA)
  
  # IndicatriceT2bis : If T2bis >= 1 ==> Indicator == 1, if T2bis < 1 ==> Indicator == T2bis 
  tab5$IndicatriceT2bis <- 1 * (tab5$T2bis >= 1) + tab5$T2bis * (tab5$T2bis < 1)			
   

  # IndicesT2bis : can see if T2bis < 1 for each patient
  CalcT2 <- function(i, tab5)
  { 
    ifelse( min( which( tab5$Ident == i & tab5$IndicatriceT2bis < 1 ) ) != Inf , min( which( tab5$Ident == i & tab5$IndicatriceT2bis < 1 ) ) , -1 )
  }

  # IndicesT2 gives, for each patient, indicator where T2bis < 1
  IndicesT2 <- unlist( sapply( 1:n , CalcT2 , tab5 ) )

  # Time of death due to other causes T2 = Number of interval (year) + T2bis
  
  T2 <- rep( NA , n )
  
  a     <- which( IndicesT2 != -1)
  T2[a] <- ( tab5$NumInt[IndicesT2[a]] + tab5$T2bis[IndicesT2[a]] )
  
  b     <- which( is.na(T2) == T )
  T2[b] <- tab2$tpspotmax[b]
  

  # Follow-up time  
  Temps       <- data.frame( T1 = T1 , T2 = T2 )
  DureeDeVie  <- apply(Temps,1,min)
  
  datasimulation <- data.frame(tab2,Temps,DureeDeVie)

  # Potential time of follow-up : tpspot 
  datasimulation$tpspot     <-  cens.admin - datasimulation$annee.entier
  datasimulation$DureeDeVie <-  pmin(datasimulation$DureeDeVie, datasimulation$tpspot)

	# Status : 1 ==> death, 0 ==> censored 
  datasimulation[, "etat"] <- (0 + 1 * (datasimulation[, "DureeDeVie"] < datasimulation[, "tpspot"]))

	# Cause of death : 1 ==> Cancer, 2 ==> Other causes 
  datasimulation$cause <-( 0
                          + 1 * (datasimulation$DureeDeVie == datasimulation$T1 & datasimulation$etat == 1)
                          + 2 * (datasimulation$DureeDeVie == datasimulation$T2 & datasimulation$etat == 1) )

	# ----- datasimulation : Final Dataframe ----- #

  return(datasimulation)

}

