
## This fonction is a first step to generate simulated survival data
  
# Parameters 
#-------------------------------------------------------
#  n  sample size (default 500)
#  beta/modele.simule:   coef and formula of true model
#  site: site label (defines the age distribution) 
#  cens.admin; year of the administrative censoring ( 01/01, default 2013)
#  temps.survie : Estimates of theoretical survival at this time (default c(1,3,5,10))
#  min.adiag /max.adiag : year of diagnosis min et max (default 1990/2010) 

#  *** WARNING ***
#-------------------------------------------------------
#  The dataframe tauxatt must be loaded 

#--------------------------------------------------------------------------------------------------------
# Output
#--------------------------------------------------------------------------------------------------------
# tab2 :        data.frame with n patients with 11 covariates 
#               (ident, sexe, annee, age, age.entier, annee.entier, agec, tpspotmax, cens.admin  )

# tab5 :		data.frame with follow-up interval of 1y for each patient (until max time of follow-up)
#               and importation of mua corresponding to each interval (for simulation time to death from other causes)
#				nrow=nb patients x nb year of follow-up, 11 covariates:
#                 5 covariates  interval-dependantes  (age.int,annee.int, NumInt,NumIntSup , MUA )
#				+ 4 variables patients-dependants (Ident, sexe, annee, age, agec )
#               + 2 variable invariantes (tpspotmax, cens.admin )

# survie.pop :  Theoretical survival for t=temps.survie=time of follow-up 



# ============================================================================================================

cDataDesign.PourDiffusion <- function ( n = 500 , beta , site , cens.admin = 2013 , temps.survie = c(1,3,5,10) , min.adiag = 1990 , max.adiag = 2010 , noeud.knot )
 
{

  tpspotmax <- cens.admin - min.adiag 
 
  Ident <- c(1:n)

  if ( site == "colon" | site == "pancreas" | site == "oeso" )
  {  
    set.seed(100)
    age   <- c(runif(n*0.25, min = 30, max = 65), runif(n*0.35, min = 65, max = 75), runif(n*0.40, min = 75, max = 85))
  }

  if (site=="sein"|site=="thyroide") 
  {  
    set.seed(100)
    age   <- c( runif(n*0.13, min = 30, max = 45) , runif(n*0.24, min = 45, max = 55), runif(n*0.25, min = 55, max = 65) , runif(n*0.23, min = 65, max = 75) , runif(n*0.15, min = 75, max = 85) )
  }

  set.seed(100)
  
  sexe <- ifelse ( runif(n) < 0.5, 2, 1 )
  if (site=="sein") sexe=rep(2,n)

  set.seed(100)
  
  annee <- runif(n, min = min.adiag, max = max.adiag)

	# ----- Creation of tables allowing to simulate survival time -----

  NumInt    <- c(0:(cens.admin - min.adiag - 1))
  NumIntSup <- NumInt + 1
    
  tab1      <- data.frame(NumInt, NumIntSup)

  tab2              <- data.frame(Ident, sexe, annee, age)  
  tab2$age.entier   <- trunc(tab2$age)
  tab2$annee.entier <- trunc(tab2$annee)
  tab2$agec         <- tab2$age.entier - round(mean(age))
  tab2$tpspotmax    <- rep(tpspotmax, n)
  tab2$cens.admin   <- rep(cens.admin, n)   
  noeud             <- noeud.knot - round(mean(age))
 
  mod.nlin.ph = list.model.simu( knot.agec = noeud )$ms.nlin.ph
  modele.simule <- mod.nlin.ph
    
  tab3 <- merge(tab1,tab2[ ,c("Ident","sexe", "annee", "age", "agec")])
  tab3$age.int   <- trunc(tab3$age) + tab3$NumInt
  tab3$age.int  <- ifelse(tab3$age.int > 99, 99, tab3$age.int)  
  tab3$annee.int <- trunc(tab3$annee) + tab3$NumInt 

  tab4 <- merge(tab3, tauxatt, by.x = c("age.int","annee.int","sexe"), by.y = c("AGEX","ANNEE","SEXE"))

  tab5 <- tab4[order(tab4$Ident,tab4$age.int),]

  # Creation of survie.tay
 
  survie.tay      = merge( data.frame(intnum=temps.survie) , data.frame(age.entier=min(tab2$age.entier):max(tab2$age.entier) ) ) 
  survie.tay      = merge( survie.tay , data.frame(annee.entier=min(tab2$annee.entier):max(tab2$annee.entier) ) )
  survie.tay$agec = survie.tay$age.entier - round(mean(age))
   
  # ----- Cumulative mortality hazard  ----- #

  TxMort.in = function( intnum , agec , beta , modele )
  { 
    d = data.frame( intnum = intnum , unsurt = 1/(intnum+1) , logt = log(intnum+1) , agec = agec )
     return( as.vector( exp( model.matrix(as.formula(modele),d)%*%beta ) ) )
  }
 
  Rescale.in <- function(gl,a,b)
  {
    gl$nodes    <- gl$nodes*(b-a)/2+(a+b)/2
    gl$weights  <- gl$weights*(b-a)/2
    return(gl)
  }
	
  surv.in <- function( x , beta , modele , GL , a = 0 )   
  {    
    gg   = Rescale.in( GL , a , x["intnum"] )
    myMat= model.matrix( as.formula(modele) , data.frame(intnum = gg$nodes , unsurt = 1/(gg$nodes+1) , logt = log(gg$nodes+1) , agec = as.numeric(x["agec"]) ) )
    return( exp(-sum(exp(myMat%*%beta)*gg$weights ) ) )
	}     
 
  # Estimates of the excess mortality hazard and of the survival 
 
  survie.tay$Taux     <- TxMort.in( intnum = survie.tay$intnum , agec = survie.tay$agec , beta = beta , modele = modele.simule )
  survie.tay$Survie   <- apply( survie.tay[,c("intnum","agec")] , 1 , FUN = surv.in , beta = beta , modele = modele.simule , GL = gauss.quad( n = 100 , kind = "legendre" ) )
  
  tab.w       = data.frame( age.entier = min(tab2$age.entier):max(tab2$age.entier) )
  tab.w       = merge( tab.w , data.frame( annee.entier = min(tab2$annee.entier):max(tab2$annee.entier) ) )
  tab.w$Poids = rep( -1 , nrow(tab.w) )
  for ( k in 1:nrow(tab.w) )
  {
    tab.w$Poids[k] = sum(tab2$age.entier == tab.w$age.entier[k]&tab2$annee.entier == tab.w$annee.entier[k])/n 
  }

  survie.tay = merge( survie.tay , tab.w , by = c("age.entier","annee.entier") )

  # Theoretical survival 
   
  temp              = survie.tay[, c("intnum","Poids","Survie") ]
  temp$SP           = temp$Survie*temp$Poids
  survie.pop        = aggregate( temp[,c("SP")] , by = list( intnum = temp$intnum ) , FUN = sum )
  names(survie.pop) = c( "temps" , "Survie" )
  
  return( list( tab5 = tab5 , tab2 = tab2 , survie.pop = survie.pop , mod.nlin.ph = mod.nlin.ph ) )
}