# ---------------------------------------------------------------------------- #
# -- FF test to check the functional form of excess mortality hazard models -- #
# ---------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------ #

# --- Initialisation --- #

# n.sim : Number of simulated Gaussian processes
# NbFichier : Number of simulated datasets on which we want to do the test
# res : dataframe containing the pvalue of the test for each simulated datasets

# --- Data --- #

# tempdonnees$IDENPAT         : Identification patient 
# tempdonnees$SEXE            : Sexe of the patient
# tempdonnees$SUIVI.A         : Follow-up in years
# tempdonnees$AGE.DIAG.A      : Age of the patient in year (rounded)
# tempdonnees$AGEC.DIAG.A     : Centered age of the patient
# tempdonnees$ANNEE.DIAG.A    : year of diagnosis (rounded)
# tempdonnees$ANNEEC.DIAG.A   : centered year of diagnosis
# tempdonnees$ETAT            : Status of the patient ( 0 = censored / 1 = dead )

# tempdonnees : Data sorted by time of follow-up

# tauxatt : life table (dataframe)

# --- Building of the process --- #

# Tps.Unique        : (ordered) different times of follow-up   
# n.tps             : Number of different times of follow-up in the population
# n.pat             : Sample size
# Modele.k2         : objects obtain from the model adjustment
# MyData            : Data with expected hazard mortality
# Covariables       : covariates matrix 
# n.cov             : Number of covariates
# Knots             : Position des noeuds interieurs du spline pr le taux de base
# Event             : Event indicator 
# Expected          : Expected hazard at censoring time
# MyTime            : ordered Follow-up time
# Time.cut          : cut points
# Idx.MyTime        : Indicate position of the follow-up time of each individual in the vector Tps.Unique
# coefficients      : Vector of the covariate effects
# Idx.Base          : Position des coefficients associes au taux de base
# Idx.Cov           : Position of the coefficients associated to covariables
# leCoef            : Number of coefficients in the model
# coef.spline       : Coefficient associated to the baseline hazard
# coef.cov          : Coefficient associated to the covariates 
# My.lambda.0       : Baseline hazard for each individual
# exp.beta.z        : Vector of exp(beta*x)
# My.lambda         : Individual excess hazard
# Mart.Obs          : Martingale residuals
# My.pD.tau         : In the paper, W_hat = n^(-0.5)( D1(t) - I(theta,t).I(theta,tau).D1(tau) ) --> My.pD.tau * G = D1(tau) for all variables, with G following a std normal law
# My.Sigma          : Information matrix
# Unique.X          : Distinct values of each covariates
# Idx.X             : Position of the individual covariate value in Unique.X
# My.MCum.x         : Observed process for each covariates
# My.List.J.x       : Matrix J components estimates
# My.List.M.x       : Estimates of P(x) for each cov
# Temp.Gauss.FF     : Simulated Gaussian processes
# Mat.max.FF        : Maximum value of each simulated gaussian processes

# --- Results : p-values --- #

# p.value.FF        : pvalue

# ------------------------------------------------------------------------------------ #

sink("S:/etude25/2679_Danieli_These/Review/test_FF_nsim1000/nlin_ph_lambda0.colon_NLIN_Moyen_1000/1000_patients/My_Test_FF_SigmaNum.lis")

library(miscTools) 

# Paths to change 

chemin     <- ("S:/etude25/2679_Danieli_These/Review/")
chemin.bis <- ("S:/etude25/2679_Danieli_These/Review/test_FF_nsim1000/")
chemin2    <- ("Simulations/nlin_ph_lambda0.colon_NLIN_Moyen_1000/1000_patients/")
chemin3    <- ("test_FF_nsim1000/nlin_ph_lambda0.colon_NLIN_Moyen_1000/1000_patients/")

sourceDir <- function(path, trace = TRUE, ...)
{
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$"))
  {
    source(file.path(path, nm), ...)
  }
}

sourceDir( paste( chemin.bis ,"Fonctions", sep = "" ) )

# Load of the data

load( paste( chemin , chemin2 , "ListDataSimulation.nlin.ph.RData" , sep="" ) )

tauxatt <- dget( file = paste( chemin , "Simulations/MUA/muaDF8933.dat",sep=""))

n.sim         <- 1000
NbFichier     <- 1000

rm(res)
res=data.frame(fichier=1:NbFichier, fit=rep(NA,NbFichier), p.value=rep(NA, NbFichier) )
 
for ( Fich in 1:NbFichier )
{
  
  cat("Numero Fichier \n")
  print(Fich)
  
  # --- Data --- #
  
  tempdonnees = ListDataSimulation.nlin.ph[[Fich]]
  
  # ----- > Durée De Suivi
  
  tempdonnees$SUIVI.J         <- tempdonnees$DureeDeVie * 365.24
  tempdonnees$SUIVI.J         <- ifelse(tempdonnees$SUIVI.J == 0, 1, tempdonnees$SUIVI.J) # ajout de 1 jour si la date de diagnostic = date de décès
  tempdonnees$SUIVI.A         <- tempdonnees$SUIVI.J/365.24
  
  tempdonnees$SUIVI.J         <- round(tempdonnees$SUIVI.J,9)
  tempdonnees$SUIVI.A         <- round(tempdonnees$SUIVI.A,9)
  
  # ----- > Age
  
  tempdonnees$AGE.DIAG.A      <- tempdonnees$age.entier
  tempdonnees$AGEC.DIAG.A     <- tempdonnees$agec
  
  # ----- > Year of diagnosis
  
  tempdonnees$ANNEE.DIAG.A    <- tempdonnees$annee.entier
  tempdonnees$ANNEEC.DIAG.A   <- tempdonnees$annee.entier - 1960
  
  # ----- >
  
  tempdonnees$IDENPAT         <- tempdonnees$Ident
  tempdonnees$SEXE            <- tempdonnees$sexe
  tempdonnees$ETAT            <- tempdonnees$etat
  
  tempdonnees$ETAT            <- ifelse( tempdonnees$SUIVI.A > 10         , 0        , tempdonnees$ETAT)
  tempdonnees$SUIVI.J         <- ifelse( tempdonnees$SUIVI.J > 10*365.24  , 10*365.24 , tempdonnees$SUIVI.J)
  tempdonnees$SUIVI.A         <- ifelse( tempdonnees$SUIVI.A > 10         , 10        , tempdonnees$SUIVI.A)
  
  tempdonnees                 <- tempdonnees[order(tempdonnees$SUIVI.A,tempdonnees$IDENPAT),]
  
  MyData <- tempdonnees 
  
  # --- Time of follow-up --- #
  
  Tps.Unique   <- round(unique(MyData$SUIVI.A),9) # Pourquoi arrondir ici ? !!!
  Tps.Unique   <- Tps.Unique[order(Tps.Unique)]
  
  n.tps <- length(Tps.Unique)
  n.pat <- dim(MyData)[1]
  
  # --- Model --- #
  
  rm(Modele.k2, temp, tempa)
  cat("\nModele ML\n")
  Modele.k2 <- Modele.v2( donnees = tempdonnees )
  
  if( Modele.k2$fit == F )
  { 
    print("Probleme : le modele n a pas converge...")
    res$fit[res$fichier==Fich]=F
    next
  }

  if( Modele.k2$fit == T )
  { 
    print("le modele est OK ")
    res$fit[res$fichier==Fich]=T
  }
  
  # --- Expected hazard mortality --- #
  
  tempTxAtt               <- matrix(0,n.pat,11)
  MyDataSlim              <- MyData[,c("IDENPAT","SUIVI.A","annee.entier","sexe","age.entier")]
  MyDataSlim$age.entier   <- MyDataSlim$age.entier - 1
  MyDataSlim$annee.entier <- MyDataSlim$annee.entier - 1
  
  for (an in 1:10)
  {
      MyDataSlim$age.entier                              <- MyDataSlim$age.entier + 1
      MyDataSlim$age.entier[MyDataSlim$age.entier >= 99] <- 99
      MyDataSlim$annee.entier                            <- MyDataSlim$annee.entier + 1
      temp                                               <- merge(MyDataSlim, tauxatt, by.x=c("annee.entier","sexe","age.entier"), by.y=c("ANNEE","SEXE","AGEX"))
      tempTxAtt[,an+1]                                   <- temp$MUA[order(temp$SUIVI.A,temp$IDENPAT)]
  }
  
  Floor.Time                  <- floor(MyData$SUIVI.A)
  Remain.Time                 <- MyData$SUIVI.A - Floor.Time
  Remain.Time[Remain.Time==0] <- 1
  Floor.Time                  <- Floor.Time - 1*(Remain.Time==1)
  CumTxAtt                    <- t(apply(tempTxAtt,1,cumsum))
  
  Expected                    <- sapply(1:n.pat,function(i){tempTxAtt[i,Floor.Time[i]+2]})
  Temp.Lambda.Expect          <- sapply(1:n.pat,function(i){CumTxAtt[i,Floor.Time[i]+1]})
  My.Lambda.Expect            <- Temp.Lambda.Expect + Remain.Time*Expected
  
  # --- Construction of the process --- # 

  #cat("\nPreparation (1)\n")
  
  names.cov <- c("AGEC.DIAG.A")
  Covariables <- model.matrix(~-1+AGEC.DIAG.A,data=MyData)
  n.cov <- dim(Covariables)[2]
  
  Knots <- c(1,5)
  
  Event <- MyData$ETAT
  
  MyTime <- MyData$SUIVI.A
  Time.cut <- cut(MyTime,breaks=c(Tps.Unique[1]-1,Tps.Unique,Tps.Unique[n.tps]+1))
  Idx.MyTime <- as.numeric(Time.cut) 
  
  coefficients <- as.vector(Modele.k2$theta)
  
  Idx.Base <- c(1:5)
  
  Idx.Cov <- 5+1:n.cov
  
  leCoef <- length(coefficients)
  
  coef.spline <- round(coefficients[Idx.Base],8) # Pourquoi arrondir ici ? !!!
  coef.cov    <- round(coefficients[Idx.Cov],8) # Pourquoi arrondir ici ? !!!

  #cat("\nPreparation (2)\n")

  # For the score, we need the following quantities
  d1.Part1 <- matrix(0,n.tps,length(coef.spline))
  d1.Part2 <- matrix(0,n.pat,length(coef.spline))
  
  d1.Part1[,1] <- sapply(Tps.Unique,FUN=function(x)return(integrate(Baseline.Haz,lower=0,upper=x,coef=coef.spline,p=Knots)$value)) 
  d1.Part1[,2] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d1.Baseline.t,lower=0,upper=x,coef=coef.spline,p=Knots)$value))
  d1.Part1[,3] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d1.Baseline.t2,lower=0,upper=x,coef=coef.spline,p=Knots)$value)) 
  d1.Part1[,4] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d1.Baseline.t2p,lower=0,upper=x,coef=coef.spline,p=Knots,wp=1)$value))
  d1.Part1[,5] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d1.Baseline.t2p,lower=0,upper=x,coef=coef.spline,p=Knots,wp=2)$value)) 
  
  d1.Part2[,1] <- 1
  d1.Part2[,2] <- MyTime
  d1.Part2[,3] <- MyTime^2
  d1.Part2[,4] <- Time.2.Knot(MyTime,Knots[1])
  d1.Part2[,5] <- Time.2.Knot(MyTime,Knots[2])
  
  # For the Hessian matrix, we need the following quantities
  d2.Part1 <- matrix(0,n.tps,length(coef.spline)*(length(coef.spline)+1)/2)
  d2.Part2 <- matrix(0,n.pat,length(coef.spline)*(length(coef.spline)+1)/2)
  
  d2.Part1[,1] <- d1.Part1[,1] 
  d2.Part1[,2] <- d1.Part1[,2] 
  d2.Part1[,3] <- d1.Part1[,3]
  d2.Part1[,4] <- d1.Part1[,3]
  d2.Part1[,5] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t.t2,lower=0,upper=x,coef=coef.spline,p=Knots)$value))
  d2.Part1[,6] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t2.t2,lower=0,upper=x,coef=coef.spline,p=Knots)$value))
  d2.Part1[,7] <- d1.Part1[,4]
  d2.Part1[,8] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t.t2p,lower=0,upper=x,coef=coef.spline,p=Knots,wp=1)$value))
  d2.Part1[,9] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t2.t2p,lower=0,upper=x,coef=coef.spline,p=Knots,wp=1)$value)) 
  d2.Part1[,10] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t2p.t2p,lower=0,upper=x,coef=coef.spline,p=Knots,p1=1,p2=1)$value)) 
  d2.Part1[,11] <- d1.Part1[,5] 
  d2.Part1[,12] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t.t2p,lower=0,upper=x,coef=coef.spline,p=Knots,wp=2)$value)) 
  d2.Part1[,13] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t2.t2p,lower=0,upper=x,coef=coef.spline,p=Knots,wp=2)$value)) 
  d2.Part1[,14] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t2p.t2p,lower=0,upper=x,coef=coef.spline,p=Knots,p1=1,p2=2)$value)) 
  d2.Part1[,15] <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t2p.t2p,lower=0,upper=x,coef=coef.spline,p=Knots,p1=2,p2=2)$value)) 
  
  d2.Part2[,1] <- 1 
  d2.Part2[,2] <- d1.Part2[,2] 
  d2.Part2[,3] <- d1.Part2[,3] 
  d2.Part2[,4] <- d1.Part2[,3] 
  d2.Part2[,5] <- d1.Part2[,3]*d1.Part2[,2] 
  d2.Part2[,6] <- d1.Part2[,3]*d1.Part2[,3] 
  d2.Part2[,7] <- d1.Part2[,4] 
  d2.Part2[,8] <- d1.Part2[,4]*d1.Part2[,2] 
  d2.Part2[,9] <- d1.Part2[,4]*d1.Part2[,3] 
  d2.Part2[,10] <- d1.Part2[,4]*d1.Part2[,4] 
  d2.Part2[,11] <- d1.Part2[,5] 
  d2.Part2[,12] <- d1.Part2[,5]*d1.Part2[,2] 
  d2.Part2[,13] <- d1.Part2[,5]*d1.Part2[,3] 
  d2.Part2[,14] <- d1.Part2[,5]*d1.Part2[,4] 
  d2.Part2[,15] <- d1.Part2[,5]*d1.Part2[,5] 
  
  My.lambda.0 <- Baseline.Haz(t=MyTime,coef=coef.spline,p=Knots)
  exp.beta.z <- exp(Covariables%*%coef.cov)
  My.lambda <- exp.beta.z*My.lambda.0

 # cat("\nMartingales and Hessian \n")

  Frac.lambda.d1 <- Event*My.lambda.0/(My.lambda+Expected)
  Frac.lambda.d2 <- Event*My.lambda.0*Expected/(My.lambda+Expected)^2
  
  Mart.Obs <- Event - exp.beta.z*d1.Part1[Idx.MyTime,1] - My.Lambda.Expect
  
  My.pD.tau       <- matrix(0,n.pat,leCoef)  
  My.pD.tau[,1:5] <- as.vector(exp.beta.z*Frac.lambda.d1)*d1.Part2 # Pour les parametres du taux de base
  
  for (j in 1:n.cov)
  { 
      My.pD.tau[,5+j] <- Covariables[,j]*exp.beta.z*Frac.lambda.d1
  }
  
  My.Hessian       <- rep(0,leCoef*(leCoef+1)/2)
  Temp1            <- Covariables[,j]*exp.beta.z
  FL2              <- as.vector(Frac.lambda.d2)
  My.Hessian[1:15] <- t(exp.beta.z)%*% (d2.Part1[Idx.MyTime,] - FL2*d2.Part2)
  Start <- 15
  
  for (j in 1:n.cov)
  {
      My.Hessian[Start+1:5] <- t(Covariables[,j]*exp.beta.z)%*%(d1.Part1[Idx.MyTime,] - FL2*d1.Part2)
      
      for (k in 1:j)
      {
          My.Hessian[Start+5+k] <- t(Covariables[,j]*Covariables[,k]*exp.beta.z)%*%(d1.Part1[Idx.MyTime,1] - FL2)
      }
      
      Start <- Start + 5 + j
  }
  
  My.Sigma <- solve(symMatrix(My.Hessian,byrow=TRUE))
  
  Unique.X    <- list() 
  Idx.X       <- list() 
  My.MCum.x   <- list() 
  My.List.J.x <- list() 
  My.List.M.x <- list() 
  lUX         <- rep(0,n.cov)

  # J
  
  J.0 <- matrix(0,n.pat,leCoef)
  J.0[,1:5] <- as.vector(exp.beta.z)*d1.Part1[Idx.MyTime,]
  
  for (j in 1:n.cov)
  {
      J.0[,5+j] <- Covariables[,j]*exp.beta.z*d1.Part1[Idx.MyTime,1]
  }
  
  for (j in 1:n.cov)
  {
  
      Temp.Unique.X <- sort(unique(Covariables[,j]))
      lUX[j] <- leU <- length(Temp.Unique.X)
  
      cut.X      <- cut( Covariables[,j] , breaks = c( Temp.Unique.X[1]-1 , Temp.Unique.X , Temp.Unique.X[leU]+1 ) )
      Temp.Idx.X <- as.numeric(cut.X)
  
      Temp.MCum.x <- rep(0,leU)
      Temp.M.x    <- matrix(0,n.pat,leU)
      Temp.J.x    <- matrix(0,leU,leCoef)
  
      for (k in 1:leU)
      {
          Indic          <- (Temp.Idx.X<=k)
          Temp.MCum.x[k] <- Indic%*%Mart.Obs
          Temp.M.x[,k]   <- Indic*Event
          Temp.J.x[k,]   <- Indic%*%J.0
      }
  
      My.MCum.x[[j]]   <- Temp.MCum.x
      My.List.M.x[[j]] <- Temp.M.x
      My.List.J.x[[j]] <- Temp.J.x
      Unique.X[[j]]    <- Temp.Unique.X
      Idx.X[[j]]       <- Temp.Idx.X
  }
 
  # --- Simulation processus --- # 
  
  cat("\nSimulation \n")
  
  My.Gaussian.Process.FF <- function( n.obs , pD.tau , M.x , RowJx.Sigma )
  {
      G     <- rnorm(n.obs,0,1)
      D.tau <- t(G%*%pD.tau)
      P.x   <- G%*%M.x
      pW    <- P.x - as.vector(RowJx.Sigma%*%D.tau)
      return(pW)
  }

  Mat.max.FF    <- matrix(0,n.sim,n.cov)
  p.value.FF    <- rep(0,n.cov)
  Gauss.Proc.FF <- list()
  Cst.Mult      <- 1/sqrt(n.pat)
  set.seed(Fich)

  for (j in 1:n.cov)
  {
  
      Temp.Gauss.FF <- matrix(0,lUX[j],n.sim)
  
      My.M.x         <- My.List.M.x[[j]]
      My.RowJx.Sigma <- My.List.J.x[[j]]%*%My.Sigma
  
      for (k in 1:n.sim)
      {
          Temp.Gauss.FF[,k] <- My.Gaussian.Process.FF( n.obs = n.pat , pD.tau = My.pD.tau , M.x = My.M.x , RowJx.Sigma = My.RowJx.Sigma )
      }
  
      Temp.Gauss.FF      <- Cst.Mult*Temp.Gauss.FF
      Mat.max.FF[,j]     <- apply(abs(Temp.Gauss.FF),2,max)
      p.value.FF[j]      <- sum( Mat.max.FF[,j] >= max(abs(My.MCum.x[[j]]))*Cst.Mult)/n.sim
      Gauss.Proc.FF[[j]] <- Temp.Gauss.FF
  }

  print(p.value.FF)
  
   res$p.value[res$fichier==Fich]=p.value.FF


}

save( res , file = paste( chemin , chemin3 , "res.RData" , sep = "" ) )
cat("\n resultas finaux res ........!! \n")
print(dim(res))
print(summary(res))
print(table(res$fit))
print(res[res$fit==F,])
cat("\n p final .......!!! \n")
print( table(res[res$fit==T, c("p.value")]<0.05)/nrow(res[res$fit==T,]) )

sink()

# -----> Graphes

bmp(paste( chemin , chemin3 , "Graphe.bmp" , sep = "" ))

 plot( Temp.Unique.X                   , 
      (My.MCum.x[[j]]) * Cst.Mult , 
      col = "red"                   , 
      ylim = c(-1,1)              , 
      type = "l"                    , 
      main = "Test FF"  ,
      xlab = "Centered age",
      ylab = "Processes"                      )
      
 for ( i in 1:1000 )
 {
   lines( Temp.Unique.X , Gauss.Proc.FF[[j]][,i] , col = "grey" )
 }

 lines( Temp.Unique.X , (My.MCum.x[[j]]) * Cst.Mult , col = "red" )

 dev.off()



