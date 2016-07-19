# ------------------------------------------------------------------------------------ #
# -- PH test to check the proportional assumption of excess mortality hazard models -- #
# ------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------ #

# --- Initialisation --- #

# n.sim : Number of simulated Gaussian processes
# NbFichier : Number of simulated datasets on which we want to do the test
# res : dataframe containing the pvalue of the test for each simulated datasets

# --- Data --- #

# tempdonnees$IDENPAT         : Identification patient 
# tempdonnees$SEXE            : Sexe of the patient
# tempdonnees$SUIVI.J         : We add 1 day of follow-up if the date of diagnosis = date of death
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
# My.Score.X.t      : Score vector for the covariates of interest (observed process)
# My.pD.tau         : In the paper, W_hat = n^(-0.5)( D1(t) - I(theta,t).I(theta,tau).D1(tau) ) --> My.pD.tau * G = D1(tau) for all variables, with G following a std normal law
# My.List.pDX.t     : In the paper, W_hat = n^(-0.5)( D1(t) - I(theta,t).I(theta,tau).D1(tau) ) --> My.List.pDX.t * G = D1(t) for the covariate of interest at each time of follow-up for each individual, with G following a std normal law 
# My.List.Hess.X.t  : - seconde derivative of logV 
# My.Sigma          : Information matrix
# Temp.Gauss.PH     : Simulated Gaussian processes
# Mat.max.PH        : Maximum value of each simulated gaussian processes

# --- Results : p-values --- #

# p.value.PH        : pvalue

# ------------------------------------------------------------------------------------ #

library(miscTools) 

# Paths to change 

sink("S:/etude25/2679_Danieli_These/Review/Test_PH_nsim1000/lin_nph_lambda0.colon_NPH_Faible_1000/2000_patients/My_Test_PH_SigmaNum.lis")
         
chemin     <- ("S:/etude25/2679_Danieli_These/Review/")
chemin.bis <- ("S:/etude25/2679_Danieli_These/Review/Test_PH_nsim1000/")
chemin2    <- ("Simulations/lin_nph_lambda0.colon_NPH_Faible_1000/2000_patients/")
chemin3    <- ("Test_PH_nsim1000/lin_nph_lambda0.colon_NPH_Faible_1000/2000_patients/")

sourceDir <- function(path, trace = TRUE, ...)
{
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$"))
    {
        source(file.path(path, nm), ...)
    }
}

sourceDir( paste( chemin.bis ,"Fonctions", sep = "" ) )

# Load of the data

load( paste( chemin , chemin2 , "ListDataSimulation.lin.nph.RData" , sep="" ) )

tauxatt <- dget( file = paste( chemin , "Simulations/MUA/muaDF8933.dat",sep=""))

n.sim         <- 1000
NbFichier     <- 1000

rm(res)
res=data.frame(fichier=1:NbFichier, fit=rep(NA,NbFichier), p.value=rep(NA, NbFichier) )

for ( Fich in 1:NbFichier )
{
  rm(Modele.k2, temp, tempa)

  cat("Numero Fichier \n")
  print(Fich)

  # --- Data --- #
  
  tempdonnees = ListDataSimulation.lin.nph[[Fich]]
  
  tempdonnees$SUIVI.J         <- tempdonnees$DureeDeVie * 365.24
  tempdonnees$SUIVI.J         <- ifelse(tempdonnees$SUIVI.J == 0, 1, tempdonnees$SUIVI.J) 
  tempdonnees$SUIVI.A         <- tempdonnees$SUIVI.J/365.24

  tempdonnees$SUIVI.J         <- round(tempdonnees$SUIVI.J,9)
  tempdonnees$SUIVI.A         <- round(tempdonnees$SUIVI.A,9)

  # ----- > Age

  tempdonnees$AGE.DIAG.A      <- tempdonnees$age.entier
  tempdonnees$AGEC.DIAG.A     <- tempdonnees$agec

  # ----- > Year of diagnosis

  tempdonnees$ANNEE.DIAG.A    <- tempdonnees$annee.entier
  tempdonnees$ANNEEC.DIAG.A   <- tempdonnees$annee.entier - 1960

  # ----- > Identpat/Sex/Status

  tempdonnees$IDENPAT         <- tempdonnees$Ident
  tempdonnees$SEXE            <- tempdonnees$sexe
  tempdonnees$ETAT            <- tempdonnees$etat

  # ----- > 10 years of follow-up

  tempdonnees$ETAT            <- ifelse( tempdonnees$SUIVI.A > 10         , 0        , tempdonnees$ETAT)
  tempdonnees$SUIVI.J         <- ifelse( tempdonnees$SUIVI.J > 10*365.24  , 10*365.24 , tempdonnees$SUIVI.J)
  tempdonnees$SUIVI.A         <- ifelse( tempdonnees$SUIVI.A > 10         , 10        , tempdonnees$SUIVI.A)
 
  tempdonnees                 <- tempdonnees[order(tempdonnees$SUIVI.A,tempdonnees$IDENPAT),]
  
  # --- Time of follow-up --- #

  Tps.Unique    <- round(unique(tempdonnees$SUIVI.A),9) 
  Tps.Unique    <- Tps.Unique[order(Tps.Unique)]

  n.tps         <- length(Tps.Unique)
  n.pat         <- dim(tempdonnees)[1]
  
  # --- Model --- #

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
  
  MyOrder               <- tempdonnees$Ident
  MyTemp                <- tempdonnees
  MyTemp$AGE.exit       <- MyTemp$AGE.DIAG.A + MyTemp$SUIVI.A - 0.00001
  MyTemp$AGE.exit       <- ifelse( MyTemp$AGE.exit >= 99 , 99 , MyTemp$AGE.exit )
  MyTemp$AGE.exit.r     <- floor(MyTemp$AGE.exit)
  MyTemp$ANNEE.exit.r   <- floor( MyTemp$ANNEE.DIAG.A + MyTemp$SUIVI.A - 0.00001 )
  
  MyTemp2 <- merge( MyTemp , tauxatt , by.x = c( "ANNEE.exit.r" , "AGE.exit.r" , "SEXE" ) , by.y = c( "ANNEE" , "AGEX" , "SEXE" ) , sort = FALSE )
  MyData  <- MyTemp2[order(MyTemp2$Ident),][MyOrder,]
  rm(MyTemp,MyTemp2)

  # --- Construction of the process --- # 

  # cat("\nPreparation (1)\n")

  Covariables   <- model.matrix( ~ -1 + AGEC.DIAG.A , data = MyData ) 
  n.cov         <- dim(Covariables)[2]

  Knots         <- c(1,5)
  
  Event         <- MyData$ETAT
  
  Expected      <- MyData$MUA
  
  MyTime        <- MyData$SUIVI.A
  Time.cut      <- cut( MyTime , breaks = c( Tps.Unique[1] - 1 , Tps.Unique , Tps.Unique[n.tps]+1 ) )
  Idx.MyTime    <- as.numeric( Time.cut ) 
  
  coefficients  <- as.vector(Modele.k2$theta)
  
  Idx.Base      <- c(1:5)
  
  Idx.Cov       <- 5+1:n.cov
  
  leCoef        <- length(coefficients)
  
  coef.spline   <- round(coefficients[Idx.Base],8) 
  coef.cov      <- round(coefficients[Idx.Cov],8) 

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
  
  d2.Part1[,1]  <- d1.Part1[,1] 
  d2.Part1[,2]  <- d1.Part1[,2] 
  d2.Part1[,3]  <- d1.Part1[,3] 
  d2.Part1[,4]  <- d1.Part1[,3] 
  d2.Part1[,5]  <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t.t2,lower=0,upper=x,coef=coef.spline,p=Knots)$value)) 
  d2.Part1[,6]  <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t2.t2,lower=0,upper=x,coef=coef.spline,p=Knots)$value)) 
  d2.Part1[,7]  <- d1.Part1[,4] 
  d2.Part1[,8]  <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t.t2p,lower=0,upper=x,coef=coef.spline,p=Knots,wp=1)$value)) 
  d2.Part1[,9]  <- sapply(Tps.Unique,FUN=function(x)return(integrate(d2.Baseline.t2.t2p,lower=0,upper=x,coef=coef.spline,p=Knots,wp=1)$value)) 
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

  My.lambda.0 <- Baseline.Haz( t= MyTime , coef = coef.spline , p = Knots )
  exp.beta.z <- exp(Covariables%*%coef.cov)
  My.lambda <- exp.beta.z*My.lambda.0

  # cat("\nScore and Hessian \n")
  
  Frac.lambda.d1 <- Event*My.lambda.0/(My.lambda+Expected)
  Frac.lambda.d2 <- Event*My.lambda.0*Expected/(My.lambda+Expected)^2
  
  My.pD.tau       <- matrix( 0 , n.pat , leCoef )
  My.pD.tau[,1:5] <- as.vector( exp.beta.z * Frac.lambda.d1 ) * d1.Part2 
  
  for (j in 1:n.cov)
  { 
    My.pD.tau[,5+j] <- Covariables[,j]*exp.beta.z*Frac.lambda.d1
  }
  
  My.Score.X.t     <- matrix(0,n.tps,n.cov) 
  My.List.Hess.X.t <- list()
  My.List.pDX.t    <- list()
  
  for (j in 1:n.cov)
  {
    My.List.Hess.X.t[[j]] <- matrix(0,n.tps,leCoef)
    My.List.pDX.t[[j]]    <- matrix(0,n.pat,n.tps)
  }
  
  for (i in 1:n.tps)
  {
    Which <- sapply(Idx.MyTime,function(x) min(x,i))
    FL1   <- as.vector((Idx.MyTime<=i)*Frac.lambda.d1) 
    FL2   <- as.vector((Idx.MyTime<=i)*Frac.lambda.d2) 
    
    for (j in 1:n.cov)
    {
        Temp1 <- Covariables[,j]*exp.beta.z
        My.Score.X.t[i,j] <- t(Temp1)%*%(FL1*d1.Part2[,1] - d1.Part1[Which,1])
        My.List.Hess.X.t[[j]][i,1:5] <- t(Temp1)%*%(d1.Part1[Which,] - FL2*d1.Part2)
        for (k in 1:n.cov)
        {
            My.List.Hess.X.t[[j]][i,5+k] <- t(Covariables[,k]*Temp1)%*%(d1.Part1[Which,1] - FL2)
        }
        My.List.pDX.t[[j]][,i] <- Temp1*FL1
    }
  }
  
  # Hessian evaluated at Tau 
  My.Hessian <- rep(0,(leCoef*(leCoef+1)/2))
  FL2.tau    <- as.vector(Frac.lambda.d2)
  
  My.Hessian[1:15] <- t(exp.beta.z)%*% (d2.Part1[Idx.MyTime,] - FL2.tau*d2.Part2)
  Start <- 15
  for (j in 1:n.cov)
  {
      My.Hessian[Start+1:5] <- t(Covariables[,j]*exp.beta.z)%*%(d1.Part1[Idx.MyTime,] - FL2.tau*d1.Part2)
      
      for (k in 1:j)
      {
          My.Hessian[Start+5+k] <- t(Covariables[,j]*Covariables[,k]*exp.beta.z)%*%(d1.Part1[Idx.MyTime,1] - FL2.tau)
      }
      Start <- Start + 5 + j
  }
  
  My.Sigma <- solve(symMatrix(My.Hessian,byrow=TRUE))

  # --- Simulation processus --- # 
  
  My.Gaussian.Process.PH <- function( n.obs , pD.tau , pDX.t , RowHess.Xt.Sigma )
  {
      G     <- rnorm( n.obs , 0 , 1 )
      D.tau <- t(G%*%pD.tau)
      DX.t  <- G%*%pDX.t
      pW    <- DX.t - as.vector(RowHess.Xt.Sigma%*%D.tau)
      return(pW)
  }
  
  Mat.max.PH    <- matrix(0,n.sim,n.cov)
  p.value.PH    <- rep(0,n.cov)
  Gauss.Proc.PH <- list()
  My.Row.Hess.t <- matrix(0,leCoef,n.tps)
  Cst.Mult      <- 1/sqrt(n.pat)
  
  set.seed(Fich)
  
  for (j in 1:n.cov)
  {
      Temp.Gauss.PH <- matrix(0,n.tps,n.sim)
  
      My.pDX.t <- My.List.pDX.t[[j]]
      
      My.RowHess.Xt.Sigma <- My.List.Hess.X.t[[j]]%*%My.Sigma
  
      for (k in 1:n.sim)
      {
        Temp.Gauss.PH[,k] <- My.Gaussian.Process.PH( n.obs = n.pat , pD.tau = My.pD.tau , pDX.t = My.pDX.t , RowHess.Xt.Sigma = My.RowHess.Xt.Sigma )
      }
  
      Temp.Gauss.PH      <- Cst.Mult * Temp.Gauss.PH
      Mat.max.PH[,j]     <- apply(abs(Temp.Gauss.PH),2,max)
      p.value.PH[j]      <- sum( Mat.max.PH[,j] >= max(abs(My.Score.X.t[,j]))*Cst.Mult)/n.sim
      Gauss.Proc.PH[[j]] <- Temp.Gauss.PH
  }
  
  print(p.value.PH)

  res$p.value[res$fichier==Fich]=p.value.PH
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

plot( Tps.Unique                   ,
    (My.Score.X.t[,j]) * Cst.Mult ,
    col = "red"                   ,
    ylim = c(-20,20)              ,
    type = "l"                    ,
    main = "Test PH"  ,
    xlab = "Time since diagnosis",
    ylab = "Processes"                      )

for ( i in 1:1000 )
{
 lines( Tps.Unique , Temp.Gauss.PH[,i] , col = "grey" )
}

lines( Tps.Unique , (My.Score.X.t[,j]) * Cst.Mult , col = "red" , ylim = c(-20,20) )

dev.off()


