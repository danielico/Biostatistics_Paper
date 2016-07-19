
#
# Adjustment of the model linear and proportional for the covariate age
#-----------------------------------------------
# 

Modele.v2 <- function( donnees = tempdonnees )
{
  #-------------------------------------------------*
  # specification of the model                       |
  #-------------------------------------------------*

  my.model = "~ intnum + I(intnum^2) + I( (intnum-1)^2 * (intnum > 1) ) + I( (intnum-5)^2 * (intnum > 5) ) + agec"

  #-------------------------------------------------*
  # Ajustement of model                       |
  #-------------------------------------------------*

  if(exists("splitdon")){rm(splitdon)}

  splitdon        = split.data.v2( donnees , relative = T , bands = c(seq(0,1,by=0.05),seq(1.1,10,by=0.1)) )
  splitdon$dcatt  = splitdon$tik * splitdon$MUA
  splitdon        = splitdon[,c("AGE.DIAG.A","ANNEE.DIAG.A","Fail","tik","intnum","dcatt","MUA","Entry","Exit")]
  splitdon$unsurt = 1/(splitdon$intnum + 1 )
  splitdon$logt   = log( splitdon$intnum + 1 )
  splitdon$agec   = trunc(splitdon$AGE.DIAG.A) - 70
  splitdon$adiagc = trunc(splitdon$ANNEE.DIAG.A) - 2000

  tempa <- try( glm( paste("Fail" , my.model) , offset = log(tik) , scale = 1 , family = POISS.RS.SPLIT.R.GLM( MyData.dcatt = splitdon$dcatt ) ,                data = splitdon , control = list( trace = T ) ) )

 if (class(tempa)[1]=="try-error") temp = excess.rate.fp.ml.v4.vagec(splitdon, modele = my.model)

 if (class(tempa)[1]!="try-error"){
 if (tempa$converged==F){temp = excess.rate.fp.ml.v4.vagec(splitdon, modele = my.model                        )
                         if (temp$fit==F) temp = excess.rate.fp.ml.v4.vagec(splitdon, modele = my.model, theta.init = tempa$coef)
                        }

 if (tempa$converged==T) {temp = excess.rate.fp.ml.v4.vagec(splitdon, modele = my.model,theta.init = tempa$coef)
                         if (temp$fit==F) temp = excess.rate.fp.ml.v4.vagec(splitdon, modele = my.model) 
                         }
                                  }
                                       
  return( list( fit = temp$fit , loglik = temp$loglik , theta = temp$theta , sigma = temp$sigma ) )
}

















