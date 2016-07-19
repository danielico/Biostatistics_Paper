
# Function used to split the data
#-------------------------------------------------------------------------



split.data.v2 <- function(jeudata, bands, relative )
{

  jeudata$AGE.exit      <- jeudata$AGE.DIAG.A + jeudata$SUIVI.A
  jeudata$AGE.exit      <- ifelse(jeudata$AGE.exit >= 99, 99, jeudata$AGE.exit)
  jeudata$AGE.exit.r    <- floor(jeudata$AGE.exit)

  jeudata$ANNEE.exit.r  <- floor(jeudata$ANNEE.DIAG.A + jeudata$SUIVI.A)

  jeudata <- merge(jeudata,tauxatt,by.x=c("ANNEE.exit.r","AGE.exit.r","SEXE"),by.y=c("ANNEE","AGEX","SEXE"))

  splitjeudata          <- lexis( entry  = 0       ,
                                  exit   = SUIVI.A ,
                                  fail   = ETAT    ,
                                  breaks = bands   ,
                                 include = list(IDENPAT, SEXE, AGE.DIAG.A, AGEC.DIAG.A, AGE.exit.r, ANNEE.exit.r, ANNEEC.DIAG.A, ANNEE.DIAG.A, ETAT , MUA),
                                 data = jeudata)

  splitjeudata$tik      <-  splitjeudata$Exit - splitjeudata$Entry


  splitjeudata$intnum <- (splitjeudata$Entry + splitjeudata$Exit) / 2

  splitjeudata[splitjeudata$Fail == 1, c("intnum")] <- splitjeudata[splitjeudata$Fail == 1, c("Exit")]

  return(splitjeudata)

}
