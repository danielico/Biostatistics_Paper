#-------------------------------------------------------------------------------------------------------------
# Generation of the data 
#-------------------------------------------------------------------------------------------------------------

# Paths to change !!

chemin  <- c("S:/etude25/2679_Danieli_These/Review/")
chemin2 <- c("Simulations/nlin_ph_lambda0.colon_NLIN_Moyen_1000/500_patients/")

sink( paste( chemin , chemin2 , "Generation_Donnees_nlin_ph.lis" , sep = "") )

library(statmod)

set.seed(1352)

print(system("hostname", intern=T))
heure.debut=Sys.time()
print(heure.debut)

source( paste( chemin , "Simulations/Fonctions/cDataDesign.PourDiffusion.NLIN.R"     , sep = "") )
source( paste( chemin , "Simulations/Fonctions/cdatasimulation.PourDiffusion.R" , sep = "") )
source( paste( chemin , "Simulations/Fonctions/ListModeles.R"                   , sep = "") )

###### ----- > Expected mortality hazard

tauxatt = dget( file = paste( chemin , "Simulations/MUA/muaDF8933.dat" , sep = "") )

###### ----- > Parameter for the simulation 

load( file = paste( chemin , "Simulations/Parametres_Theoriques/colon_TxBase/parametre_theorique/coef.colon.lin.nph.RData"  , sep = "" ) )
load( file = paste( chemin , "Simulations/Parametres_Theoriques/nlin_ph_col/parametre_theorique/coef.col.nlin.ph.Rdata" , sep = "" ) )

coef.col.nlin.ph[1:5] <- coef.colon.lin.nph[1:5]

NbFichier = 1000
n         = 500

rm( mod.nlin.ph )

###### ----- > Data simulation

# DataDesign.lin.nph : Contains all the constant information 
# ListDataSimulation : Generation of time to death du to cancer and time to death due to other causes

cat("Heure debut \n")
print(date())

cat("Heure debut \n")
print(date())

DataDesign.nlin.ph = cDataDesign.PourDiffusion( n             = n                     , 
                                                beta          = coef.col.nlin.ph , 
                                                temps.survie  = c(1,3,5,10)           , 
                                                site          = "colon"               , 
                                                min.adiag     = 1990                  , 
                                                max.adiag     = 2010                  ,
                                                noeud.knot    = 54                    )


ListDataSimulation.nlin.ph = list()

for (i in 1:NbFichier)
{
  print(paste(cat("\n simulation du fichier i "),i ))
  print(Sys.time())
  ListDataSimulation.nlin.ph[[i]] <- cdatasimulation.PourDiffusion( donnees = DataDesign.nlin.ph , beta = coef.col.nlin.ph , modele.simule = DataDesign.nlin.ph$mod.nlin.ph , cens.admin = 2013 , min.adiag = 1990 )
}

#========================save ===================

save( DataDesign.nlin.ph         , file = paste( chemin , chemin2 , "DataDesign.nlin.ph.RData"         , sep = "") )
save( ListDataSimulation.nlin.ph , file = paste( chemin , chemin2 , "ListDataSimulation.nlin.ph.RData" , sep = "") )

print(cat("\n Fin de la simulation des donnees \n"))
print(date())

print(cat("\n SURVIE THEORIQUE DE LA POPULATION \n"))
print(DataDesign.nlin.ph$survie.pop)

heure.fin=Sys.time()
print(heure.fin)

print(cat("\n", "\n", "Durée du programme en heure","\n"))
print(difftime(heure.fin, heure.debut,units = c("hours")))

sink()








