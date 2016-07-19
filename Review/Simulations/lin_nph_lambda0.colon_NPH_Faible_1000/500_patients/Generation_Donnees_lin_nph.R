#-------------------------------------------------------------------------------------------------------------
# Generation of the data 
#-------------------------------------------------------------------------------------------------------------

# Paths to change !!

chemin  <- c("C:/Users/coraline.danieli/Desktop/HCL/V/2679/Review/")
chemin2 <- c("Simulations/lin_nph_lambda0.colon_NPH_Faible_1000/500_patients/")

sink( paste( chemin , chemin2 , "Generation_Donnees_lin_nph.lis" , sep = "") )

library(statmod)

set.seed(1351)

print(system("hostname", intern=T))
heure.debut=Sys.time()
print(heure.debut)

source( paste( chemin , "Simulations/Fonctions/cDataDesign.PourDiffusion.R"     , sep = "") )
source( paste( chemin , "Simulations/Fonctions/cdatasimulation.PourDiffusion.R" , sep = "") )
source( paste( chemin , "Simulations/Fonctions/ListModeles.R"                   , sep = "") )

###### ----- > Expected mortality hazard

tauxatt = dget( file = paste( chemin , "Simulations/MUA/muaDF8933.dat" , sep = "") )

###### ----- > Parameter for the simulation 

load( file = paste( chemin , "Simulations/Parametres_Theoriques/colon_TxBase/parametre_theorique/coef.colon.lin.nph.RData"      , sep = "" ) )
load( file = paste( chemin , "Simulations/Parametres_Theoriques/lin_nph_estomac/parametre_theorique/coef.estomac.lin.nph.Rdata" , sep = "" ) )

coef.estomac.lin.nph[1:5] <- coef.colon.lin.nph[1:5]

NbFichier = 1000
n         = 500

rm( mod.lin.nph )

mod.lin.nph = list.model.simu()$ms.lin.nph

###### ----- > Simulation de données

# DataDesign.lin.nph : Contains all the constant information 
# ListDataSimulation : Generation of time to death du to cancer and time to death due to other causes

cat("Heure debut \n")
print(date())

DataDesign.lin.nph = cDataDesign.PourDiffusion( n             = n                    , 
                                                beta          = coef.estomac.lin.nph , 
                                                modele.simule = mod.lin.nph          , 
                                                temps.survie  = c(1,3,5,10)          , 
                                                site          = "colon"              , 
                                                min.adiag     = 1990                 , 
                                                max.adiag     = 2010                 )

ListDataSimulation.lin.nph = list()

for (i in 1:NbFichier)
{
  print(paste(cat("\n simulation du fichier i "),i ))
  print(Sys.time())
  ListDataSimulation.lin.nph[[i]] <- cdatasimulation.PourDiffusion( donnees = DataDesign.lin.nph , beta = coef.estomac.lin.nph , modele.simule = mod.lin.nph , cens.admin = 2013 , min.adiag = 1990 )
}

#========================save ===================

save( DataDesign.lin.nph         , file = paste( chemin , chemin2 , "DataDesign.lin.nph.RData"         , sep = "") )
save( ListDataSimulation.lin.nph , file = paste( chemin , chemin2 , "ListDataSimulation.lin.nph.RData" , sep = "") )

print(cat("\n Fin de la simulation des donnees \n"))
print(date())

print(cat("\n SURVIE THEORIQUE DE LA POPULATION \n"))
print(DataDesign.lin.nph$survie.pop)

heure.fin=Sys.time()
print(heure.fin)

print(cat("\n", "\n", "Durée du programme en heure","\n"))
print(difftime(heure.fin, heure.debut,units = c("hours")))

sink()








