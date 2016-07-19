#-------------------------------------------------------------------------------------------------------------
# Generation of the data 
#-------------------------------------------------------------------------------------------------------------

# Paths to change !!

chemin  <- c("C:/Users/coraline.danieli/Desktop/HCL/V/2679/Review/")
chemin2 <- c("Simulations/lin_ph_lambda0.colon_0.02_2000/2000_patients/")

sink( paste( chemin , chemin2 , "Generation_Donnees_lin_ph.lis" , sep = "") )

library(statmod)

set.seed(1000)

print(system("hostname", intern=T))
heure.debut=Sys.time()
print(heure.debut)

source( paste( chemin , "Simulations/Fonctions/cDataDesign.PourDiffusion.R"     , sep = "") )
source( paste( chemin , "Simulations/Fonctions/cdatasimulation.PourDiffusion.R" , sep = "") )
source( paste( chemin , "Simulations/Fonctions/ListModeles.R"                   , sep = "") )

###### ----- > Expected mortality hazard

tauxatt = dget( file = paste( chemin , "Simulations/MUA/muaDF8933.dat" , sep = "") )

###### ----- > Parameter for the simulation 

load( file = paste( chemin , "Simulations/Parametres_Theoriques/lin_ph_colon/parametre_theorique/coef.colon.lin.ph.Rdata" , sep = "" ) )

coef.colon.lin.ph[6] <- 0.02

NbFichier = 2000
n         = 2000

rm( mod.lin.ph )

mod.lin.ph = list.model.simu()$ms.lin.ph

###### ----- > Simulation de données

# DataDesign.lin.nph : Contains all the constant information 
# ListDataSimulation : Generation of time to death du to cancer and time to death due to other causes

cat("Heure debut \n")
print(date())

DataDesign.lin.ph = cDataDesign.PourDiffusion(  n             = n                    , 
                                                beta          = coef.colon.lin.ph    , 
                                                modele.simule = mod.lin.ph           , 
                                                temps.survie  = c(1,3,5,10)          , 
                                                site          = "colon"              , 
                                                min.adiag     = 1990                 , 
                                                max.adiag     = 2010                 )

ListDataSimulation.lin.ph = list()

for (i in 1:NbFichier)
{
  print(paste(cat("\n simulation du fichier i "),i ))
  print(Sys.time())
  ListDataSimulation.lin.ph[[i]] <- cdatasimulation.PourDiffusion( donnees = DataDesign.lin.ph , beta = coef.colon.lin.ph , modele.simule = mod.lin.ph , cens.admin = 2013 , min.adiag = 1990 )
}

#========================save ===================

save( DataDesign.lin.ph         , file = paste( chemin , chemin2 , "DataDesign.lin.ph.RData"         , sep = "") )
save( ListDataSimulation.lin.ph , file = paste( chemin , chemin2 , "ListDataSimulation.lin.ph.RData" , sep = "") )

print(cat("\n Fin de la simulation des donnees \n"))
print(date())

print(cat("\n SURVIE THEORIQUE DE LA POPULATION \n"))
print(DataDesign.lin.ph$survie.pop)

heure.fin=Sys.time()
print(heure.fin)

print(cat("\n", "\n", "Durée du programme en heure","\n"))
print(difftime(heure.fin, heure.debut,units = c("hours")))

sink()








