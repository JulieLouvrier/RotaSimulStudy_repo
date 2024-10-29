##Script to run the simulation 

## -----------------------------------------------------------------------------
#library(mipfp)
library(purrr) 
library(unmarked)
library(extraDistr)
library(here)
set.seed(2022)
Nsim <- 100


## -----------------------------------------------------------------------------
source("r/SimulSevenLevelsPredation_FitUnmarkedPr0.6_08_bothPsiA_6LevelsSites.r")

#the function to run the simulation is "SimulFitUnrmakedThreeLevels <- function(psiA, N, J, Nsim, Pr1, Pr2)" 
#running the simulations
PsiA = c(0.5, 0.7)
rep <- c(3,5,10)
site <- c(30, 50, 100, 150, 200, 250)
Pr1 = 0.6
Pr2 = 0.8

for(a in PsiA){
  for(i in site) {
    for(j in rep) {
      results <- SimulFitUnrmakedThreeLevels(psiA = a, N = i, J = j, Nsim = Nsim, Pr1 = Pr1, Pr2 = Pr2)
      saveRDS(results, file = paste0(getwd(),"/output","/results_7Levels_sims_Pr0.6_0.8_PsiA",a,"_",i,"N",j,"J", "Pr1.1_Pr2.1",Nsim, "Sims.RDS" ))
    }
  }
}
