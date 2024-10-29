SimulFitUnrmakedThreeLevels <- function(psiA, N, J, Nsim, Pr1, Pr2){ 
  
  # function to run the simul and fit them - static occupancy model w/ 2 species
  # N is the number of sites, J the number of occasions, and Nsim the number of simulations
  
  runsimuls <- function(true_psi){
    psiA <- psiA
    psiBA <- true_psi[2]
    psiBa <- true_psi[3]
    psiB <- psiA *psiBA + psiBa * (1 - psiA)
    psi11 <- psiBA * psiA # psiAB = psiBA * psiA
    psi10 <- psiA - psi11
    psi01 <- psiB - psi11
    psi00 <- 1 - (psi10 + psi01 + psi11)
    #probability that the species is at least detected once Pr1 : 0.6 and 0.8 
    # so then we extract the value of the general detection probability as Pr1 = 1-(1-p)^j
    
    ## Set nb surveys and sites----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    nsurveys <- J
    nsites <- N
    
    p1 <- 1-(1-Pr1)^(1/nsurveys) # detection prob species A
    p2 <- 1-(1-Pr2)^(1/nsurveys) # detection prob species B
    
    
    ## Simulate data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # true occupancy states; z_cat = 1 if both species present, 2 if only species 1, 3 if only species 2, and 4 if both absent
    z_cat <- extraDistr::rcat(nsites, c(psi11, psi10, psi01, psi00)) # 11, 10, 01, 00
    # species detections and non-detections
    y <- list()
    y[[1]] <- matrix(NA, nsites, nsurveys)
    y[[2]] <- matrix(NA, nsites, nsurveys)
    for (k in 1:nsurveys){
      y[[1]][,k] <- rbinom(nsites, 1, p1 * ((z_cat == 1) | (z_cat == 2)))
      y[[2]][,k] <- rbinom(nsites, 1, p2 * ((z_cat == 1) | (z_cat == 3)))
    }
    names(y) <- c('sp1','sp2')
    
    
    ## Fit both models with and without interactions----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    data <- unmarkedFrameOccuMulti(y = y)
    occFormulas_mod1 <- c('~1','~1','~0') # with no interactions
    occFormulas_mod2 <- c('~1','~1','~1') # with interactions
    detFormulas <- c('~1','~1')
    mod1 <- occuMulti(detFormulas,occFormulas_mod1,data,se=T)
    mod2 <- occuMulti(detFormulas,occFormulas_mod2,data,se=T)
    
    
    ## Get estimated quantities from model 1 with no interections:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    mod1_pA <- predict(mod1,'det',species=1)[1,1] # detection spA
    mod1_pA_CI <- predict(mod1,'det',species=1)[1,c(3:4)]
    
    mod1_pB <- predict(mod1,'det',species=2)[1,1] # detection spB
    mod1_pB_CI <- predict(mod1,'det',species=2)[1,c(3:4)]
    
    mod1_psiA <- predict(mod1,'state',species=1)[1,1] # psiA
    mod1_psiA_CI <- predict(mod1,'state',species=1)[1,c(3:4)]
    
    mod1_psiB <- predict(mod1,'state',species=2)[1,1] # psiB
    mod1_psiB_CI <- predict(mod1,'state',species=2)[1,c(3:4)]
    
    mod1_psiBA <- predict(mod1,'state',species='sp2',cond='sp1')[1,1] # psiBA
    mod1_psiBA_CI <- predict(mod1,'state',species='sp2',cond='sp1')[1,c(3:4)]
    
    mod1_psiBa <- predict(mod1,'state',species='sp2',cond='-sp1')[1,1] # psiBa
    mod1_psiBa_CI <- predict(mod1,'state',species='sp2',cond='-sp1')[1,c(3:4)]
    
    
    ## Estimated quantities from model 2 with interactions ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    mod2_pA <- predict(mod2,'det',species=1)[1,1] # detection spA
    mod2_pA_CI <- predict(mod2,'det',species=1)[1,c(3:4)]
    
    mod2_pB <- predict(mod2,'det',species=2)[1,1] # detection spB
    mod2_pB_CI <- predict(mod2,'det',species=2)[1,c(3:4)]
    
    mod2_psiA <- predict(mod2,'state',species=1)[1,1] # psiA
    mod2_psiA_CI <- predict(mod2,'state',species=1)[1,c(3:4)]
    
    mod2_psiB <- predict(mod2,'state',species=2)[1,1] # psiB
    mod2_psiB_CI <- predict(mod2,'state',species=2)[1,c(3:4)]
    
    mod2_psiBA <- predict(mod2,'state',species='sp2',cond='sp1')[1,1] # psiBA
    mod2_psiBA_CI <- predict(mod2,'state',species='sp2',cond='sp1')[1,c(3:4)]
    
    mod2_psiBa <- predict(mod2,'state',species='sp2',cond='-sp1')[1,1] # psiBa
    mod2_psiBa_CI <- predict(mod2,'state',species='sp2',cond='-sp1')[1,c(3:4)]
    
    
    
    ## Compare results----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    res <- data.frame(param = c('psiA', 
                                'psiBA', 
                                'psiBa',
                                'pA',
                                'pB'),
                      true = c(psiA, 
                               psiBA,  
                               psiBa, 
                               Pr1, 
                               Pr2), 
                      mod_without = c(mod1_psiA, 
                                      mod1_psiBA, 
                                      mod1_psiBa, 
                                      mod1_pA, 
                                      mod1_pB),
                      mod_without_lower = c(mod1_psiA_CI$lower,
                                            mod1_psiBA_CI$lower,
                                            mod1_psiBa_CI$lower,
                                            mod1_pA_CI$lower,
                                            mod1_pB_CI$lower), 
                      mod_without_upper = c(mod1_psiA_CI$upper,
                                            mod1_psiBA_CI$upper,
                                            mod1_psiBa_CI$upper,
                                            mod1_pA_CI$upper,
                                            mod1_pB_CI$upper),
                      mod_with = c(mod2_psiA, 
                                   mod2_psiBA, 
                                   mod2_psiBa, 
                                   mod2_pA, 
                                   mod2_pB),
                      mod_with_lower = c(mod2_psiA_CI$lower,
                                         mod2_psiBA_CI$lower,
                                         mod2_psiBa_CI$lower,
                                         mod2_pA_CI$lower,
                                         mod2_pB_CI$lower),
                      mod_with_upper = c(mod2_psiA_CI$upper,
                                         mod2_psiBA_CI$upper,
                                         mod2_psiBa_CI$upper,
                                         mod2_pA_CI$upper,
                                         mod2_pB_CI$upper))
    
    params <- data.frame( 
      AIC_without = mod1@AIC,
      negLogLike_without = mod1@negLogLike,
      nparams_without = length(mod1@opt$par),
      AIC_with = mod2@AIC,
      negLogLike_with = mod2@negLogLike,
      nparams_with = length(mod2@opt$par))
    
    freq <- table(z_cat)/nsites #frequencies of 11, 10, 01, 00 
    simul_psiA = freq[2] + freq[1] #psiA = psi10 + psi11
    simul_psiB = freq[3] + freq[1] #psiB = psi01 + psi11
    simul_psiBA = freq[1] / simul_psiA # psiBA = psi11/psiA
    simul_psiBa = (simul_psiB - (simul_psiA*simul_psiBA))/(1 - simul_psiA)
    
    simuls <- data.frame(
      freq = c(simul_psiA, simul_psiB, simul_psiBa, simul_psiBA),
      true = c(psiA, psiB, psiBa, psiBA),
      names = c("psiA", "psiB", "psiBa", "psiBA")
    )
    
    estimatesmode2 <- mod2@estimates
    ListResults <- list()
    
    ListResults[[1]] <- res
    ListResults[[2]] <- params
    ListResults[[3]] <- simuls
    ListResults[[4]] <- estimatesmode2
    
    return(ListResults)
    
  }
  
  ## Set occupancy and detection parameters ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # occupancy prob of species A e.g. predator
  # occupancy CONDITIONAL prob of species B e.g. prey: 
  # - psiBA is prob B present GIVEN A present
  # - psiBa is prob B present GIVEN A absent
  
  
  ################################################################################
  ## No level of avoidance ## ----------------------------------------------------
  ################################################################################
  
  true_psi_none <- c(psiA,.6,.6) # psiA, psiBA, psiBa
  
  ################################################################################
  ## Low level of avoidance ## ---------------------------------------------------
  ################################################################################
  
  true_psi_low1 <- c(psiA,.5,.6) # psiA, psiBA, psiBa
  true_psi_low2 <- c(psiA,.4,.6) # psiA, psiBA, psiBa
  true_psi_low3 <- c(psiA,.3,.6) # psiA, psiBA, psiBa
  
  true_psi_low4 <- c(psiA,.2,.6) 
  true_psi_low5 <- c(psiA,.1,.6)
  ################################################################################
  ## High level of avoidance ## --------------------------------------------------
  ################################################################################
  
  true_psi_high <- c(psiA,.1,.8) # psiA, psiBA, psiBa
  
  
  resultsnone <- list()
  resultslow1 <- list()
  resultslow2 <- list()
  resultslow3 <- list()
  resultslow4 <- list()
  resultslow5 <- list()
  resultshigh <- list()
  
  for (i in 1:Nsim) {
    resultsnone[[i]] <- tryCatch({runsimuls(true_psi_none)}, error = function(e){})
    resultslow1[[i]] <- tryCatch({runsimuls(true_psi_low1)}, error = function(e){})
    resultslow2[[i]] <- tryCatch({runsimuls(true_psi_low2)}, error = function(e){})
    resultslow3[[i]] <- tryCatch({runsimuls(true_psi_low3)}, error = function(e){})
    resultslow4[[i]] <- tryCatch({runsimuls(true_psi_low4)}, error = function(e){})
    resultslow5[[i]] <- tryCatch({runsimuls(true_psi_low5)}, error = function(e){})
    resultshigh[[i]] <- tryCatch({runsimuls(true_psi_high)}, error = function(e){})
    print(i)
  }
  
  
  listSevenlevels <- list(resultsnone, resultslow1, resultslow2, resultslow3, 
                          resultslow4, resultslow5, resultshigh)
  
  return(listSevenlevels)
  
}


