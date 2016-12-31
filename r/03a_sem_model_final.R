###### Load Libraries for analysis
library(nlme)
library(lsmeans)
library(piecewiseSEM)
library(car)

###### Load data
source("./02_data_load_merge.R")
exp_data_func <- as.data.frame(exp_data_func)

###some useful functions for piecewiseSEM use
resp <- function(i) Reduce(paste, deparse(formula(i)[[2]]))


sem.anova <- function(modlist){
  data.frame(purrr::map_df(modlist, ~{
    a <- as.data.frame(Anova(.))
    a$predictor <- rownames(a)
    a$response <- resp(.)
    a %>%
      select(response, predictor, Chisq, Df, `Pr(>Chisq)`) %>%
      rename(p_value = `Pr(>Chisq)`) %>%
      mutate(Chisq = round(Chisq,4), p_value = round(p_value, 4))
  }))
}

status_lsmeans <- function(x, ...){
  r <- resp(x)
  ret <- lsmeans(x, list(pairwise ~ status), ...)
  names(ret) <- r
  ret
}

######################### THIS IS WHAT I RAN FOR THE MODEL
exp_data_func$rd <- exp_data_func$`RNA/DNA`
exp_data_func$belowC <- exp_data_func$`below%C`
exp_data_func$belowN <- exp_data_func$`below%N`
gen_div_status <- lme(observed_otus~  status, 
                      random =~ 1|Genotype, data=exp_data_func, method="ML")
gen_activity_status <- lme(rd ~ status + observed_otus, 
                           random =~ 1|Genotype, data=exp_data_func, method="ML")
gen_meta_status <- lme(`Metabolism`~ status+ observed_otus, 
                       random =~ 1|Genotype, data=exp_data_func, method="ML")
phen_mod_status <- lme(belowgallic_uM ~ rd + observed_otus + Metabolism + status, 
                       random =~ 1|Genotype, data=exp_data_func, method="ML")
c_mod_status <- lme(belowC ~ rd + observed_otus + Metabolism + status + belowgallic_uM +
                      belowbiomass_g + abovebiomass_g, 
                     random =~ 1|Genotype, data=exp_data_func, method="ML")
n_mod_status <- lme(belowN ~ rd + observed_otus + Metabolism + status +
                      belowbiomass_g + abovebiomass_g, 
                     random =~ 1|Genotype, data=exp_data_func, method="ML")

biomass_mod_status <- lme(belowbiomass_g ~ rd + belowgallic_uM +  
                            observed_otus + Metabolism + status, 
                          random =~ 1|Genotype, data=exp_data_func, method="ML")
Abiomass_mod_status <- lme(abovebiomass_g ~ rd + belowgallic_uM + 
                             observed_otus + Metabolism + status, 
                           random =~ 1|Genotype, data=exp_data_func, method="ML")

sem_mod_nlme <- list(
  gen_div_status,
  gen_activity_status,
  gen_meta_status,
  phen_mod_status,
  c_mod_status,
  n_mod_status,
  biomass_mod_status,
  Abiomass_mod_status
  
)

#get SEM fit information
sem.fit(sem_mod_nlme, data=exp_data_func, 
        corr.errors=c("belowbiomass_g~~abovebiomass_g", "belowN ~~ belowC"))
        
#Get coefficients from the fit model
sem.coefs(sem_mod_nlme, data=exp_data_func, 
          corr.errors=c("belowbiomass_g~~abovebiomass_g", "belowN ~~ belowC"),
          intercept=T)

#Evaluate Chi Square tests of parameter significance
sem.anova(sem_mod_nlme)
write.csv(sem.anova(sem_mod_nlme), "./sem_chisq.csv", row.names=FALSE)

ph_tests <- lapply(sem_mod_nlme, status_lsmeans, adjust="none")
ph_comp <- lapply(ph_tests, function(x) x[[2]])
names(ph_comp) <- sapply(ph_tests, function(x) names(x)[1])


ph_comp_tab <- do.call(rbind, lapply(ph_comp, function(x) as.data.frame(print(x))))
write.csv(ph_comp_tab, "./post_hocs.csv", row.names=TRUE)

#### Get coefficients, including categorical means
sem_mod_coefs <- lapply(sem_mod_nlme, function(x) update(x, . ~ . - 1))

ctab <- sem.coefs(sem_mod_coefs, data=exp_data_func, 
                  corr.errors=c("belowbiomass_g~~abovebiomass_g", "belowN ~~ belowC"),
                      intercept=T)
ctab <- data.frame(ctab) %>%
  arrange(response, predictor)
write.csv(ctab, "./sem_coefs.csv", row.names=FALSE)



std_ctab <- sem.coefs(sem_mod_coefs, data=exp_data_func, 
                      corr.errors=c("belowbiomass_g~~abovebiomass_g", "belowN ~~ belowC"),
                      standardize="range", intercept=T)
std_ctab <- data.frame(std_ctab) %>%
  arrange(response, predictor)
write.csv(std_ctab, "./sem_std_coefs.csv", row.names=FALSE)
