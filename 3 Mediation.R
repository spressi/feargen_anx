if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
#requirePackage("openintro", load=F)
requirePackage("mediation", load=F)
#requirePackage("lme4")
#requirePackage("lm.beta")

data.wide.z <- data.wide %>%
  mutate(Gen_all_ms.z = scale(Gen_all_ms)[, 1],
         Gen_all_ms.non.z = scale(Gen_all_ms.non)[, 1],
         Gen_all_lds.z = scale(Gen_all_lds)[, 1],
         SPAI.z = scale(SPAI)[, 1],
         STAI.z = scale(STAI)[, 1],)

## Mediation via attention to diagnostic regions

# direct effect x to m (path a)
path_a <-lm(Gen_all_ms.z~SPAI.z,data.wide.z)
summary(path_a)
confint(path_a, level = 0.95)
#lm.beta(path_a)

# direct effect x to y (path c) and effect m to y (path b)
path_b_c <-lm(Gen_all_lds.z~Gen_all_ms.z+SPAI.z,data.wide.z)
summary(path_b_c)
confint(path_b_c, level = 0.95)
#lm.beta(path_b_c)

# mediation (which part of the effect x on y runs through mediator, i.e. path a and b)
set.seed(3681253)
results <- mediation::mediate(path_a, path_b_c, 
                              treat = "SPAI.z", mediator = "Gen_all_ms.z", 
                              boot = TRUE)
summary(results)
# ACME = average causal mediation effect (indirect effect, path a and b)
# ADE = average direct effect (direct effect controlling for effect of mediator)
# total effect = ACME + ADE


## Mediation via attention to non-diagnostic regions

# direct effect x to m (path a)
path_a <-lm(Gen_all_ms.non.z ~SPAI.z,data.wide.z)
summary(path_a)
confint(path_a, level = 0.95)
#lm.beta(path_a)

# direct effect x to y (path c) and effect m to y (path b)
path_b_c <-lm(Gen_all_lds.z~Gen_all_ms.non.z+SPAI.z,data.wide.z)
summary(path_b_c)
confint(path_b_c, level = 0.95)
#lm.beta(path_b_c)

# mediation (which part of the effect x on y runs through mediator, i.e. path a and b)
set.seed(7894237)
results <- mediation::mediate(path_a, path_b_c, 
                              treat = "SPAI.z", mediator = "Gen_all_ms.non.z", 
                              boot = TRUE)
summary(results)


## Mediation via attention to diagnostic regions with STAI ------------------------------

# direct effect x to m (path a)
path_a <-lm(Gen_all_ms.z~STAI.z,data.wide.z)
summary(path_a)
confint(path_a, level = 0.95)
#lm.beta(path_a)

# direct effect x to y (path c) and effect m to y (path b)
path_b_c <-lm(Gen_all_lds.z~Gen_all_ms.z+STAI.z,data.wide.z)
summary(path_b_c)
confint(path_b_c, level = 0.95)
#lm.beta(path_b_c)

# mediation (which part of the effect x on y runs through mediator, i.e. path a and b)
set.seed(34993253)
results <- mediation::mediate(path_a, path_b_c, 
                              treat = "STAI.z", mediator = "Gen_all_ms.z", 
                              boot = TRUE)
summary(results)
# ACME = average causal mediation effect (indirect effect, path a and b)
# ADE = average direct effect (direct effect controlling for effect of mediator)
# total effect = ACME + ADE


## Mediation via attention to non-diagnostic regions

# direct effect x to m (path a)
path_a <-lm(Gen_all_ms.non.z ~STAI.z,data.wide.z)
summary(path_a)
confint(path_a, level = 0.95)
#lm.beta(path_a)

# direct effect x to y (path c) and effect m to y (path b)
path_b_c <-lm(Gen_all_lds.z~Gen_all_ms.non.z+STAI.z,data.wide.z)
summary(path_b_c)
confint(path_b_c, level = 0.95)
#lm.beta(path_b_c)

# mediation (which part of the effect x on y runs through mediator, i.e. path a and b)
set.seed(5677612)
results <- mediation::mediate(path_a, path_b_c, 
                              treat = "STAI.z", mediator = "Gen_all_ms.non.z", 
                              boot = TRUE)
