if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
requirePackage("openintro")
requirePackage("mediation")
requirePackage("lm.beta")

## Mediation via attention to diagnostic regions

# direct effect x to m (path a)
path_a <-lm(Gen_all_ms~SPAI,data.wide)
summary(path_a)
lm.beta(path_a)

# direct effect x to y (path c) and effect m to y (path b)
path_b_c <-lm(Gen_all_lds~Gen_all_ms+SPAI,data.wide)
summary(path_b_c)
lm.beta(path_b_c)

# mediation (which part of the effect x on y runs through mediator, i.e. path a and b)

results <- mediate(path_a, path_b_c, 
                   treat = "SPAI", mediator = "Gen_all_ms", 
                   boot = TRUE)
summary(results)
# ACME = average causal mediation effect (indirect effect, path a and b)
# ADE = average direct effect (direct effect controlling for effect of mediator)
# total effect = ACME + ADE


## Mediation via attention to non-diagnostic regions

# direct effect x to m (path a)
path_a <-lm(Gen_all_ms.non~SPAI,data.wide)
summary(path_a)
lm.beta(path_a)

# direct effect x to y (path c) and effect m to y (path b)
path_b_c <-lm(Gen_all_lds~Gen_all_ms.non+SPAI,data.wide)
summary(path_b_c)
lm.beta(path_b_c)

# mediation (which part of the effect x on y runs through mediator, i.e. path a and b)

results <- mediate(path_a, path_b_c, 
                   treat = "SPAI", mediator = "Gen_all_ms.non", 
                   boot = TRUE)
summary(results)


