if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
#source("0 General.R" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/C10 - Gamer Fear/3 Diagnostic Generalization/3 Analysis/main Anx/", .))

trials.n.discrimination = 20
path.discrimination = "Data/discrimination/" %>% paste0(path, .) #discrimination task

discrim = list.files(path.discrimination, pattern = ".csv", full.names = T) %>% 
  lapply(read_tsv) %>% lapply(function(x) x %>% mutate(subject = subject %>% as.character())) %>% 
  bind_rows() %>% select(-expName) %>% select(-block) %>% 
  mutate(picNum = picName %>% gsub("stim\\", "", ., fixed=T) %>% str_sub(end=1),
         diagnostic = if_else(picNum %in% c(1, 2, 5, 6), "eyes", "mouth/nose")) %>% 
  mutate(subject = subject %>% gsub("vp", "", .) %>% as.integer()) %>% arrange(subject) %>% 
  filter(subject %in% exclusions == F)

discrim.counts = discrim %>% count(subject, pairs)
discrim.counts %>% filter(n %% trials.n.discrimination != 0) #subjects that didn't finish each block
discrim.counts %>% pull(subject) %>% {function(x) x[x %>% duplicated()]}(.) #subjects with both pairs


# Last 20 trials of each participant --------------------------------------
discrim.filt = discrim %>% mutate(.by = subject, #for each subject
                                  trial = 1:n(), #recode trial number
                                  trialhelper = trial - max(trial) + trials.n.discrimination) %>% 
  filter(trialhelper > 0) %>% #for each subject, only keep last 20 trials
  select(-trialhelper)

# discrim.counts2 = discrim.filt %>% count(subject, pairs)
# discrim.counts2 %>% filter(n %% trials.n.discrimination != 0) #subjects that didn't finish each block
# discrim.counts2 %>% pull(subject) %>% {function(x) x[x %>% duplicated()]}(.) #subjects with both pairs

discrim.correct = discrim.filt %>% summarize(.by = c(subject, diagnostic, pairs),
                                             correct = mean(correct),
                                             n = n()) %>% arrange(correct)
discrim.problem = discrim.correct %>% filter(correct < .7) %>% pull(subject) %>% unique() %>% sort()

#discrim %>% filter(subject %in% discrim.problem) %>% View() #looks fine

#discrim.correct = discrim.correct %>% filter(subject %in% discrim.problem == F)

discrim.correct %>% ggplot(aes(x = correct, fill = diagnostic)) +
  geom_histogram(color = "black")

with(discrim.correct %>% pivot_wider(names_from = diagnostic, values_from = correct), 
     t.test(eyes, `mouth/nose`, paired=T)) %>% apa::t_apa(es_ci=T)
discrim.correct %>% summarize(.by = diagnostic, M = mean(correct), SD = sd(correct))


# Anxiety -----------------------------------------------------------------
discrim.correct = left_join(discrim.correct, questionnaires) %>% 
  mutate(correct = scale(correct), SPAI = scale(SPAI), STAI = scale(STAI), #z-transform
         diagnostic = ifelse(diagnostic=="eyes", -1, 1))

requirePackage("lme4")
requirePackage("lmerTest")
discrim.correct %>% lmer(correct ~ SPAI*diagnostic + (1|subject), .) %>% summary() %>% print()
discrim.correct %>% lmer(correct ~ STAI*diagnostic + (1|subject), .) %>% summary() %>% print()
