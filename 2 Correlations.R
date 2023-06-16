if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
#source("0 General.R" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/C10 - Gamer Fear/3 Diagnostic Generalization/3 Analysis/main Anx/", .))

requirePackage("lme4")
requirePackage("lmerTest")

#ratings.wide = read_rds("ratings.wide.rds" %>% paste0(path.rds, .)) #or path
#eyes.wide = read_rds("eyes.wide.rds" %>% paste0(path.rds, .)) #or path
#heart.wide = read_rds("heart.wide.rds" %>% paste0(path.rds, .)) #or path
#eda.wide = read_rds("eda.wide.rds" %>% paste0(path.rds, .)) #or path

data.wide = full_join(questionnaires, ratings.wide, by="subject") %>% 
  full_join(eyes.wide, by=c("subject", "SPAI", "STAI")) %>% 
  full_join(heart.wide, by="subject") %>% 
  full_join(eda.wide, by="subject")

# data.wide.first.half = full_join(questionnaires, ratings.first.wide, by="subject") %>% 
#   full_join(eyes.wide, by="subject") %>% 
#   full_join(heart.wide, by="subject") #%>% 
# #full_join(eda.wide, by="subject")

#write_rds(data.wide, "data.wide.rds" %>% paste0(path.rds, .))


# Screening Quality -------------------------------------------------------
questionnaires %>% with(cor.test(Screening_vorher, Screening_nachher, alternative="greater")) %>% correlation_out()
#questionnaires %>% ggplot(aes(x=Screening_vorher, y=Screening_nachher)) + geom_smooth(method="lm", color="black") + geom_point(size=2) + myGgTheme

questionnaires %>% with(cor.test(SPAI, Screening_vorher, alternative="greater")) %>% correlation_out()
#questionnaires %>% ggplot(aes(x=SPAI, y=Screening_vorher)) + geom_smooth(method="lm", color="black") + geom_point(size=2) + myGgTheme

questionnaires %>% with(cor.test(SPAI, Screening_nachher, alternative="greater")) %>% correlation_out()
#questionnaires %>% ggplot(aes(x=SPAI, y=Screening_nachher)) + geom_smooth(method="lm", color="black") + geom_point(size=2) + myGgTheme

questionnaires %>% with(cor.test(SPAI, STAI, alternative="greater")) %>% correlation_out()
#questionnaires %>% ggplot(aes(x=SPAI, y=STAI)) + geom_smooth(method="lm", color="black") + geom_point(size=2) + myGgTheme


# Significance of gradients -----------------------------------------------
data.wide %>% pull(Gen_all_lds) %>% t.test(mu = 0, alternative="greater") %>% apa::t_apa(es_ci=T)
#gradients %>% pull(Gen_all_lds) %>% mean() %>% signif(3) %>% paste0("M = ", .)
data.wide %>% summarise(M = mean(Gen_all_lds, na.rm = TRUE),
                        SD = sd(Gen_all_lds, na.rm = TRUE)) %>% tibble()

# Correlation LDS & SPAI/STAI
data.wide %>% with(cor.test(Gen_all_lds, SPAI, alternative="greater")) %>% correlation_out()  #apa::cor_apa(r_ci=T)
data.wide %>% with(cor.test(Gen_all_lds, STAI, alternative="greater")) %>% correlation_out()

# Correlation Diff & SPAI/STAI
data.wide %>% with(cor.test(Gen_all_diff, SPAI, alternative="greater")) %>% correlation_out()  #apa::cor_apa(r_ci=T)
data.wide %>% with(cor.test(Gen_all_diff, STAI, alternative="greater")) %>% correlation_out()

# Correlation Level & SPAI/STAI
data.wide %>% with(cor.test(Gen_all_level, SPAI, alternative="greater")) %>% correlation_out()  #apa::cor_apa(r_ci=T)
data.wide %>% with(cor.test(Gen_all_level, STAI, alternative="greater")) %>% correlation_out()

# # Significance of gradients for first half -----------------------------------------------
# data.wide.first.half %>% pull(Gen_all_lds) %>% t.test(mu = 0, alternative="greater") %>% apa::t_apa(es_ci=T)
# #gradients %>% pull(Gen_all_lds) %>% mean() %>% signif(3) %>% paste0("M = ", .)
# data.wide.first.half %>% summarise(M = mean(Gen_all_lds, na.rm = TRUE),
#                         SD = sd(Gen_all_lds, na.rm = TRUE)) %>% tibble()
# 
# # Correlation LDS & SPAI/STAI
# data.wide.first.half %>% with(cor.test(Gen_all_lds, SPAI, alternative="greater")) %>% correlation_out()  #apa::cor_apa(r_ci=T)
# data.wide.first.half %>% with(cor.test(Gen_all_lds, STAI, alternative="greater")) %>% correlation_out()
# 
# # Correlation Diff & SPAI/STAI
# data.wide.first.half %>% with(cor.test(Gen_all_diff, SPAI, alternative="greater")) %>% correlation_out()  #apa::cor_apa(r_ci=T)
# data.wide.first.half %>% with(cor.test(Gen_all_diff, STAI, alternative="greater")) %>% correlation_out()
# 
# # Correlation Level & SPAI/STAI
# data.wide.first.half %>% with(cor.test(Gen_all_level, SPAI, alternative="greater")) %>% correlation_out()  #apa::cor_apa(r_ci=T)
# data.wide.first.half %>% with(cor.test(Gen_all_level, STAI, alternative="greater")) %>% correlation_out()


# Replication Anxiety -----------------------------------------------------
#data.wide = read_rds("data.wide.rds" %>% paste0(path.rds, .))
#data.wide %>% select(subject:STAI, Gen_all_dwell:Gen_mn_roiSwitch) %>% View()

#dwell time
data.wide %>% with(cor.test(Gen_eyes_dwell, SPAI, alternative="two.sided")) %>% correlation_out()
data.wide %>% with(cor.test(Gen_mn_dwell, SPAI, alternative="two.sided")) %>% correlation_out()
data.wide %>% with(cor.test(Gen_eyes_dwell.non, SPAI, alternative="two.sided")) %>% correlation_out()
#data.wide %>% with(cor.test(Gen_mn_dwell, Gen_eyes_dwell.non, alternative="less")) %>% correlation_out()
data.wide %>% with(cor.test(Gen_mn_dwell.non, SPAI, alternative="two.sided")) %>% correlation_out()

#this should be practically the same as the "small" eye ANOVA with SPAI as covariate (but results differ!?)
# reg.spai = data.wide %>% gather("dwell.type", "dwell", c("Gen_eyes_dwell", "Gen_mn_dwell", "Gen_eyes_dwell.non", "Gen_mn_dwell.non")) %>% 
#   select(subject:STAI, dwell.type:dwell) %>% 
#   #sample_n(15) %>% #for testing of mutate
#   mutate(#Diagnosticity = ifelse(grepl(".non", dwell.type), "Non-Diagnostic", "Diagnostic") %>% as.factor(),
#          #ROI = ifelse(grepl("_eyes_", dwell.type), "Eyes", "Mouth/Nose") %>% as.factor(),
#          Diagnosticity = ifelse(grepl(".non", dwell.type), -1, 1),
#          ROI = ifelse(grepl("_eyes_", dwell.type), 1, -1),
#          dwell = scale(dwell), SPAI = scale(SPAI)) %>% 
#   lmer(dwell ~ Diagnosticity*ROI*SPAI + (1|subject), .)
#   
# reg.spai %>% summary() %>% print()
# reg.spai %>% lmer.ci()
# #reg.spai %>% lmer.ci(twotailed = F)

#latency
data.wide %>% with(cor.test(Gen_eyes_ms, SPAI, alternative="two.sided")) %>% correlation_out()
data.wide %>% with(cor.test(Gen_mn_ms, SPAI, alternative="two.sided")) %>% correlation_out()
data.wide %>% with(cor.test(Gen_eyes_ms.non, SPAI, alternative="two.sided")) %>% correlation_out()
data.wide %>% with(cor.test(Gen_mn_ms.non, SPAI, alternative="two.sided")) %>% correlation_out()


# Ratings -----------------------------------------------------------------
#data.wide = read_rds("data.wide.rds" %>% paste0(path.rds, .))

#Ratings: LDS & Square Root of Time to diagnostic ROI
reg.lds.ms.sqr = data.wide %>% filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
  gather("lds.type", "lds", c("Gen_eyes_lds", "Gen_mn_lds")) %>% 
  gather("ms.diag.type", "ms.diag", c("Gen_eyes_ms", "Gen_mn_ms")) %>%
  select(subject:STAI, lds.type:ms.diag) %>% 
  mutate(ms.diag = sqrt(ms.diag),
         lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
         ms.diag.type = ifelse(grepl("_eyes_", ms.diag.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
  filter(lds.type == ms.diag.type) %>% 
  mutate(ms.diag.type = ifelse(ms.diag.type=="Eyes", -1, 1)) %>% 
  mutate(lds = scale(lds), ms.diag = scale(ms.diag), SPAI = scale(SPAI), STAI = scale(STAI)) %>% #z-transform
  #lmer(lds ~ ms.diag*ms.diag.type*SPAI + (1|subject), .)
  lmer(lds ~ ms.diag*ms.diag.type*STAI + (1|subject), .)
  #lmer(lds ~ ms.diag*ms.diag.type + (1|subject), .)

reg.lds.ms.sqr %>% summary() %>% print()
reg.lds.ms.sqr %>% lmer.ci() 
#reg.lds.ms.sqr %>% lmer.ci(twotailed = F)

#ms.diag main effect
data.wide %>% with(cor.test(Gen_all_lds, Gen_all_ms, alternative="less")) %>% correlation_out()  #apa::cor_apa(r_ci=T)

print(correl.ms.plot <- data.wide %>% filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
        gather("lds.type", "lds", c("Gen_eyes_lds", "Gen_mn_lds")) %>% 
        gather("ms.diag.type", "ms.diag", c("Gen_eyes_ms", "Gen_mn_ms")) %>%
        select(subject:STAI, lds.type:ms.diag) %>% 
        mutate(lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
               ms.diag.type = ifelse(grepl("_eyes_", ms.diag.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
        filter(lds.type == ms.diag.type) %>% 
        ggplot(aes(x=sqrt(ms.diag), y=lds, color=ms.diag.type, fill=ms.diag.type, shape=ms.diag.type)) +
        geom_smooth(method="lm", size=1.5, alpha = .2) +
        #stat_smooth(method="lm", size=1.5, alpha = .2, fullrange = T) +
        geom_point(size=4, alpha=.8) + 
        ylab("Linear Deviation Score") + xlab(expression("Time to Diagnostic ROI (" * sqrt(ms) * ")")) + labs(color="Diagnostic", fill="Diagnostic", shape="Diagnostic") +
        scale_shape_manual(values=c(16, 15)) +
        myGgTheme +
        theme(
          #aspect.ratio = 1,
          #legend.position = "none",
          legend.position = c(.87, 1-.87))
          #legend.position = c(.5, 1-.87))
)


#Ratings: LDS & Time to diagnostic ROI
reg.lds.ms = data.wide %>% filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
  gather("lds.type", "lds", c("Gen_eyes_lds", "Gen_mn_lds")) %>% 
  gather("ms.diag.type", "ms.diag", c("Gen_eyes_ms", "Gen_mn_ms")) %>%
  select(subject:STAI, lds.type:ms.diag) %>% 
  mutate(lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
         ms.diag.type = ifelse(grepl("_eyes_", ms.diag.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
  filter(lds.type == ms.diag.type) %>% 
  mutate(ms.diag.type = ifelse(ms.diag.type=="Eyes", -1, 1)) %>% 
  mutate(lds = scale(lds), ms.diag = scale(ms.diag), SPAI = scale(SPAI), STAI = scale(STAI)) %>% #z-transform
  #lmer(lds ~ ms.diag*ms.diag.type*SPAI + (1|subject), .)
  #lmer(lds ~ ms.diag*ms.diag.type*STAI + (1|subject), .)
  lmer(lds ~ ms.diag*ms.diag.type + (1|subject), .)

reg.lds.ms %>% summary() %>% print()
reg.lds.ms %>% lmer.ci() 
#reg.lds.ms %>% lmer.ci(twotailed = F)

data.wide %>% filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
  gather("lds.type", "lds", c("Gen_eyes_lds", "Gen_mn_lds")) %>% 
  gather("ms.diag.type", "ms.diag", c("Gen_eyes_ms", "Gen_mn_ms")) %>%
  select(subject:STAI, lds.type:ms.diag) %>% 
  mutate(lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
         ms.diag.type = ifelse(grepl("_eyes_", ms.diag.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
  filter(lds.type == ms.diag.type) %>% 
  ggplot(aes(x=ms.diag, y=lds, color=ms.diag.type, fill=ms.diag.type, shape=ms.diag.type)) +
  geom_smooth(method="lm", size=1.5, alpha = .2) +
  #stat_smooth(method="lm", size=1.5, alpha = .2, fullrange = T) +
  geom_point(size=4, alpha=.8) + 
  ylab("Linear Deviation Score") + xlab("Time to Diagnostic ROI (ms)") + labs(color="Diagnostic", fill="Diagnostic", shape="Diagnostic") +
  scale_shape_manual(values=c(16, 15)) +
  myGgTheme +
  theme(
    #aspect.ratio = 1,
    #legend.position = "none",
    legend.position = c(.87, 1-.87))
    #legend.position = c(.5, 1-.87))



#Ratings: LDS & Fixations
reg.lds.dwell = data.wide %>% filter(subject %in% exclusions.eye.dwell == F) %>% #manual exclusion because of extreme dwell
  gather("lds.type", "lds", c("Gen_eyes_lds", "Gen_mn_lds")) %>% 
  gather("dwell.type", "dwell", c("Gen_eyes_dwell", "Gen_mn_dwell")) %>%
  select(subject:STAI, lds.type:dwell) %>% 
  mutate(lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
         dwell.type = ifelse(grepl("_eyes_", dwell.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
  filter(lds.type == dwell.type) %>% 
  #mutate(dwell.type = factor(dwell.type)) %>% #dummy coding
  mutate(dwell.type = ifelse(dwell.type=="Eyes", -1, 1)) %>% #effect coding
  mutate(lds = scale(lds), dwell = scale(dwell), SPAI = scale(SPAI), STAI = scale(STAI)) %>% #z-transform
  lmer(lds ~ dwell*dwell.type*SPAI + (1|subject), .)
  #lmer(lds ~ dwell*dwell.type*STAI + (1|subject), .)
  #lmer(lds ~ dwell*dwell.type + (1|subject), .)

reg.lds.dwell %>% summary() %>% print()
reg.lds.dwell %>% lmer.ci()
#reg.lds.dwell %>% lmer.ci(twotailed = F)

#dwell*dwell.type
print(correl.dwell.plot <- data.wide %>% filter(subject %in% exclusions.eye.dwell == F) %>% #manual exclusion because of extreme dwell
        gather("lds.type", "lds", c("Gen_eyes_lds", "Gen_mn_lds")) %>% 
        gather("dwell.type", "dwell", c("Gen_eyes_dwell", "Gen_mn_dwell")) %>%
        select(subject:STAI, lds.type:dwell) %>% 
        mutate(lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
               dwell.type = ifelse(grepl("_eyes_", dwell.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
        filter(lds.type == dwell.type) %>% 
        ggplot(aes(x=dwell, y=lds, color=dwell.type, fill=dwell.type, shape=dwell.type)) +
        #geom_line(aes(group=subject), color="black", linetype="dashed") + 
        geom_smooth(method="lm", size=1.5, alpha = .2) +
        geom_point(size=4, alpha=.8) + 
        ylab("Linear Deviation Score") + xlab("Diagnostic Dwell (%)") + labs(color="Diagnostic", fill="Diagnostic", shape="Diagnostic") +
        scale_shape_manual(values=c(16, 15)) +
        myGgTheme +
        theme(
          #legend.position = "none",
          #legend.position = c(.87, .87))
          #legend.position = c(1-.87, 1-.87))
          legend.position = c(1-.87, .87))
)

#STAI*dwell.type
# print(correl.stai.plot <- data.wide %>% gather("lds.type", "lds", c("Gen_eyes_lds", "Gen_mn_lds")) %>% 
#         gather("dwell.type", "dwell", c("Gen_eyes_dwell", "Gen_mn_dwell")) %>%
#         select(subject:STAI, lds.type:dwell) %>% 
#         mutate(lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
#                dwell.type = ifelse(grepl("_eyes_", dwell.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
#         filter(lds.type == dwell.type) %>% 
#         ggplot(aes(x=STAI, y=lds, color=dwell.type, fill=dwell.type)) +
#         #geom_line(aes(group=subject), color="black", linetype="dashed") + 
#         geom_smooth(method="lm", size=1.5, alpha = .2) +
#         geom_point(size=2) +
#         ylab("Linear Deviation Score") + labs(color="Diagnostic", fill="Diagnostic") +
#         myGgTheme +
#         theme(
#           #legend.position = "none",
#           #legend.position = c(.87, .87),
#           legend.position = c(1-.87, 1-.87),
#         )
# )
# data.wide %>% with(cor.test(Gen_mn_lds, STAI, alternative="two.sided")) %>% correlation_out()
# data.wide %>% with(cor.test(Gen_eyes_lds, STAI, alternative="two.sided")) %>% correlation_out()

#cowplot::plot_grid(correl.dwell.plot, correl.ms.plot, ncol=1, labels="auto") #Figure 6 (preprint) or 5 (manuscript)


# Heartrate -----------------------------------------------------------------
#Heartrate: LDS & Fixations
heart.dwell = data.wide %>% filter(subject %in% exclusions.eye.dwell == F) %>% #manual exclusion because of extreme dwell
  gather("lds.type", "lds", c("HR_Gen_eyes_lds", "HR_Gen_mn_lds")) %>% 
  gather("dwell.type", "dwell", c("Gen_eyes_dwell", "Gen_mn_dwell")) %>%
  select(subject:STAI, lds.type:dwell) %>% 
  mutate(lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
         dwell.type = ifelse(grepl("_eyes_", dwell.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
  filter(lds.type == dwell.type) %>% 
  mutate(dwell.type.dummy = ifelse(dwell.type=="Eyes", -1, 1)) %>% 
  mutate(lds.z = scale(lds), dwell.z = scale(dwell), SPAI.z = scale(SPAI), STAI.z = scale(STAI)) #z-transform
reg.heart.dwell = heart.dwell %>% 
  #lmer(lds.z ~ dwell.z*dwell.type.dummy*SPAI.z + (1|subject), .)
  #lmer(lds.z ~ dwell.z*dwell.type.dummy*STAI.z + (1|subject), .)
  lmer(lds.z ~ dwell.z*dwell.type.dummy + (1|subject), .)

reg.heart.dwell %>% summary()
reg.heart.dwell %>% lmer.ci() 
#reg.heart.dwell %>% lmer.ci(twotailed = F)

#main effect dwell
print(reg.dwell.hr.plot <- heart.dwell %>% 
        ggplot(aes(x=dwell, y=lds, color=dwell.type, fill=dwell.type, shape=dwell.type)) +
        #geom_line(aes(group=subject), color="black", linetype="dashed") + 
        geom_smooth(method="lm", size=1.5, alpha = .2) +
        geom_point(size=4, alpha=.8) + 
        ylab("LDS of Heart Change") + xlab("Diagnostic Dwell (%)") + labs(color="Diagnostic", fill="Diagnostic", shape="Diagnostic") +
        scale_shape_manual(values=c(16, 15)) +
        myGgTheme +
        theme(
          #legend.position = "none",
          #legend.position = c(.87, .87),
          legend.position = c(.87, 1-.87),
        )
)

#dwell.type x SPAI
heart.dwell %>% group_by(dwell.type) %>% summarise(cor.test(lds, SPAI) %>% apa::cor_apa(r_ci=T, print=F))
heart.dwell %>% ggplot(aes(x=SPAI, y=lds, color=dwell.type, fill=dwell.type, shape=dwell.type)) +
  facet_wrap(vars(dwell.type)) +
  geom_smooth(method="lm", size=1.5, color="black", alpha = .2) +
  #geom_point(size=4, alpha=.8) +
  geom_point(size=4) + 
  ylab("LDS of Heart Change") + xlab("SPAI") + labs(color="Diagnostic", fill="Diagnostic", shape="Diagnostic") +
  scale_shape_manual(values=c(16, 15)) +
  scale_color_viridis_d() + scale_fill_viridis_d() +
  myGgTheme +
  theme(
    #legend.position = "none",
    #legend.position = c(.87, .87),
    #legend.position = c(.87, 1-.87),
  )

#dwell x SPAI
# SDs = -2:2
# spai.descr = with(questionnaires, data.frame(SPAI=mean(SPAI)+sd(SPAI)*SDs, SPAI.z=SDs))
# reg.heart.spaiSD.predict = spai.descr %>%
#   expand.grid.df(data.frame(dwell=unique(data.wide$Gen_all_dwell %>% na.omit())) %>% mutate(dwell.z=scale(dwell))) %>%
#   expand.grid.df(data.frame(subject=0, dwell.type.dummy=0)) %>% #include mockup subject for population-level prediction (could avoid this and add "re.form=NA" to predict but less flexibility with additional subject-level predictions or based on non-lmer model objects, i.e. predict.merMod)
#   mutate(lds = predict(reg.heart.dwell, newdata=., allow.new.levels=T)) #predict values from dwell time & z-scaled SPAI
# 
# print(reg.dwell.hr.spai.plot <- heart.dwell %>% group_by(subject) %>% 
#         summarise(dwell=mean(dwell, na.rm=T), lds=mean(lds, na.rm=T), SPAI=mean(SPAI, na.rm=T)) %>% 
#         ggplot(aes(x=dwell, y=lds, color=SPAI, group=SPAI)) +
#         geom_hline(yintercept = 0, linetype="dashed") +
#         geom_line(data=reg.heart.spaiSD.predict, size=4) +
#         geom_point(size=2) +
#         ylab("LDS of Heart Change") + xlab("Diagnostic Dwell (%)") +
#         myGgTheme + scale_color_viridis_c()
# )


#Heartrate: LDS & Square Root of Time to diagnostic ROI
heart.ms = data.wide %>% filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
  gather("lds.type", "lds", c("HR_Gen_eyes_lds", "HR_Gen_mn_lds")) %>% 
  gather("ms.diag.type", "ms.diag", c("Gen_eyes_ms", "Gen_mn_ms")) %>%
  select(subject:STAI, lds.type:ms.diag) %>% 
  mutate(ms.diag = sqrt(ms.diag),
         lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
         ms.diag.type = ifelse(grepl("_eyes_", ms.diag.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
  filter(lds.type == ms.diag.type) %>% 
  mutate(ms.diag.type = ifelse(ms.diag.type=="Eyes", -1, 1)) %>% 
  mutate(lds.z = scale(lds), ms.diag.z = scale(ms.diag), SPAI.z = scale(SPAI), STAI.z = scale(STAI)) #z-transform
reg.heart.ms = heart.ms %>% 
  lmer(lds.z ~ ms.diag.z*ms.diag.type*SPAI.z + (1|subject), .)
  #lmer(lds.z ~ ms.diag.z*ms.diag.type*STAI.z + (1|subject), .)
  #lmer(lds.z ~ ms.diag.z*ms.diag.type + (1|subject), .)

reg.heart.ms %>% summary() %>% print()
reg.heart.ms %>% lmer.ci()
#reg.heart.ms %>% lmer.ci(twotailed = F)

#SPAI x ms.diag (sqrt)
SDs = -2:2
spai.descr = with(questionnaires, data.frame(SPAI=mean(SPAI)+sd(SPAI)*SDs, SPAI.z=SDs))
reg.heart.ms.spaiSD.predict = spai.descr %>%
  expand.grid.df(data.frame(ms.diag=unique(data.wide$Gen_all_ms %>% na.omit() %>% sqrt())) %>% mutate(ms.diag.z=scale(ms.diag))) %>%
  expand.grid.df(data.frame(subject=0, ms.diag.type=0)) %>% #include mockup subject for population-level prediction (could avoid this and add "re.form=NA" to predict but less flexibility with additional subject-level predictions or based on non-lmer model objects, i.e. predict.merMod)
  mutate(lds = predict(reg.heart.ms, newdata=., allow.new.levels=T)) #predict values from dwell time & z-scaled SPAI

print(reg.ms.hr.spai.plot <- heart.ms %>% group_by(subject) %>% 
        summarise(ms.diag=mean(ms.diag, na.rm=T), lds=mean(lds, na.rm=T), SPAI=mean(SPAI, na.rm=T)) %>% 
        ggplot(aes(x=ms.diag, y=lds, color=SPAI, group=SPAI)) +
        geom_line(data=reg.heart.ms.spaiSD.predict, size=4) +
        geom_hline(yintercept = 0, linetype="dashed") +
        geom_point(size=2) +
        ylab("LDS of Heart Change") + xlab(expression("Time to Diagnostic ROI (" * sqrt(ms) * ")")) + 
        myGgTheme + scale_color_viridis_c()
)


#SPAI x ms.diag.type (n.s. in dwell analysis despite same dummy variable -> suppressed/enhanced by other predictors)
print(data.wide %>% filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
        gather("lds.type", "lds", c("HR_Gen_eyes_lds", "HR_Gen_mn_lds")) %>%
        gather("ms.diag.type", "ms.diag", c("Gen_eyes_ms", "Gen_mn_ms")) %>%
        select(subject:STAI, lds.type:ms.diag) %>% 
        mutate(lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
               ms.diag.type = ifelse(grepl("_eyes_", ms.diag.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>%
        filter(lds.type == ms.diag.type) %>%
        ggplot(aes(x=SPAI, y=lds, color=ms.diag.type, fill=ms.diag.type)) +
        #geom_line(aes(group=subject), color="black", linetype="dashed") +
        geom_smooth(method="lm", size=1.5, alpha = .2) +
        geom_point(size=2) +
        ylab("LDS of Heart Change") + labs(color="Diagnostic", fill="Diagnostic") +
        myGgTheme +
        theme(
          #legend.position = "none",
          #legend.position = c(.87, .87),
          legend.position = c(1-.87, 1-.87),
        )
)
data.wide %>% with(cor.test(HR_Gen_eyes_lds, SPAI, alternative="two.sided")) %>% correlation_out()
data.wide %>% with(cor.test(HR_Gen_mn_lds, SPAI, alternative="two.sided")) %>% correlation_out()



# EDA --------------------------------------------------
reg.eda.dwell = data.wide %>% gather("lds.type", "lds", c("EDA_Gen_eyes_lds", "EDA_Gen_mn_lds")) %>% 
  gather("dwell.type", "dwell", c("Gen_eyes_dwell", "Gen_mn_dwell")) %>%
  select(subject:STAI, lds.type:dwell) %>% 
  mutate(lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
         dwell.type = ifelse(grepl("_eyes_", dwell.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
  filter(lds.type == dwell.type) %>% 
  #mutate(dwell.type = factor(dwell.type)) %>% 
  mutate(dwell.type = ifelse(dwell.type=="Eyes", -1, 1)) %>% 
  mutate(lds = scale(lds), dwell = scale(dwell), SPAI = scale(SPAI), STAI = scale(STAI)) %>% #z-transform
  #lmer(lds ~ dwell*dwell.type*SPAI + (1|subject), .) 
  lmer(lds ~ dwell*dwell.type*STAI + (1|subject), .) 
  #lmer(lds ~ dwell*dwell.type + (1|subject), .) 

reg.eda.dwell %>% summary() %>% print()
reg.eda.dwell %>% lmer.ci()
#reg.eda.dwell %>% lmer.ci(twotailed = F)

data.wide %>% gather("lds.type", "lds", c("EDA_Gen_eyes_lds", "EDA_Gen_mn_lds")) %>% 
  gather("dwell.type", "dwell", c("Gen_eyes_dwell", "Gen_mn_dwell")) %>%
  select(lds.type:dwell) %>% 
  mutate(lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
         dwell.type = ifelse(grepl("_eyes_", dwell.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
  filter(lds.type == dwell.type) %>% 
  ggplot(aes(x=dwell, y=lds, color=dwell.type, fill=dwell.type)) +
  geom_smooth(method="lm", size=1.5, alpha = .2) +
  geom_point(size=2) +
  ylab("Linear Deviation Score (EDA)") + xlab("Diagnostic Dwell (%)") + labs(color="Diagnostic", fill="Diagnostic") +
  myGgTheme + theme(
    #legend.position = "none",
    legend.position = c(1-.87, .87)
  )


#EDA: LDS & Time to diagnostic ROI
reg.eda.ms = data.wide %>% gather("lds.type", "lds", c("EDA_Gen_eyes_lds", "EDA_Gen_mn_lds")) %>% 
  gather("ms.diag.type", "ms.diag", c("Gen_eyes_ms", "Gen_mn_ms")) %>%
  select(subject:STAI, lds.type:ms.diag) %>% 
  mutate(ms.diag = sqrt(ms.diag),
         lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
         ms.diag.type = ifelse(grepl("_eyes_", ms.diag.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
  na.omit() %>% 
  filter(lds.type == ms.diag.type) %>% 
  mutate(ms.diag.type = ifelse(ms.diag.type=="Eyes", -1, 1)) %>% 
  mutate(lds = scale(lds), ms.diag = scale(ms.diag), SPAI = scale(SPAI), STAI = scale(STAI)) %>% #z-transform
  #lmer(lds ~ ms.diag*ms.diag.type*SPAI + (1|subject), .) 
  lmer(lds ~ ms.diag*ms.diag.type*STAI + (1|subject), .) 
  #lmer(lds ~ ms.diag*ms.diag.type + (1|subject), .) 

reg.eda.ms %>% summary() %>% print()
reg.eda.ms %>% lmer.ci()
reg.eda.ms %>% lmer.ci(twotailed = F)

data.wide %>% gather("lds.type", "lds", c("EDA_Gen_eyes_lds", "EDA_Gen_mn_lds")) %>% 
  gather("ms.diag.type", "ms.diag", c("Gen_eyes_ms", "Gen_mn_ms")) %>%
  select(lds.type:ms.diag) %>% 
  mutate(lds.type = ifelse(grepl("_eyes_", lds.type), "Eyes", "Mouth/Nose") %>% as.factor(),
         ms.diag.type = ifelse(grepl("_eyes_", ms.diag.type), "Eyes", "Mouth/Nose") %>% as.factor()) %>% 
  filter(lds.type == ms.diag.type) %>% 
  ggplot(aes(x=sqrt(ms.diag), y=lds, color=ms.diag.type, fill=ms.diag.type)) +
  geom_smooth(method="lm", size=1.5, alpha = .2) +
  #stat_smooth(method="lm", size=1.5, alpha = .2, fullrange = T) +
  geom_point(size=2) +
  ylab("Linear Deviation Score (EDA)") + xlab(expression("Time to Diagnostic ROI (" * sqrt(ms) * ")")) + labs(color="Diagnostic", fill="Diagnostic") +
  myGgTheme + theme(
    #legend.position = "none",
    legend.position = c(.87, .87))
