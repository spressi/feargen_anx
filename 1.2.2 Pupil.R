if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
#source("0 General.R" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/C10 - Gamer Fear/3 Diagnostic Generalization/3 Analysis/main no fMRI/", .))

#exclusions.eye.num
sample.rate.pupil = 100

pupilToNum = function(pupilFile, path.sep="/", file.ext="\\.") {
  pupilFile %>% pathToCode(path.sep, file.ext) %>% gsub("[^0-9]", "", .) %>% as.integer()
}

files.pupil = list.files(path.pupil, pattern = ".csv", full.names = T)
pupil.avg = list() #vector("list", length(files.pupil))
pupil.trial = list()
for (filename in files.pupil) {
  #filename = files.pupil %>% sample(1)
  df = read.csv2(filename) %>% group_by(trial) %>% nest() %>%
    merge(readRDS("ratings.rds" %>% paste0(path.rds, .)) %>% filter(subject==filename %>% pupilToNum()) %>% select(trial, shock, threat, pair, diagnostic, sex, phase, threat_num, threat_both), by.x = "trial") %>% 
    mutate(data = data %>% map(function(x) {
      signal = x %>% as.numeric()
      signal = signal - mean(signal[1:sample.rate.pupil]) #1 sec baseline correction
      data.frame(mmChange = signal, 
                 change.z = scale(signal, center=F), #trial-level standardization (cf. Greenberg et al. 2013)
                 time = 0:(length(signal)-1)/sample.rate.pupil - 1)
    }))
  
  #subject.avg = df %>% group_by(phase, threat) %>%
  subject.avg = df %>% group_by(phase, threat, diagnostic, sex) %>% 
    summarise(signal = data %>% bind_rows() %>% group_by(time) %>% summarise(
      mmChange.se = se(mmChange, na.rm=T), mmChange = mean(mmChange, na.rm=T),
      change.z.se = se(change.z, na.rm=T), change.z = mean(change.z, na.rm=T)))
  
  pupil.avg[[filename %>% pathToCode()]] = subject.avg
  pupil.trial[[filename %>% pathToCode()]] = df
}


# Exploration -------------------------------------------------------------------
pupil = pupil.avg %>% bind_rows(.id="subject") %>% mutate(time = signal$time, mmChange = signal$mmChange, change.z = signal$change.z) %>% select(-signal) %>% 
  group_by(phase, threat, time, subject) %>% summarise(mmChange = mean(mmChange, na.rm=T), change.z = mean(change.z, na.rm=T)) #collapse across unwanted grouping variables

pupil.diff = pupil %>% pivot_wider(names_from = threat, values_from = mmChange, id_cols = c("phase", "time", "subject")) %>% 
  mutate(`diff.CS+` = `CS+` - `CS-`, `diff.GS1` = `GS1` - `CS-`, `diff.GS2` = `GS2` - `CS-`, `diff.GS3` = `GS3` - `CS-`, `diff.GS4` = `GS4` - `CS-`) %>%
  select(subject, phase, time, starts_with("diff")) %>% 
  pivot_longer(cols = starts_with("diff."), names_to = "threat", values_to = "diff") %>% 
  na.omit() %>% mutate(threat = threat %>% gsub("diff.", "", .)) %>% 
  merge(pupil, ., by=c("subject", "phase", "time", "threat"), all.x=T)

pupil.diff %>% #filter(time <= 4) %>%
  group_by(phase, threat, time) %>% 
  summarise(mmChange.se = se(mmChange, na.rm=T), mmChange = mean(mmChange, na.rm=T), 
            diff.se = se(diff, na.rm=T), diff = mean(diff, na.rm=T)) %>% 
  ggplot(aes(x=time, y=mmChange, color=threat, group=threat)) + 
  facet_wrap(vars(phase)) +
  #geom_ribbon(aes(ymin=mmChange-mmChange.se*1.96, ymax=mmChange+mmChange.se*1.96, fill=threat), color=NA, alpha=.1) + #replace 1.96 with qnorm(.975)?
  geom_ribbon(aes(ymin=mmChange-diff.se*1.96, ymax=mmChange+diff.se*1.96, fill=threat), color=NA, alpha=.1) + #within-CIs relative to CS-
  #geom_ribbon(aes(ymin=mmChange-diff, ymax=mmChange+diff, fill=threat), color=NA, alpha=.1) + #sanity check of within-CI calculation (CI should exactly touch reference wave)
  geom_line(size=2) +
  scale_color_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_shape_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_fill_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + 
  ylab("Pupil Size Change (mm)") + xlab("Trial Time (sec)") + labs(color="Threat", fill="Threat") +
  theme_bw() + theme(
    #aspect.ratio = 1,
    panel.background = element_rect(fill="white", color="white"),
    legend.background = element_rect(fill="white", color="grey"),
    legend.key=element_rect(fill='white'),
    axis.text = element_text(color="black"),
    axis.ticks.x = element_line(color="black"),
    axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black"),
    legend.text = element_text(size=14, color="black"),
    legend.title = element_text(size=14, color="black"),
    strip.text.x = element_text(size=12, color="black"),
    axis.text.x = element_text(size=16, color="black"),
    axis.text.y = element_text(size=16, color="black"),
    axis.title = element_text(size=16, color="black"))


#grand averages by diagnostic region
pupil.diag = pupil.avg %>% bind_rows(.id="subject") %>% mutate(time = signal$time, mmChange = signal$mmChange) %>% select(-signal) %>% 
  group_by(phase, diagnostic, threat, time) %>% summarise(mmChange.se = se(mmChange, na.rm=T), 
                                                          mmChange = mean(mmChange, na.rm=T))
pupil.diag %>% ggplot(aes(x=time, y=mmChange, color=threat, group=threat)) + 
  facet_grid(rows=vars(diagnostic), cols=vars(phase)) +
  geom_ribbon(aes(ymin=mmChange-mmChange.se*1.96, ymax=mmChange+mmChange.se*1.96, fill=threat), color=NA, alpha=.1) +
  geom_line(size=2) +
  scale_color_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_shape_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_fill_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + 
  ylab("Pupil Size Change (mm)") + xlab("Trial Time (sec)") + labs(color="Threat", fill="Threat") +
  theme_bw() + theme(
    #aspect.ratio = 1,
    panel.background = element_rect(fill="white", color="white"),
    legend.background = element_rect(fill="white", color="grey"),
    legend.key=element_rect(fill='white'),
    axis.text = element_text(color="black"),
    axis.ticks.x = element_line(color="black"),
    axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black"),
    legend.text = element_text(size=14, color="black"),
    legend.title = element_text(size=14, color="black"),
    strip.text.x = element_text(size=12, color="black"),
    axis.text.x = element_text(size=16, color="black"),
    axis.text.y = element_text(size=16, color="black"),
    axis.title = element_text(size=16, color="black"))


# ANOVA -------------------------------------------------------------------
pupil.trial.avg = pupil.trial %>% bind_rows(.id="subject") %>% group_by(subject, trial, phase, threat, diagnostic, sex) %>% 
  summarise(pupil = data %>% bind_rows() %>% filter(time > .5, time <= 4) %>% summarise(
    mmChange.se = se(mmChange, na.rm=T), mmChange = mean(mmChange, na.rm=T),
    change.z.se = se(change.z, na.rm=T), change.z = mean(change.z, na.rm=T))) %>% 
  mutate(subject = subject %>% gsub("Look", "", .) %>% gsub("_pupil", "", .) %>% as.numeric(), #subject to number
         mmChange = pupil$mmChange, mmChange.se = pupil$mmChange.se, mmChange.z = pupil$change.z, mmChange.z.se = pupil$change.z.se) %>% select(-pupil) 
#pupil.trial.avg %>% write_rds("pupil.trial.avg.rds" %>% paste0(path.rds, .))

pupil.subject = pupil.avg %>% bind_rows(.id="subject") %>% mutate(time = signal$time, mmChange = signal$mmChange, change.z = signal$change.z) %>% select(-signal) %>% 
  filter(time > .5, time <= 4) %>% group_by(subject, phase, threat, diagnostic, sex) %>% 
  summarise(mmChange = mean(mmChange, na.rm=T), change.z = mean(change.z, na.rm=T)) %>% 
  group_by(subject, phase) %>% mutate(mmChange.z = scale(mmChange))

pupil.subject %>% filter(phase=="PreAcq") %>% group_by(subject, threat) %>% summarise(mmChange = mean(mmChange, na.rm=T)) %>% 
  t.test(mmChange ~ threat, ., paired=T) %>% apa::t_apa(es_ci=T)
pupil.subject %>% filter(phase=="Acq") %>% group_by(subject, threat) %>% summarise(mmChange = mean(mmChange, na.rm=T)) %>% 
  t.test(mmChange ~ threat, ., paired=T) %>% apa::t_apa(es_ci=T)

pupil.subject %>% filter(phase=="Acq") %>% group_by(threat, subject) %>% summarise(mmChange = mean(mmChange, na.rm=T)) %>% 
  summarise(mmChange.se = se(mmChange, na.rm=T), mmChange = mean(mmChange, na.rm=T))
  

anova.pupil = pupil.subject %>% filter(phase=="Gen") %>% 
  ez::ezANOVA(dv=.(mmChange), wid=.(subject),
              #within=.(threat), within_full=.(diagnostic),
              within=.(threat, diagnostic),
              detailed=T, type=2)
anova.pupil %>% apa::anova_apa()
anova.pupil %>% ez.ci()

# t-tests of threat-levels vs. CS-
pupil.ga.gen.subj = pupil.subject %>% filter(phase=="Gen") %>% group_by(threat, subject) %>% summarise(mmChange = mean(mmChange, na.rm=T)) %>% 
  mutate(threat_num = threat %>% as.numeric())
for (i in 2:6) {
  levels = c(1, i)
  cat(paste0("\n\nComparing levels: ", paste(levels, collapse=" vs. "), "\n"))
  pupil.ga.gen.subj %>% filter(threat_num %in% levels) %>% 
    t.test(mmChange ~ threat, ., alternative="less", 
           paired=T) %>% apa::t_apa(es_ci=T) #schoRsch::t_out()
}

#mean scores generalization phase
pupil.ga.gen = pupil.ga.gen.subj %>% group_by(threat) %>% 
  summarise(mmChange.se = se(mmChange, na.rm=T), mmChange = mean(mmChange, na.rm=T))
print(pupil.gradient.plot <- pupil.ga.gen %>% ggplot(aes(x=threat, y=mmChange, color=threat, group=NA)) + 
        #geom_dotplot(data=pupil.ga.gen.subj, mapping=aes(group=threat, fill=threat), binaxis="y", alpha=.25, color="black", stackratio=1, stackdir="centerwhole", dotsize=.5) +
        geom_point() + geom_path(data=pupil.ga.gen %>% filter(threat %in% c("CS-", "CS+")), color = "black", size=1.5) + #generalization line (geom_point first for order of x-axis)
        geom_line(size=1) + geom_point(size=4.5) + 
        #geom_errorbar(aes(ymin=mmChange-mmChange.se*1.96, ymax=mmChange+mmChange.se*1.96), size=1.5) +
        geom_errorbar(aes(ymin=mmChange-mmChange.se, ymax=mmChange+mmChange.se), size=1.5) +
        scale_color_manual(values=colors, guide=guide_legend(reverse=T)) +
        #scale_fill_manual(values=colors, guide=guide_legend(reverse=T)) +
        scale_fill_manual(values=rep("grey", 6), guide=guide_legend(reverse=T)) +
        scale_x_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) +
        ylab("Pupil Size Change (mm)") + xlab("Threat") + labs(color="Threat") +
        theme_bw() + theme(
          #aspect.ratio = 1,
          legend.position = "none",
          panel.background = element_rect(fill="white", color="white"),
          legend.background = element_rect(fill="white", color="grey"),
          legend.key=element_rect(fill='white'),
          legend.text = element_text(size=14, color="black"),
          legend.title = element_text(size=14, color="black"),
          axis.text = element_text(color="black"),
          axis.text.x = element_text(size=16, color="black"),
          axis.text.y = element_text(size=16, color="black"),
          strip.text.x = element_text(size=12, color="black"),
          axis.ticks.x = element_line(color="black"),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          axis.title = element_text(size=16, color="black"),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, "pt"))))
#pupil.gradient.plot %>% ggsave("5.2 Pupil Generalization.png", ., device="png", path=paste0(path.rds, "../../4 Results/1 no MRI/"), width=1920, height=1080, units="px")

#pupil grand average for generalization & only between errorbars
print(pupil.grandAverage.plot <- pupil %>% filter(phase=="Gen") %>% #filter(time <= 4) %>%
        group_by(phase, threat, time, subject) %>% summarise(mmChange = mean(mmChange, na.rm=T)) %>% 
        summarise(mmChange.se = se(mmChange, na.rm=T), mmChange = mean(mmChange, na.rm=T)) %>% 
        ggplot(aes(x=time, y=mmChange, color=threat, group=threat)) + 
        #geom_rect(xmin=.5, xmax=4, ymin=-Inf, ymax=Inf, fill="grey", colour=NA) +
        annotate("rect", xmin=.5, xmax=4, ymin=-Inf, ymax=Inf, fill="black", alpha=.2, colour=NA) +
        geom_hline(yintercept=0) +
        #geom_ribbon(aes(ymin=mmChange-mmChange.se*1.96, ymax=mmChange+mmChange.se*1.96, fill=threat), color=NA, alpha=.1) + #replace 1.96 with qnorm(.975)?
        geom_ribbon(aes(ymin=mmChange-mmChange.se, ymax=mmChange+mmChange.se, fill=threat), color=NA, alpha=.1) + #replace 1.96 with qnorm(.975)?
        geom_line(size=2) +
        scale_color_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_shape_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_fill_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + 
        guides(colour=guide_legend(reverse=T), fill=guide_legend(reverse=T)) + 
        ylab("Pupil Size Change (mm)") + xlab("Trial Time (sec)") + labs(color="Threat", fill="Threat") +
        theme_bw() + theme(
          #aspect.ratio = 1,
          panel.background = element_rect(fill="white", color="white"),
          legend.background = element_rect(fill="white", color="grey"),
          legend.key=element_rect(fill='white'),
          axis.text = element_text(color="black"),
          axis.ticks.x = element_line(color="black"),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          legend.text = element_text(size=14, color="black"),
          legend.title = element_text(size=14, color="black"),
          strip.text.x = element_text(size=12, color="black"),
          axis.text.x = element_text(size=16, color="black"),
          axis.text.y = element_text(size=16, color="black"),
          axis.title = element_text(size=16, color="black"),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, "pt"))))
#pupil.grandAverage.plot %>% ggsave("5.1 Pupil Generalization over trial time.png", ., device="png", path=paste0(path.rds, "../../4 Results/1 no MRI/"), width=1920, height=1080, units="px")
#cowplot::plot_grid(pupil.grandAverage.plot, pupil.gradient.plot, ncol=1, labels="auto") %>% ggsave("Figure 3. Pupil.png", ., device="png", path=paste0(path.rds, "../../5 Output/1 Paper - Fear Generalization x Attention/"), width=16.26, height=21, units="cm") #Figure 3

# Reliability -------------------------------------------------------------
pupil.reliability = pupil.trial.avg %>% filter(phase == "Gen") %>% select(subject, threat, diagnostic, mmChange) %>% 
  group_by(subject, threat, diagnostic) %>% mutate(trial = 1:n()) %>% arrange(subject, trial, threat, diagnostic)
alpha.pupil = tibble()
for (t in pupil.reliability %>% pull(threat) %>% unique() %>% sort()) {
  for (d in pupil.reliability %>% pull(diagnostic) %>% unique() %>% sort()) {
    alpha.temp = pupil.reliability %>% filter(threat == t, diagnostic == d) %>% 
      pivot_wider(names_from = trial, names_prefix = "trial", values_from = mmChange) %>% ungroup() %>% select(contains("trial")) %>% 
      ltm::cronbach.alpha(CI=T, na.rm=T)
    cat(t, " & ", d, ": ", alpha.temp$alpha %>% signif(3), " [", paste0(alpha.temp$ci %>% signif(3), collapse=", "), "]\n", sep="")
    alpha.pupil = alpha.pupil %>% bind_rows(tibble(Threat = t, Diagnostic = d, alpha = alpha.temp$alpha, alpha.low = min(alpha.temp$ci), alpha.high = max(alpha.temp$ci)))
  }
}
alpha.pupil %>% mutate(alpha.z = alpha %>% psych::fisherz()) %>% 
  summarise(alpha.m = alpha.z %>% mean() %>% psych::fisherz2r(), alpha.se = alpha.z %>% sd() %>% {. / sqrt(n())} %>% psych::fisherz2r())

# Gradient Analysis -------------------------------------------------------
individualPlots = F
pupil.subject.gen = pupil.subject %>% filter(phase=="Gen")

#subject-level means generalization phase
pupil.gradients.simple = data.frame(subject=character(), lds=numeric(), diff=numeric(), level=numeric())
for (s in pupil.subject.gen$subject %>% unique()) {
  pupil.temp = pupil.subject.gen %>% filter(subject==s) %>% 
    group_by(threat, subject) %>% summarise(mmChange.se = se(mmChange, na.rm=T), mmChange = mean(mmChange, na.rm=T),
                                            mmChange.z.se = se(mmChange.z, na.rm=T), mmChange.z = mean(mmChange.z, na.rm=T))
  
  gradient = gradient.analysis(pupil.temp$mmChange)
  gradient.z = gradient.analysis(pupil.temp$mmChange.z)
  pupil.gradients.simple = bind_rows(pupil.gradients.simple, data.frame(subject=s, lds=gradient[1], diff=gradient[2], level=gradient[3], lds.z=gradient.z[1], diff.z=gradient.z[2], level.z=gradient.z[3]))
  
  if (individualPlots) {
    gradient = gradient %>% signif(3)
    print(pupil.temp %>% ggplot(aes(x=threat, y=mmChange, color=threat)) + #generalization line
            geom_point() + #do this first to establish correct order of x-axis
            geom_path(data=pupil.temp %>% filter(threat %in% c("CS-", "CS+")), aes(group=NA), color = "black", size=1.5) +
            geom_point() + scale_color_manual(values=colors) + 
            geom_errorbar(aes(ymin=mmChange-mmChange.se*1.96, ymax=mmChange+mmChange.se*1.96)) +
            ggtitle(paste0(s, ": lds = ", gradient[1], ", diff = ", gradient[2], ", level = ", gradient[3])))
  }
}
names(pupil.gradients.simple)[-1] = "Gen_all_" %>% paste0(names(pupil.gradients.simple)[-1])

#subject-level means generalization phase by diagnostic
pupil.gradients.diagnostic = data.frame()
for (s in pupil.subject.gen$subject %>% unique()) {
  pupil.temp = pupil.subject.gen %>% filter(subject==s)
  
  gradient.eyes = "Eyes" %>% {filter(pupil.temp, diagnostic==.)} %>% .$mmChange %>% gradient.analysis()
  gradient.mn = "Mouth/Nose" %>% {filter(pupil.temp, diagnostic==.)} %>% .$mmChange %>% gradient.analysis()
  gradient.eyes.z = "Eyes" %>% {filter(pupil.temp, diagnostic==.)} %>% .$mmChange.z %>% gradient.analysis()
  gradient.mn.z = "Mouth/Nose" %>% {filter(pupil.temp, diagnostic==.)} %>% .$mmChange.z %>% gradient.analysis()
  pupil.gradients.diagnostic = bind_rows(pupil.gradients.diagnostic, 
                                         data.frame(subject=s, 
                                                    Gen_eyes_lds = gradient.eyes[1], Gen_eyes_diff = gradient.eyes[2], Gen_eyes_level = gradient.eyes[3],
                                                    Gen_mn_lds = gradient.mn[1], Gen_mn_diff = gradient.mn[2], Gen_mn_level = gradient.mn[3], 
                                                    Gen_eyes_lds.z = gradient.eyes.z[1], Gen_eyes_diff.z = gradient.eyes.z[2], Gen_eyes_level.z = gradient.eyes.z[3],
                                                    Gen_mn_lds.z = gradient.mn.z[1], Gen_mn_diff.z = gradient.mn.z[2], Gen_mn_level.z = gradient.mn.z[3]))
  
  if (individualPlots) {
    print(pupil.temp %>% ggplot(aes(x=threat, y=mmChange, color=diagnostic, shape=diagnostic)) + 
            #geom_path(data=ratings.ga.m.diagnostic %>% filter(threat %in% c(1, 6)), aes(group=NA), color = "black", size=1.5) + #generalization line
            geom_point(position=dodge) + 
            #geom_errorbar(aes(ymin=rating.m-rating.se*1.96, ymax=rating.m+rating.se*1.96), position=dodge, width=dodge.width) + 
            ggtitle(s))
  }
}
#names(pupil.gradients.diagnostic) = c("lds", "diff", "level") %>% c(., paste0(., ".z")) %>% {c("subject", paste0("Gen_eyes_", .), paste0("Gen_mn_", .))}
pupil.gradients = merge(pupil.gradients.simple, pupil.gradients.diagnostic, by="subject") %>% 
  mutate(subject = subject %>% pupilToNum())
names(pupil.gradients)[-1] = "Pup_" %>% paste0(names(pupil.gradients)[-1])

#saveRDS(pupil.gradients, "pupil.rds" %>% paste0(path.rds, .))

# Significance of gradients -----------------------------------------------
pupil.gradients %>% pull(Pup_Gen_all_lds) %>% t.test(mu = 0, alternative="greater") %>% apa::t_apa(es_ci=T)
#pupil.gradients %>% pull(Pup_Gen_all_lds) %>% mean() %>% signif(3) %>% paste0("M = ", .)
pupil.gradients %>% summarise(M = Pup_Gen_all_lds %>% mean(), SD = Pup_Gen_all_lds %>% sd()) %>% tibble()
