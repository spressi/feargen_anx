if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
#source("0 General.R" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/C10 - Gamer Fear/3 Diagnostic Generalization/3 Analysis/main Anx/", .))

exclusions.hr = c(
  40:43, #FFP2 mask (interacts with breathing more heavily)
  54, #temporary exclusion until workaround for missing ratings :D
  76, #too many ectopic beats
  95 #too many ectopic beats
) %>% c(exclusions) %>% unique() %>% sort()

minhr =  50 #minimum plausible heart rate
maxhr = 120 #maximum plausible heart rate

scaling.window = c(-1, seq(0, 5.5, by=.5)) # Scoring bins in seconds (real time scaling; may be non-integer)

scaleHR = function(hr_t, hr, st, en) {
  wsum = 0 #weighted sum (will be divided by valid ranges below)
  naRanges = 0
  i = which(hr_t > st)[1] #first r peak after marker start (might be NA)
  while (!is.na(i) && i <= length(hr_t) && hr_t[i] < en) { #build up weighted sum until last r peak in interval (exclusively! see below)
    #determine range
    lower = ifelse(i > 1 && hr_t[i-1] > st, hr_t[i-1], st) #find lower end of scoring interval (might be st if previous beat does not exist or out or scoring range)
    range = hr_t[i] - lower
    
    #determine value
    if (i > 1 && !is.na(hr[i])) value = hr[i]
    else {
      value = 0
      naRanges = naRanges + range #lower total range by invalid ranges
    }
    
    wsum = wsum + range * value
    i = i + 1
  } #wsum and naRanges are computed except for last beat
  
  #last beat
  if (!is.na(i) && i > 1 && i <= length(hr_t)) {
    range = en - max(c(hr_t[i-1], st)) #it might happen that there is only one beat before st and one after en => use hr[i] (implicitly done by setting range = en - st via max(c(...)))
    
    value = hr[i]
    if (is.na(value)) {
      value = 0
      naRanges = naRanges + range #lower total range by invalid ranges
    } else value = hr[i]
    
    wsum = wsum + range * value
  }
  
  result = wsum / (en - st - naRanges)
  if (is.na(result) || result ==  0) result = NA
  return(result)
}


# Read & Score HR --------------------------------------------------------------------
vpn.ecg.rpeaks = list.files(path.rpeaks, pattern=path.rpeaks.postfix, full.names=TRUE)
#ratings.all = read_rds("ratings.rds" %>% paste0(path.rds, .))
#rawfiles = vpn.ecg.rpeaks %>% gsub("rpeaks/", "", .) %>% gsub(path.rpeaks.postfix, "", .) %>% paste0(".txt")

hr.list = list() #vector("list", length(vpn.ecg.rpeaks))
for (vpi in seq(vpn.ecg.rpeaks)) {
  vp = vpn.ecg.rpeaks[vpi]
  if (vp %>% pathToCode() %>% codeToNum() %in% exclusions.hr) next
  print(vp %>% pathToCode())
  
  code = vp %>% pathToCode() %>% pathToCode(file.ext = "_")
  ratingfile = ratings.all %>% filter(subject == code %>% codeToNum())
  
  #load r peaks
  allrpeak = read.csv2(vp)[, 1] #only first column
  hr = 60/diff(allrpeak) #convert to bpm
  
  #load triggers
  trigger = read.phys(paste0(path.phys, code, ".txt")) %>% .$Trigger %>% recode.triggers()
  
  if (exclusions.phys.trials[[code]] %>% is.null() == F) { #exclude trials manually (cp. triggerCheck)
    triggers.only = trigger[trigger != 0]
    toExclude = exclusions.phys.trials[[code]]
    triggers.only[toExclude] = 0
    trigger[trigger != 0] = triggers.only
  }
  triggers.n = sum(trigger!=0)
  
  timeline = seq(trigger) / sample.rate
  marker = timeline[{trigger != 0} %>% which() %>% tail(trials.n)]
  
  missingBegin = min(marker) + min(scaling.window); missingBegin = missingBegin[missingBegin <= 0] # <= 0 because first sample point is 1 / sample.rate, i.e. 0 is out of range
  missingEnd = max(timeline) - (max(marker) + max(scaling.window)); missingEnd = missingEnd[missingEnd < 0] # < 0 because if difference exactly 0, then last sample point in rage
  
  #check if data is missing at the edges of data
  if (!is_empty(missingBegin)) warning(paste(code, ": Lacking data at BEGINNING", -missingBegin, "sec"))
  if (!is_empty(missingEnd)) warning(paste(code, ": Lacking data at END", -missingEnd, "sec"))
  
  #check for plausiblilty and issue warning
  if (min(hr) < minhr || max(hr) > maxhr) warning(paste(code, ": Implausible heart rate", hr %>% min() %>% round(1), "-", hr %>% max() %>% round(1)))
  
  hr = c(NA, hr) #matching allrpeak[i] to hr[i] (heart rate values only valid if PREVIOUS r peak exists)
  
  # Real time scoring
  allhr = numeric()
  for (trial in seq(marker)) {
    mtime = marker[trial]
    hr_t = allrpeak - mtime #heart rate time (relative to marker)
    
    # Real time scaling
    hrtrial = numeric()
    for (j in seq(scaling.window)[-1]) { #for all indices except for the first (due to j-1 indexing)
      current = ifelse(mtime + scaling.window[j] < 0, NA, #skip marker time points that refer to negative times (i.e. out of data)
                       scaleHR(hr_t, hr, scaling.window[j-1], scaling.window[j]))
      hrtrial = c(hrtrial, current)
    }
    
    allhr = rbind(allhr, hrtrial)
  }
  
  deltaval = allhr - matrix(allhr[, 1], nrow=nrow(allhr), ncol=ncol(allhr))
  ratings.conditions = ratingfile$condition %>% tail(marker %>% length())
  shocks = ratingfile$shock == "True" %>% tail(n=trials.n)
  out = data.frame(trial = 1:trials.n, condition = ratings.conditions, 
                   shock = shocks, shockPrior = c(FALSE, lag(shocks)[-1]),
                   hrbl = allhr[, 1], hr = deltaval[, 2:ncol(deltaval)])
  
  hr.list[[code]] = out
  #write.csv2(out, paste(savepath.ecg, code,"_task.csv",sep=""), row.names=FALSE, quote=FALSE)
}

#list to one giant dataframe
heart.wide = hr.list %>% bind_rows(.id="subject") %>% 
  mutate(subject = subject %>% gsub("\\D", "", .) %>% as.integer()) %>% tibble()
rm(hr.list); row.names(heart.wide) = NULL

heart = heart.wide %>% subset(subject %in% exclusions.hr == F) %>% gather(key="time", value="HRchange", matches("hr\\.\\d+")) %>% tibble() %>% 
  left_join(questionnaires, by="subject") %>% 
  mutate(threat = condition %% 10, 
         pair = condition %/% 10,
         diagnostic = as.factor(ifelse(pair %% 2 == 0, "mouth/nose", "eyes")),
         sex = as.factor(ifelse(pair > 2, "female", "male")),
         phase = c("PreAcq", "Acq", "Gen") %>% cut(trial, breaks=c(0, preAcqEnd, acqEnd, Inf), labels=.),
         pair = as.factor(pair),
         pairs = ifelse(pair %in% c(1, 4), 1, 2) %>% as.factor(), 
         threat = as.factor(ifelse(pair %in% 2:3 & threat %in% 2:5, 7 - threat, threat)),
         time = time %>% gsub("hr.", "", .) %>% as.integer() %>% {. / 2} %>% as.factor(),
         SPAI.z = scale(SPAI)[,1], STAI.z = scale(STAI)[,1]) %>% select(-condition)

#all(heart == read_rds("heart.rds" %>% paste0(path.rds, .)), na.rm=T) #check equivalence of processing
#heart %>% write_rds("heart.rds" %>% paste0(path.rds, .))

# Plots -------------------------------------------------------------
#heart = read_rds("heart.rds" %>% paste0(path.rds, .))

#plot hr change over trial time
# heart.ga.gen.timeOnly = heart %>% filter(phase == "Gen") %>% group_by(time) %>% 
#   summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T)) %>% 
#   mutate(time = time %>% as.character() %>% as.numeric()) %>% 
#   bind_rows(data.frame(time=0, HRchange=0, HRchange.se=0)) #add origin
# print(heart.trialtimeOnly.plot <- heart.ga.gen.timeOnly %>% ggplot(aes(x=time, y=HRchange)) + 
#         geom_hline(yintercept = 0, linetype="dashed") +
#         geom_ribbon(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96), color=NA, alpha=.1) +
#         #geom_errorbar(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96)) +
#         geom_point(size=3) + geom_line() + 
#         ylab("Heart Rate Change (bpm)") + xlab("Trial Time (sec)") + myGgTheme)

#plot hr change over trial time with factor threat
heart.ga.gen.time = heart %>% filter(phase == "Gen") %>% group_by(threat, time) %>% 
  summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T)) %>% 
  mutate(time = time %>% as.character() %>% as.numeric()) %>% 
  bind_rows(data.frame(threat=unique(heart$threat), time=0, HRchange=0, HRchange.se=0)) #add origin
print(heart.trialtime.plot <- heart.ga.gen.time %>% ggplot(aes(x=time, y=HRchange, color=threat, group=threat, shape=threat)) + 
        geom_hline(yintercept = 0, linetype="dashed") +
        geom_ribbon(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96, fill=threat), color=NA, alpha=.1) +
        #geom_errorbar(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96)) +
        geom_point(size=3) + geom_line() + 
        scale_color_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_shape_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_fill_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + 
        ylab(expression(Delta ~ "Heart Rate (bpm)")) + xlab("Trial Time (sec)") + labs(color="Threat", shape="Threat", fill="Threat") + myGgTheme)
#ggsave("plots/Heart Time.png", plot=heart.trialtime.plot, scale=1, device="png", dpi=300, units="in", width=1920/300, height = 1080/300)

# Split in HSA/LSA 
heart.ga.gen.time.spai = heart %>% filter(phase == "Gen") %>% 
  mutate(spai.split = ifelse(SPAI.z >= 0, "HSA", "LSA")) %>%
  group_by(threat, time, spai.split) %>% 
  summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T)) %>% 
  mutate(time = time %>% as.character() %>% as.numeric()) %>% 
  bind_rows(data.frame(threat=unique(heart$threat), time=0, spai.split = "HSA", HRchange=0, HRchange.se=0)) %>%
  bind_rows(data.frame(threat=unique(heart$threat), time=0, spai.split = "LSA", HRchange=0, HRchange.se=0)) %>% #add origin
  mutate(spai.split = spai.split %>% factor(levels=c("LSA", "HSA"))) #low anxiety first instead of alphabetical order
print(heart.trialtime.spai.plot <- heart.ga.gen.time.spai %>% ggplot(aes(x=time, y=HRchange, color=threat, group=threat, shape=threat)) + 
        geom_hline(yintercept = 0, linetype="dashed") +
        geom_ribbon(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96, fill=threat), color=NA, alpha=.1) +
        #geom_errorbar(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96)) +
        geom_point(size=3) + geom_line() + 
        facet_wrap(vars(spai.split), labeller = labeller(spai.split = c(LSA = "Low Social Anxiety", HSA = "High Social Anxiety"))) +
        scale_color_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_shape_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_fill_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + 
        ylab(expression(Delta ~ "Heart Rate (bpm)")) + xlab("Trial Time (sec)") + labs(color="Threat", shape="Threat", fill="Threat") + myGgTheme)
#ggsave("plots/Heart Time SPAI.png", plot=heart.trialtime.spai.plot, scale=1, device="png", dpi=300, units="in", width=1920*1.7/300, height = 1080/300)


#this is an oversimplification now (cf. exploration)
heart.simple = heart %>% #filter(shockPrior==F) %>% 
  group_by(subject, threat, pair, diagnostic, sex, phase, pairs, SPAI, SPAI.z, STAI, STAI.z) %>% #group_by_at(vars(-trial, -time, -shock, -shockPrior, -hrbl, -HRchange)) %>% 
  summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T), hrbl = mean(hrbl, na.rm=T))

heart.simple.gen = heart.simple %>% filter(phase == "Gen")

#plot across all trials (i.e. experiment time)
# heart.ga.trials = heart %>% group_by(trial, threat, subject) %>% summarise(HRchange.se=se(HRchange, na.rm=T), HRchange=mean(HRchange, na.rm=T)) %>% 
#   group_by(trial, threat) %>% summarise(HRchange.se=se(HRchange, na.rm=T), HRchange=mean(HRchange, na.rm=T))
# heart.ga.trials %>% ggplot(aes(x=trial, y=HRchange, color=threat, shape=threat, group=threat)) +
#   #geom_errorbar(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96)) + 
#   #geom_ribbon(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96, fill=threat), color=NA, alpha=.1) +
#   geom_point() + geom_line() +
#   geom_vline(xintercept = c(preAcqEnd, acqEnd), linetype=2) +
#   scale_color_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_shape_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_fill_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + 
#   #scale_y_reverse() +
#   guides(colour=guide_legend(reverse=T), shape=guide_legend(reverse=T), fill=guide_legend(reverse=T)) + 
#   ylab("Heart Rate Change (bpm)") + xlab("Trial Count") + labs(color="Threat", shape="Threat", fill="Threat")

#hr change generalization
heart.ga.gen.subj = heart.simple.gen %>% group_by(threat, subject) %>% summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T))
heart.ga.gen = heart.ga.gen.subj %>% summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T))
print(heart.gradient.plot <- heart.ga.gen %>% ggplot(aes(x=threat, y=HRchange, color=threat, group=NA)) + 
        #geom_dotplot(data=heart.ga.gen.subj, mapping=aes(group=threat, fill=threat), binaxis="y", alpha=.25, color="black", stackratio=1, stackdir="centerwhole", dotsize=.5) +
        geom_path(data=heart.ga.gen %>% filter(threat %in% c(1, 6)), aes(group=NA), color = "black", size=1.5) + #generalization line
        geom_errorbar(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96), size=1.5) +
        geom_line(size=1) + 
        geom_point(size=4.5) +
        scale_color_manual(values=colors)+#, guide=guide_legend(reverse=T)) + scale_y_reverse() +
        #scale_fill_manual(values=colors, guide=guide_legend(reverse=T)) +
        scale_fill_manual(values=rep("grey", 6), guide=guide_legend(reverse=T)) +
        scale_x_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) +
        ylab(expression(Delta ~ "Heart Rate (bpm)")) + xlab("Threat") + labs(color="Threat") +
        myGgTheme + theme(
          #aspect.ratio = 1,
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, "pt")),
          legend.position = "none"))
#ggsave("plots/Heart Gradient.png", plot=heart.gradient.plot, scale=1, device="png", dpi=300, units="in", width=1920/300, height = 1080/300)

heart.ga.gen %>% select(threat, HRchange, HRchange.se) #descriptive values

#Figure Heart
#cowplot::plot_grid(heart.trialtime.plot, heart.gradient.plot, ncol=1, labels="auto") %>% ggsave("figures/Figure Heart (old).png", plot=., scale=.7, device="png", dpi=300, units="in", width=6.5, height = 6.5 * 2 / sqrt(2))
{(heart.gradient.plot + (heart.trialtime.plot + theme(axis.title.y=element_blank())))/heart.trialtime.spai.plot + plot_annotation(tag_levels = 'a') + 
    plot_layout(axis_titles = "collect_y", guides = "collect")} %>% 
  ggsave("figures/Figure Heart.png", plot=., scale=1.45, device="png", dpi=300, units="in", width=6.5, height = 6.5 / sqrt(2))

#mean scores generalization phase by pairs
# heart.simple.gen %>% group_by(threat, pairs, subject) %>% summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T)) %>% 
#   group_by(threat, pairs) %>% summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T)) %>% 
#   ggplot(aes(x=threat, y=HRchange, color=pairs, group=pairs)) +
#   #geom_line(position=dodge) + 
#   geom_point(position=dodge) + 
#   geom_errorbar(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96), position=dodge, width=dodge.width) +
#   scale_y_reverse() + myGgTheme

#mean scores generalization phase by diagnostic and sex
# heart.simple.gen %>% group_by(threat, pair, subject) %>% summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T)) %>%
#   group_by(threat, pair) %>% summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T)) %>% 
#   ggplot(aes(x=threat, y=HRchange, color=pair, group=pair)) +
#   #geom_line(position=dodge) + 
#   geom_point(position=dodge) + 
#   geom_errorbar(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96), position=dodge, width=dodge.width) +
#   scale_y_reverse() + myGgTheme

# Inference tests ---------------------------------------------------------
ez::ezANOVA(data=heart.simple.gen, 
            dv=.(HRchange), wid=.(subject), 
            within=.(threat, diagnostic), 
            between=.(SPAI.z), observed=.(SPAI.z),
            #between=.(pairs),
            detailed=T, type=2) %>% apa::anova_apa(force_sph_corr=T)

ez::ezANOVA(data=heart.simple.gen, 
            dv=.(HRchange), wid=.(subject), 
            within=.(threat, diagnostic), 
            between=.(STAI.z), observed=.(STAI.z),
            #between=.(pairs),
            detailed=T, type=2) %>% apa::anova_apa(force_sph_corr=T)

# intercept
heart.simple.gen %>% group_by(subject) %>% summarise(HRchange = mean(HRchange)) %>% summarise(m = mean(HRchange), sd = sd(HRchange), se = se(HRchange))

# threat: t-tests of threat-levels vs. CS-
heart.simple.threat.gen = heart.simple %>% filter(phase == "Gen") %>% group_by(subject, threat) %>% 
  summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T))
for (i in 2:6) {
  levels = c(1, i)
  cat(paste0("\n\nComparing levels: ", paste(levels, collapse=" vs. "), "\n"))
  heart.simple.threat.gen %>% filter(threat %in% levels) %>% 
    t.test(HRchange ~ threat, ., #alternative="greater", 
           paired=T) %>% apa::t_apa(es_ci=T)
}

# SPAI
heart.simple.gen.lvl = heart.simple.gen %>% group_by(subject, SPAI, STAI) %>% summarise(HRchange.se = se(HRchange, na.rm=T), HRchange.m = mean(HRchange, na.rm=T))
heart.simple.gen.lvl %>% with(cor.test(HRchange.m, SPAI, alternative="less")) %>% correlation_out()
heart.simple.gen.lvl %>% with(cor.test(HRchange.m, SPAI)) %>% correlation_out()
print(heart.spai.plot <- heart.simple.gen.lvl %>% ggplot(aes(x=SPAI, y=HRchange.m, color=SPAI, fill=SPAI)) +
        geom_errorbar(aes(ymin=HRchange.m-HRchange.se*1.96, ymax=HRchange.m+HRchange.se*1.96), width=spai.width) +
        stat_smooth(method="lm", color = "black") +
        #geom_point(size=4, shape=21, color="black") +
        geom_point(size=4) + 
        ylab("Heart Rate Change (bpm)") +
        scale_color_viridis_c() + scale_fill_viridis_c() + myGgTheme + theme(legend.position = "none"))
#ggsave("plots/Heart SPAI.png", plot=heart.spai.plot, scale=1, device="png", dpi=300, units="in", width=1920/300, height = 1080/300)

#SPAI x threat x diagnostic
heart.simple.spaiXthreatXdiagnostic = heart.simple.gen %>% group_by(subject, SPAI, STAI, threat, diagnostic) %>% summarise(HRchange.se = se(HRchange, na.rm=T), HRchange.m = mean(HRchange, na.rm=T))
print(heart.spaiInteraction.plot <- heart.simple.spaiXthreatXdiagnostic %>% ggplot(aes(x=SPAI, y=HRchange.m, color=SPAI, fill=SPAI)) +
        #facet_grid(rows=vars(threat), cols=vars(diagnostic),  labeller = label_both) +
        facet_grid(rows=vars(diagnostic), cols=vars(threat),  labeller = label_both) +
        geom_errorbar(aes(ymin=HRchange.m-HRchange.se*1.96, ymax=HRchange.m+HRchange.se*1.96), width=spai.width) +
        stat_smooth(method="lm", color = "black") +
        #geom_point(size=4, shape=21, color="black") +
        geom_point(size=4) + 
        ylab("Heart Rate Change (bpm)") +
        scale_color_viridis_c() + scale_fill_viridis_c() + myGgTheme + theme(legend.position = "none"))
heart.simple.spaiXthreatXdiagnostic.r = heart.simple.spaiXthreatXdiagnostic %>% group_by(threat, diagnostic) %>% summarise(rtest = cor.test(HRchange.m, SPAI) %>% apa::cor_apa(r_ci=T, print=F), r = cor(HRchange.m, SPAI)) %>% arrange(desc(r)) %>% 
  mutate(helper = rtest %>% gsub("\\[", "; ", .) %>% gsub("\\]", "; ", .)) %>% separate_wider_delim(helper, delim="; ", names=c(NA, "CI95.lo", "CI95.hi", NA)) %>% mutate(across(contains("CI"), as.numeric))
heart.simple.spaiXthreatXdiagnostic.r
heart.simple.spaiXthreatXdiagnostic.r %>% ggplot(aes(x = threat, y = r, color = diagnostic)) + 
  facet_wrap(vars(diagnostic)) +
  geom_hline(yintercept = 0) + geom_errorbar(aes(ymin = CI95.lo, ymax = CI95.hi)) + geom_line() + geom_point() + myGgTheme


#hr change acquisition
heart.ga.acq = heart %>% filter(phase == "Acq") %>% subset((subject %in% exclusions.onlyGen & phase != "Gen") == F) %>% #exclude habituation & acquisition for onlyGen subjects (see "0 General.R")
  group_by(threat, subject) %>% summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T)) %>%
  group_by(threat) %>% summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T))
heart.ga.acq %>% ggplot(aes(x=threat, y=HRchange, color=threat, group=NA)) + 
  geom_point() + #geom_line() + 
  geom_errorbar(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96)) +
  scale_color_manual(values=colors[c(1, 6)]) + #scale_y_reverse() +
  myGgTheme
heart.simple %>% filter(phase == "Acq") %>% subset((subject %in% exclusions.onlyGen & phase != "Gen") == F) %>% #exclude habituation & acquisition for onlyGen subjects (see "0 General.R")
  group_by(subject, threat) %>% summarise(HRchange = mean(HRchange, na.rm=T)) %>% 
  t.test(HRchange ~ threat, ., alternative="greater", 
         paired=T) %>% apa::t_apa(es_ci=T) #schoRsch::t_out()



# Exploration -------------------------------------------------------------
heart.gen = heart %>% filter(phase == "Gen") #%>% merge(questionnaires, by="subject")

# heart.anova.big = ez::ezANOVA(data=heart.gen %>% mutate(time = time %>% as.factor()), #treat time as factor to allow for non-linear dynamics
#                               dv=.(HRchange), wid=.(subject), 
#                               within=.(threat, diagnostic, time), 
#                               #within_full = .(trial),
#                               #between=.(pairs),
#                               detailed=T, type=2)
# heart.anova.big %>% apa::anova_apa(force_sph_corr=T)

heart.anova.big.spai = ez::ezANOVA(data=heart.gen %>% mutate(time = time %>% as.factor()), #treat time as factor to allow for non-linear dynamics
                                   dv=.(HRchange), wid=.(subject), 
                                   within=.(threat, diagnostic, time), 
                                   between=.(SPAI.z), observed=SPAI.z,
                                   #between=.(pairs),
                                   detailed=T, type=2)
heart.anova.big.spai %>% apa::anova_apa(force_sph_corr=T)

#SPAI main effect
heart.spai.plot
with(heart.simple.gen.lvl, cor.test(SPAI, HRchange.m)) %>% apa::cor_apa(r_ci=T)

#SPAI x time
requirePackage("lme4")
requirePackage("lmerTest")
hr.spai.lmm = lmer(HRchange ~ time*SPAI + (1|subject), #emulate effects of interest in LMM since function "predict" can't deal with aov.list
                   data=heart.gen %>% mutate(time = time %>% as.factor(), SPAI = scale(SPAI)))
#hr.spai.lmm %>% summary()
SDs = -2:2 #SDs = c(-1, 1)
spai.descr = with(questionnaires, data.frame(SPAI.raw=mean(SPAI)+sd(SPAI)*SDs, SPAI=SDs))
heart.spaiSD.predict = spai.descr %>%
  expand.grid.df(data.frame(time=heart$time %>% unique() %>% as.factor())) %>% #all combinations with time
  expand.grid.df(data.frame(subject=0)) %>% #include mockup subject for population-level prediction (could avoid this and add "re.form=NA" to predict but less flexibility with additional subject-level predictions or based on non-lmer model objects, i.e. predict.merMod)
  mutate(HRchange = predict(hr.spai.lmm, newdata=., allow.new.levels=T), #predict values from time & z-scaled SPAI
         time = time %>% as.character() %>% as.numeric()) %>% #time as numeric to allow adding time origin (new factor level -> error)
  rename(SPAI.z = SPAI, SPAI = SPAI.raw) %>% #rename to use raw scores for plotting
  bind_rows(., data.frame(SPAI=unique(.$SPAI), time=0, HRchange=0, subject=0)) #add origin 

print(heart.spaiXtime.plot <- heart.gen %>%
        group_by(time, SPAI, subject) %>% summarise(HRchange.se = se(HRchange), HRchange = mean(HRchange)) %>% 
        mutate(time = time %>% as.character() %>% as.numeric()) %>% 
        bind_rows(., data.frame(subject=unique(.$subject), SPAI=head(.$SPAI, n=length(unique(.$subject))), time=0, HRchange=0, HRchange.se=0)) %>% #add origin
        ggplot(aes(x=time, y=HRchange, group=subject, color=SPAI)) +
        #stat_smooth(method="lm") +
        geom_line(size=1, alpha=.5) +
        #geom_line(aes(y = predict(hr.spai.lmm), group=SPAI), size=4) +
        geom_line(data=heart.spaiSD.predict, mapping=aes(group=SPAI), size=4) +
        ylab("Heart Rate Change (bpm)") + xlab("Trial Time (sec)") +
        scale_color_viridis_c() + myGgTheme)

#threat 
heart.gradient.plot

#threat x time
heart.trialtime.plot

#SPAI x threat x diagnostic
hr.spai.3x.lmm = data=heart.gen %>% lmer(HRchange ~ SPAI*threat*diagnostic + (1|subject), data=.) #emulate effects of interest in LMM since function "predict" can't deal with aov.list
#hr.spai.3x.lmm %>% summary()
SDs = -2:2
#spai.descr = with(questionnaires, data.frame(SPAI.raw=mean(SPAI)+sd(SPAI)*SDs, SPAI=SDs))
heart.spaiSD.3x.predict = spai.descr %>%
  expand.grid.df(data.frame(threat=heart$threat %>% unique() %>% as.factor())) %>% #all combinations with threat
  expand.grid.df(data.frame(diagnostic=heart$diagnostic %>% unique() %>% as.factor())) %>% #all combinations with diagnostic
  expand.grid.df(data.frame(subject=0)) %>% #include mockup subject for population-level prediction (could avoid this and add "re.form=NA" to predict but less flexibility with additional subject-level predictions or based on non-lmer model objects, i.e. predict.merMod)
  mutate(HRchange = predict(hr.spai.3x.lmm, newdata=., allow.new.levels=T)) %>% #predict values from time & z-scaled SPAI
  rename(SPAI.z = SPAI, SPAI = SPAI.raw) #rename to use raw scores for plotting

print(heart.spai.all.plot <- heart.spaiSD.3x.predict %>%
        group_by(SPAI, threat, diagnostic, subject) %>% summarise(HRchange.se = se(HRchange), HRchange = mean(HRchange)) %>% 
        mutate(threat = threat %>% recode_factor(`1` = "CS-", `2` = "GS1", `3` = "GS2", `4` = "GS3", `5` = "GS4", `6` = "CS+")) %>% 
        ggplot(aes(x=threat, y=HRchange, group=SPAI, color=SPAI)) +
        geom_line(size=4) +
        facet_wrap(vars(diagnostic)) +
        ylab("Heart Rate Change (bpm)") + #xlab("Trial Time (sec)") +
        scale_color_viridis_c() + myGgTheme)

#SPAI x threat x diagnostic x time
hr.spai.all.lmm = lmer(HRchange ~ SPAI*threat*diagnostic*time + (1|subject), #emulate effects of interest in LMM since function "predict" can't deal with aov.list
                   data=heart.gen %>% mutate(time = time %>% as.factor(), SPAI = scale(SPAI)))
#hr.spai.all.lmm %>% summary()
SDs = -2:2
#spai.descr = with(questionnaires, data.frame(SPAI.raw=mean(SPAI)+sd(SPAI)*SDs, SPAI=SDs))
heart.spaiSD.all.predict = spai.descr %>%
  expand.grid.df(data.frame(threat=heart$threat %>% unique() %>% as.factor())) %>% #all combinations with threat
  expand.grid.df(data.frame(diagnostic=heart$diagnostic %>% unique() %>% as.factor())) %>% #all combinations with diagnostic
  expand.grid.df(data.frame(time=heart$time %>% unique() %>% as.factor())) %>% #all combinations with time
  expand.grid.df(data.frame(subject=0)) %>% #include mockup subject for population-level prediction (could avoid this and add "re.form=NA" to predict but less flexibility with additional subject-level predictions or based on non-lmer model objects, i.e. predict.merMod)
  mutate(HRchange = predict(hr.spai.all.lmm, newdata=., allow.new.levels=T), #predict values from time & z-scaled SPAI
         time = time %>% as.character() %>% as.numeric()) %>% #time as numeric to allow adding time origin (new factor level -> error)
  rename(SPAI.z = SPAI, SPAI = SPAI.raw) #rename to use raw scores for plotting
heart.spaiSD.all.predict = heart.spaiSD.all.predict %>% 
  bind_rows(heart.spaiSD.all.predict %>% filter(time==.5) %>% mutate(time=0, subject=0, HRchange=0)) #add origin 

print(heart.spai.all.plot <- heart.spaiSD.all.predict %>%
        group_by(SPAI, threat, diagnostic, time, subject) %>% summarise(HRchange.se = se(HRchange), HRchange = mean(HRchange)) %>% 
        mutate(time = time %>% as.character() %>% as.numeric(),
               threat = threat %>% recode_factor(`1` = "CS-", `2` = "GS1", `3` = "GS2", `4` = "GS3", `5` = "GS4", `6` = "CS+")) %>% 
        ggplot(aes(x=time, y=HRchange, group=SPAI, color=SPAI)) +
        geom_line(size=4) +
        facet_grid(rows=vars(threat), cols=vars(diagnostic), labeller = label_both) +
        ylab("Heart Rate Change (bpm)") + xlab("Trial Time (sec)") +
        scale_color_viridis_c() + myGgTheme)

# Gradient Analysis -----------------------------------------------------------
individualPlots = F

#subject-level means generalization phase
hr.gradients.simple = data.frame(subject=character(), lds=numeric(), diff=numeric(), level=numeric())
for (s in heart.simple.threat.gen$subject %>% unique()) {
  hr.temp = heart.simple.threat.gen %>% filter(subject==s)
  
  gradient = gradient.analysis(hr.temp$HRchange)
  hr.gradients.simple = rbind(hr.gradients.simple, c(s, gradient))
  
  if (individualPlots) {
    gradient = gradient %>% signif(3)
    print(hr.temp %>% ggplot(aes(x=threat, y=HRchange, color=threat)) + #generalization line
            geom_path(data=hr.temp %>% filter(threat %in% c(1, 6)), aes(group=NA), color = "black", size=1.5) +
            geom_point() + scale_color_manual(values=colors) + 
            geom_errorbar(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96)) +
            ggtitle(paste0(s, ": lds = ", gradient[1], ", diff = ", gradient[2], ", level = ", gradient[3])))
  }
}

names(hr.gradients.simple) = c("subject", "HR_Gen_all_lds", "HR_Gen_all_diff", "HR_Gen_all_level")

#subject-level means generalization phase by diagnostic
hr.gradients.diagnostic = data.frame()
for (s in heart.simple.gen$subject %>% unique()) {
  hr.temp = heart.simple.gen %>% filter(subject==s)
  
  gradient.eyes = "eyes" %>% {filter(hr.temp, diagnostic==.)} %>% .$HRchange %>% gradient.analysis()
  gradient.m_n = "mouth/nose" %>% {filter(hr.temp, diagnostic==.)} %>% .$HRchange %>% gradient.analysis()
  hr.gradients.diagnostic = rbind(hr.gradients.diagnostic, c(s, gradient.eyes, gradient.m_n))
  
  if (individualPlots) {
    print(hr.temp %>% ggplot(aes(x=threat, y=HRchange, color=diagnostic, shape=diagnostic)) + 
            #geom_path(data=hr.temp %>% filter(threat %in% c(1, 6)), aes(group=NA), color = "black", size=1.5) + #generalization line
            geom_point(position=dodge) + 
            geom_errorbar(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96), position=dodge, width=dodge.width) + ggtitle(s))
  }
}

names(hr.gradients.diagnostic) = c("lds", "diff", "level") %>% {c("subject", paste0("HR_Gen_eyes_", .), paste0("HR_Gen_mn_", .))}
hr.gradients = merge(hr.gradients.simple, hr.gradients.diagnostic, by="subject")

heart.wide = hr.gradients %>% tibble() #TODO add more info, e.g. subject-average of HR baseline?
#all(heart.wide == read_rds("heart.wide.rds" %>% paste0(path.rds, .)), na.rm=T) #check equivalence of processing
#heart.wide %>% write_rds("heart.wide.rds" %>% paste0(path.rds, .))
