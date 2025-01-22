if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
#source("0 General.R" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/C10 - Gamer Fear/3 Diagnostic Generalization/3 Analysis/main Anx/", .))

outlierLimit.ratings = .25 #maximum percentage of invalid trials per subject
createFmriFiles = F

# Functions ---------------------------------------------------------------
{ toPulseTime = function(x, pulses) { #x: vector of times to convert, pulses: vector of pulse times (in the same unit as x)
    result = c()
    for (xi in x) {
      indices = which.max(xi < pulses) %>% c(.-1, .) #indices of pulses before and after xi
      if (indices %>% is.na() %>% any() || any(indices < 1)) result = result %>% c(NA)
      else {
        result = result %>% c(
          min(indices) + (xi - pulses[min(indices)]) / diff(pulses[indices])
        )
      }
    }
    return(result)
  }
  
  toTrialNum = function(x, trialOnsets) { #x: vector of times to convert, trialOnsets: vector of trial onset (in the same unit as x)
    result = c()
    for (xi in x) result = result %>% c(sum(xi >= trialOnsets))
    return(result)
  }
}


# Read & tidy data --------------------------------------------------------
ratings.list = vector("list", length(files.rating)); conditions.list = vector("list", length(files.rating))
for (filename in files.rating) {
  code = pathToCode(filename, file.ext="-")
  #file = read_table2(filename, skip=3) %>% #read.delim(filename, skip=3) %>% 
  file = read_delim(filename, delim="\t", skip=3, show_col_types=F, name_repair="minimal") %>% #name_repair="minimal" to suppress warning of name_repair but error when selecting
    select(Subject, Time, `Event Type`, Code) %>% #selection before rename to avoid errors
    rename(Event.Type = `Event Type`) %>% 
    filter(Event.Type!="Quit") #%>% mutate(Time = Time + Uncertainty, TTime = TTime + Uncertainty) #correct pulse times
  
  # pulses = file %>% filter(Event.Type=="Pulse")
  # diffs = pulses %>% .$Time %>% diff()
  # #diffs %>% {round(. / 10)} %>% hist()
  # 
  # pulses = file %>% filter(Event.Type=="Pulse") %>% mutate(pulseN = 1:n())
  # pulses %>% nrow() %>% paste0(code, ": ", ., " pulses") %>% print()
  # jitter = pulses %>% .$Time %>% diff() %>% {round(. / 10)}
  # outliers = jitter > mean(jitter) + sd(jitter) | jitter < mean(jitter) - sd(jitter)
  # # pulses %>% filter(c(F, outliers)) %>% select(Subject, pulseN, Uncertainty) %>% 
  # #   mutate(TR = jitter[outliers], jitter = TR - round(mean(jitter)), Uncertainty = Uncertainty %>% {round(. / 10)}) %>% print()
  # #jitter %>% hist(main=code)
  
  fileLong = file %>% mutate(
    #PTime = toPulseTime(Time, pulses$Time),
    Event.Type = as.character(Event.Type),
    Code = as.character(Code),
    Code = case_when(Event.Type=="Response" ~ as.character(as.integer(Code) - 1), 
                     #Code %>% grepl(stimExt, .) ~ paste0(substr(Code, 1, 1), "_", Code %>% substr(3, 5) %>% as.integer() %>% {. / 20 + 1}),
                     TRUE ~ Code)) %>% 
    filter(Event.Type %in% c("Picture", "Response"), #!is.na(PTime), 
           !(Code %in% c("fixation", "rating", "noshock", "shock", "shock2", "shock3"))) %>% #noshock/shock: flag telling the kind of trial (timing not indicative of shock onset)
    select(Subject, Time, Event.Type, Code)
  
  file = fileLong %>% select(Subject, Event.Type, Code, Time) %>% mutate(Event.Type = ifelse(Code=="shock1", "shock", Event.Type))
  trialOnsets = file %>% filter(Event.Type=="Picture") %>% .$Time
  file = file %>% mutate(Trial = toTrialNum(Time, trialOnsets),
                         helper = paste(Event.Type, Trial)) %>% 
    filter(!duplicated(helper)) %>% select(Subject, Trial, Event.Type, Code) %>% 
    spread(Event.Type, Code) %>% separate(Picture, c("pair", "morph"), sep="_", remove=F) %>% 
    transmute(subject = Subject, trial = Trial, pic = Picture,
              morph = morph %>% gsub(stimExt, "", ., fixed=T) %>% as.integer(), 
              pair = as.integer(pair), 
              rating = as.integer(Response), 
              shock = ifelse(is.na(shock), FALSE, TRUE)) %>% filter(pic %>% is.na() %>% {. == FALSE})
  
  conditions = file %>% filter(shock) %>% head(2) %>% mutate(csp = ifelse(morph==100, 1, 2)) %>% select(subject, pair, csp)
  toSwitch = conditions %>% filter(csp != 1) %>% .$pair
  file = file %>% mutate(threat = morph %>% {. / 20 + 1}, #pretend 100 would always be cs+ (true if csp==1)
                         threat = ifelse(pair %in% toSwitch, threatLevels.n+1 - threat, threat)) #check if threat level has to be inverted (because csp==2)
  
  ratings.list[[code]] = file
  conditions.list[[code]] = conditions
  
  fileLong[grepl(stimExt, fileLong$Code, fixed=T), "Code"] = paste(ifelse(file$pair %in% c(1, 3), "E", "M"), file$threat, sep="_") #replace picture name by threat level from "file" dataframe
  if (createFmriFiles) fileLong %>% mutate(Event = paste(Event.Type, Code, sep="_"),
                                           Event = ifelse(Code=="shock1", Code, Event)) %>% 
    select(Subject, PTime, Event) %>% write.table(paste0(filename %>% gsub("/log/", "/MRI/", .), ".csv"), sep="\t", row.names=FALSE)
}

#merge into dataframe & tidy
ratings.all = ratings.list %>% bind_rows() %>% separate(subject, c("subject", "block"), sep="_") %>% 
  mutate(subject = subject %>% gsub("vp", "", ., fixed=T) %>% as.integer(),
         block = as.integer(block),
         trial = case_when(block==2 ~ trial + as.integer(acqEnd),
                           block==3 ~ trial + as.integer(gen1End),
                           TRUE ~ trial), 
         threat_num = threat,
         threat = case_when(threat_num==1 ~ "CS-",
                            #threat==2 ~ "GS1", threat==3 ~ "GS2", threat==4 ~ "GS3", threat==5 ~ "GS4",
                            threat_num==threatLevels.n ~ "CS+",
                            TRUE ~ paste0("GS", threat_num-1)) %>% factor(levels = c("CS-", paste0("GS", 1:(threatLevels.n-2)), "CS+")),
         threat_both = paste0(threat_num, ": ", threat) %>% as.factor(),
         phase = case_when(trial < preAcqEnd ~ "Hab",
                           trial <    acqEnd ~ "Acq",
                           TRUE ~              "Gen") %>% factor(levels=c("Hab", "Acq", "Gen")),
         diagnostic = as.factor(ifelse(pair %% 2 == 0, "Mouth/Nose", "Eyes")),
         sex = as.factor(ifelse(pair > 2, "female", "male")),
         pairs = ifelse(pair %in% 2:3, 2, 1),
         pair = as.factor(pair),
         condition = 10*as.integer(pair) + threat_num)
#TODO subject 54: ratings of first block missing => conditions NA

conditions = conditions.list %>% bind_rows() %>% separate(subject, c("subject", "block"), sep="_") %>% arrange(subject, pair) %>% 
  mutate(subject = subject %>% gsub("vp", "", ., fixed=T) %>% as.integer(), pair = as.integer(pair))
#conditions %>% group_by(subject, pair) %>% summarise(csp.sd = sd(csp)) %>% filter(csp.sd != 0) #show inconsistencies
conditions = conditions %>% group_by(subject, pair) %>% summarise(csp = mean(csp))

conditions.csp = conditions %>% mutate(pairs = ifelse(pair %in% 2:3, 2, 1),
  pair = ifelse(pair < mean(unique(conditions$pair)), "csp1", "csp2")) %>% spread(pair, csp) %>% ungroup()

rm(ratings.list, conditions.list)
#all(ratings.all == read_rds("ratings.rds" %>% paste0(path.rds, .)), na.rm=T) #check equivalence of processing
#write_rds(ratings.all, "ratings.rds" %>% paste0(path.rds, .)); write_rds(conditions, "conditions.rds" %>% paste0(path.rds, .)); write_rds(conditions.csp, "conditions.csp.rds" %>% paste0(path.rds, .))

# Exclusions --------------------------------------------------------------
#ratings.all = read_rds("ratings.rds" %>% paste0(path.rds, .)); conditions = read_rds("conditions.rds" %>% paste0(path.rds, .)); conditions.csp = read_rds("conditions.csp.rds" %>% paste0(path.rds, .))

ratings = ratings.all %>% subset(subject %in% exclusions == F) %>% #TODO rather set all ratings to NA? better for preservation of conditions. but then Eye analysis has to be changed, cp. filter(condition != 0) #no valid fixation during trial
  subset((subject %in% exclusions.onlyGen & phase != "Gen") == F) #exclude habituation & acquisition for onlyGen subjects (see "0 General.R")
#conditions.csp %>% filter(subject %in% exclusions == F) %>% .[-1] %>% apply(2, table) # every condition independently
#conditions.csp %>% filter((subject %in% c(exclusions, exclusions.onlyGen)) == F) %>% #only subjects that have all phases
conditions.csp %>% filter((subject %in% exclusions) == F) %>% #all subjects with valid generalization phase
  transmute(condition = paste0(pairs, csp1, csp2)) %>% table() #condition combinations

ratings.valid = ratings %>% group_by(phase, subject) %>% summarise(NAs = rating %>% is.na() %>% sum() / n()) %>% arrange(desc(NAs))
ratings %>% group_by(subject) %>% summarise(NAs = rating %>% is.na() %>% sum() / n()) %>% arrange(desc(NAs)) %>% summarize(M = mean(NAs), SD = sd(NAs))
ratings %>% filter(phase == "Gen") %>% group_by(subject) %>% summarise(NAs = rating %>% is.na() %>% sum() / n()) %>% arrange(desc(NAs)) %>% summarize(M = mean(NAs), SD = sd(NAs))

#hist(ratings.valid$NAs); abline(v=outlierLimit.ratings, col="red", lwd=3, lty=2) #, breaks=seq(0, outlierLimit.ratings, length.out=20+1))
ratings.valid %>%
  ggplot(aes(x=NAs, fill=phase)) + geom_histogram(boundary=outlierLimit.ratings, color="black") + 
  geom_vline(xintercept=outlierLimit.ratings, color="red") + myGgTheme + scale_y_continuous(breaks=scales::breaks_pretty())
#ratings.valid %>% arrange(desc(NAs))
print(problem <- ratings %>% group_by(subject, phase) %>% summarise(NAs = rating %>% is.na() %>% sum() / n()) %>% arrange(desc(NAs)) %>% filter(NAs > outlierLimit.ratings))

exclusions.rating = problem %>% filter(phase=="Gen") %>% .$subject %>% unique() %>% c(exclusions)
ratings = ratings %>% filter(subject %in% exclusions.rating == F)


# Grand Average Plots ------------------------------------------------------
ratings.gen = ratings %>% filter(phase=="Gen") %>% left_join(questionnaires, by="subject") %>% 
  mutate(SPAI.z = scale(SPAI)[,1], STAI.z = scale(STAI)[,1])

ratings.first.gen = ratings %>% filter(block==2) %>% left_join(questionnaires, by="subject") %>%
  mutate(SPAI.z = scale(SPAI)[,1], STAI.z = scale(STAI)[,1])

#plot across all trials
ratings.ga.trials = ratings %>% group_by(trial, threat, threat_num, threat_both, subject) %>% summarise(rating.se=se(rating, na.rm=T), rating=mean(rating, na.rm=T)) %>% 
  summarise(rating.se=se(rating, na.rm=T), rating=mean(rating, na.rm=T))
print(ratings.trials.plot <- ratings.ga.trials %>% ggplot(aes(x=trial, y=rating, color=threat, shape=threat, group=threat)) +
        geom_ribbon(aes(ymin=rating-rating.se*1.96, ymax=rating+rating.se*1.96, fill=threat), color=NA, alpha=.1) +
        geom_point() + geom_line() +
        geom_vline(xintercept = c(preAcqEnd, acqEnd), linetype=2) +
        
        #annotations near the lines
        # annotate("label", label="Habituation", x = preAcqEnd-.5, hjust = "right", y = Inf, vjust = "top") +
        # annotate("label", label="Acquisition", x = mean(c(preAcqEnd, acqEnd)), y = Inf, vjust = "top") +
        # annotate("label", label="Generalization", x = acqEnd+.5, hjust="left", y = Inf, vjust = "top") +
        
        #annotations with equal distance (good of same amount of letters per label)
        # annotate("label", label="Hab", x = mean(c(0, preAcqEnd)), y = Inf, vjust = "top") +
        # annotate("label", label="Acq", x = mean(c(preAcqEnd, acqEnd)), y = Inf, vjust = "top") +
        # annotate("label", label="Gen", x = 2*mean(c(preAcqEnd, acqEnd)) - mean(c(0, preAcqEnd)), y = Inf, vjust = "top") +
        
        #annotations in the middle of each phase
        annotate("label", label="Hab.", x = mean(c(0, preAcqEnd)), y = Inf, vjust = "top") +
        annotate("label", label="Acquisition", x = mean(c(preAcqEnd, acqEnd)), y = Inf, vjust = "top") +
        annotate("label", label="Generalization", x = mean(c(acqEnd, trials.n)), y = Inf, vjust = "top") +
        
        coord_cartesian(clip="off") +
        #scale_color_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_shape_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) + scale_fill_manual(values=colors, labels=c("CS-", paste0("GS", 1:4), "CS+")) +
        scale_color_manual(values=colors) + scale_fill_manual(values=colors) + 
        guides(colour=guide_legend(reverse=T), shape=guide_legend(reverse=T), fill=guide_legend(reverse=T)) + 
        ylab("Threat Rating (1-5)") + xlab("Trial Count") + labs(color="Threat", shape="Threat", fill="Threat") + 
        scale_x_continuous(expand=c(1/trials.n, 1/trials.n)) + 
        myGgTheme)
#ggsave("plots/Ratings Trials.png", plot=ratings.trials.plot, scale=1, device="png", dpi=300, units="in", width=1920/300, height = 1080/300)

#mean scores generalization phase
ratings.ga.gen.subj = ratings.gen %>% group_by(threat, threat_num, threat_both, subject) %>% summarise(rating.se=se(rating, na.rm=T), rating=mean(rating, na.rm=T))
ratings.ga.gen = ratings.ga.gen.subj %>% summarise(rating.se=se(rating, na.rm=T), rating=mean(rating, na.rm=T))
print(ratings.gradient.plot <- ratings.ga.gen %>% ggplot(aes(x=threat_both, y=rating, color=threat, group=NA)) + 
        geom_dotplot(data=ratings.ga.gen.subj, mapping=aes(group=threat, fill=threat), binaxis="y", alpha=.25, color="black", stackratio=1, stackdir="centerwhole", dotsize=.5) +
        #try ggbeeswarm ? or position_jitter ? or geom_quasirandom ? or violinplot?
        #geom_path(data=ratings.ga.gen %>% filter(threat_num %in% c(1, threatLevels.n)), color = "black", size=1.5) + #generalization line
        geom_line(linewidth=1) + geom_point(size=4.5) + 
        geom_errorbar(aes(ymin=rating-rating.se*1.96, ymax=rating+rating.se*1.96), size=1.5) +
        scale_color_manual(values=colors, guide=guide_legend(reverse=T)) +
        scale_fill_manual(values=colors, guide=guide_legend(reverse=T)) + #individual dots in color
        #scale_fill_manual(values=rep("grey", 6), guide=guide_legend(reverse=T)) + #individual dots grey
        scale_x_discrete(labels=levels(ratings.ga.gen$threat)) +
        ylab("Threat Rating (1-5)") + xlab("Threat") + labs(color="Threat", fill="Threat") +
        myGgTheme)
#ggsave("plots/Ratings Gradient.png", plot=ratings.gradient.plot, scale=1, device="png", dpi=300, units="in", width=1920/300, height = 1080/300)

ratings.ga.gen.SPAI.subj = ratings.gen %>% 
  mutate(SPAI = if_else(SPAI < median(SPAI, na.rm=T), "low", "high")) %>% 
  #mutate(SPAI = SPAI %>% ntile(3) %>% case_match(1 ~ "low", 2 ~ "middle", 3 ~ "high")) %>% 
  group_by(SPAI, threat, threat_num, threat_both, subject) %>% summarise(rating.se=se(rating, na.rm=T), rating=mean(rating, na.rm=T))
ratings.ga.gen.SPAI = ratings.ga.gen.SPAI.subj %>% summarise(rating.se=se(rating, na.rm=T), rating=mean(rating, na.rm=T), n=n())
print(ratings.gradient.SPAI <- ratings.ga.gen.SPAI %>% ggplot(aes(x=threat_both, y=rating, color=SPAI, group=SPAI)) +
        #geom_dotplot(data=ratings.ga.gen.SPAI.subj %>% filter(SPAI=="low"), mapping=aes(group=interaction(threat, SPAI), fill=SPAI), stackdir="up", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.5) +
        #geom_dotplot(data=ratings.ga.gen.SPAI.subj %>% filter(SPAI=="high"), mapping=aes(group=interaction(threat, SPAI), fill=SPAI), stackdir="down", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.5) +
        geom_line(linewidth=1, position=dodge) + geom_point(size=4.5, position=dodge) + 
        geom_errorbar(aes(ymin=rating-rating.se*1.96, ymax=rating+rating.se*1.96), size=1.5, position=dodge, width=dodge.width) +
        scale_color_viridis_d(direction=-1) + scale_fill_viridis_d(direction=-1) +
        scale_x_discrete(labels=levels(ratings.ga.gen$threat)) +
        ylab("Threat Rating (1-5)") + xlab("Threat") +
        myGgTheme)
#ggsave("plots/Ratings Gradient SPAI.png", plot=ratings.gradient.SPAI, scale=1, device="png", dpi=300, units="px", width=1920, height = 1080)

ratings.ga.gen.STAI.subj = ratings.gen %>% mutate(STAI = if_else(STAI < median(STAI, na.rm=T), "low", "high")) %>% 
  group_by(STAI, threat, threat_num, threat_both, subject) %>% summarise(rating.se=se(rating, na.rm=T), rating=mean(rating, na.rm=T))
ratings.ga.gen.STAI = ratings.ga.gen.STAI.subj %>% summarise(rating.se=se(rating, na.rm=T), rating=mean(rating, na.rm=T), n=n())
print(ratings.gradient.STAI <- ratings.ga.gen.STAI %>% ggplot(aes(x=threat_both, y=rating, color=STAI, group=STAI)) +
        #geom_dotplot(data=ratings.ga.gen.STAI.subj %>% filter(STAI=="low"), mapping=aes(group=interaction(threat, STAI), fill=STAI), stackdir="up", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.5) +
        #geom_dotplot(data=ratings.ga.gen.STAI.subj %>% filter(STAI=="high"), mapping=aes(group=interaction(threat, STAI), fill=STAI), stackdir="down", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.5) +
        geom_line(linewidth=1, position=dodge) + geom_point(size=4.5, position=dodge) + 
        geom_errorbar(aes(ymin=rating-rating.se*1.96, ymax=rating+rating.se*1.96), size=1.5, position=dodge, width=dodge.width) +
        scale_color_viridis_d(direction=-1) + scale_fill_viridis_d(direction=-1) +
        scale_x_discrete(labels=levels(ratings.ga.gen$threat)) +
        ylab("Threat Rating (1-5)") + xlab("Threat") +
        myGgTheme)
#ggsave("plots/Ratings Gradient STAI.png", plot=ratings.gradient.STAI, scale=1, device="png", dpi=300, units="px", width=1920, height = 1080)

#Figure Ratings
#cowplot::plot_grid(ratings.trials.plot, ratings.gradient.plot + theme(legend.position="none"), ncol=1, labels="auto") %>% ggsave("figures/Figure Ratings (old).png", plot=., scale=1, device="png", dpi=300, units="in", width=6.5, height = 6.5 * 2 / sqrt(2))
#{ratings.trials.plot / ratings.gradient.plot / ratings.gradient.STAI + plot_annotation(tag_levels = 'a')} %>% ggsave("figures/Figure Ratings.png", plot=., scale=1.8, device="png", dpi=300, units="in", width=6.5/2, height = 6.5/2 * 3 / sqrt(2))
{ratings.trials.plot / ((ratings.gradient.plot + theme(legend.position="none")) + (ratings.gradient.SPAI + ylab("")) + plot_layout(widths=c(1.25, 1))) + plot_annotation(tag_levels = 'a')} %>% 
  ggsave("figures/Figure Ratings (SPAI).png", plot=., scale=1.45, device="png", dpi=300, units="in", width=6.5, height = 6.5 / sqrt(2))
{ratings.trials.plot / ((ratings.gradient.plot + theme(legend.position="none")) + (ratings.gradient.STAI + ylab("")) + plot_layout(widths=c(1.25, 1))) + plot_annotation(tag_levels = 'a')} %>% 
  ggsave("figures/Figure Ratings.png", plot=., scale=1.45, device="png", dpi=300, units="in", width=6.5, height = 6.5 / sqrt(2))


# #mean scores first generalization phase
# ratings.ga.gen.subj = ratings.first.gen %>% group_by(threat, threat_num, threat_both, subject) %>% summarise(rating.se=se(rating, na.rm=T), rating=mean(rating, na.rm=T))
# ratings.ga.gen = ratings.ga.gen.subj %>% summarise(rating.se=se(rating, na.rm=T), rating=mean(rating, na.rm=T))
# print(ratings.gradient.plot <- ratings.ga.gen %>% ggplot(aes(x=threat_both, y=rating, color=threat, group=NA)) + 
#         geom_dotplot(data=ratings.ga.gen.subj, mapping=aes(group=threat, fill=threat), binaxis="y", alpha=.25, color="black", stackratio=1, stackdir="centerwhole", dotsize=.5) +
#         #try ggbeeswarm ? or position_jitter ? or geom_quasirandom ?
#         #geom_path(data=ratings.ga.gen %>% filter(threat_num %in% c(1, threatLevels.n)), color = "black", size=1.5) + #generalization line
#         geom_line(size=1) + geom_point(size=4.5) + 
#         geom_errorbar(aes(ymin=rating-rating.se*1.96, ymax=rating+rating.se*1.96), size=1.5) +
#         scale_color_manual(values=colors, guide=guide_legend(reverse=T)) +
#         #scale_fill_manual(values=colors, guide=guide_legend(reverse=T)) + #individual dots in color
#         scale_fill_manual(values=rep("grey", 6), guide=guide_legend(reverse=T)) + #individual dots grey
#         scale_x_discrete(labels=levels(ratings.ga.gen$threat)) +
#         ylab("Threat Rating (1-5)") + xlab("Threat") + labs(color="Threat", fill="Threat") +
#         myGgTheme + theme(legend.position="none"))

#mean scores generalization phase by diagnostic
ratings.gen %>% group_by(threat, threat_num, threat_both, diagnostic, subject) %>% summarise(rating.m=mean(rating, na.rm=T), rating.se=se(rating, na.rm=T)) %>%
  summarise(rating.se=se(rating.m, na.rm=T), rating.m=mean(rating.m, na.rm=T)) %>% 
  ggplot(aes(x=threat, y=rating.m, color=diagnostic, shape=diagnostic, group=diagnostic)) +
  geom_point(position=dodge, size=3) + geom_line(position=dodge) +
  geom_errorbar(aes(ymin=rating.m-rating.se*1.96, ymax=rating.m+rating.se*1.96), position=dodge, width=dodge.width) +
  ylab("Threat Rating (1-5)") + xlab("Threat") + labs(color="Diagnostic Region", shape="Diagnostic Region") + myGgTheme

# #mean scores generalization phase by sex
# ratings.gen %>% group_by(threat, sex, subject) %>% summarise(rating.m=mean(rating, na.rm=T), rating.se=se(rating, na.rm=T)) %>%
#   group_by(threat, sex) %>% summarise(rating.se=se(rating.m, na.rm=T), rating.m=mean(rating.m, na.rm=T)) %>% 
#   ggplot(aes(x=threat, y=rating.m, color=sex, shape=sex)) +
#   geom_point(position=dodge) +
#   geom_errorbar(aes(ymin=rating.m-rating.se*1.96, ymax=rating.m+rating.se*1.96), position=dodge, width=dodge.width) + myGgTheme

#mean scores generalization phase by diagnostic and sex (i.e. pair variable)
ratings.gen %>% group_by(threat, pair, subject) %>% summarise(rating.m=mean(rating, na.rm=T), rating.se=se(rating, na.rm=T)) %>%
  summarise(rating.se=se(rating.m, na.rm=T), rating.m=mean(rating.m, na.rm=T)) %>% 
  ggplot(aes(x=threat, y=rating.m, color=pair, shape=pair, group=pair)) +
  geom_point(position=dodge, size=3) + geom_line(position=dodge) +
  geom_errorbar(aes(ymin=rating.m-rating.se*1.96, ymax=rating.m+rating.se*1.96), position=dodge, width=dodge.width) +
  ylab("Threat Rating (1-5)") + xlab("Threat") + labs(color="Pair", shape="Pair") +
  scale_shape_manual(values=c(0:2, 6), labels=c("male eyes", "male mouth/nose", "female eyes", "female mouth/nose")) +
  scale_color_discrete(labels=c("male eyes", "male mouth/nose", "female eyes", "female mouth/nose")) + myGgTheme


# Inference Tests ---------------------------------------------------------
ratings.subject = ratings %>% group_by(subject, threat, phase) %>% 
  summarise(rating.se = se(rating, na.rm=T), rating = mean(rating, na.rm=T))
ratings.subject %>% filter(phase == "Hab") %>% t.test(rating ~ threat, ., paired=T) %>% apa::t_apa(es_ci=T)
ratings.subject %>% filter(phase == "Hab") %>% group_by(threat) %>% summarise(rating.m = mean(rating, na.rm=T), rating.sd = sd(rating, na.rm=T), rating.se = se(rating, na.rm=T))
# d = .5 #effect size that constitutes equivalence bounds
# ratings.subject %>% filter(phase == "Hab") %>% select(-rating.se) %>% spread(threat, rating) %>%
#   TOSTER::dataTOSTpaired(pairs=list(c(i1="CS+", i2="CS-")), low_eqbound=-d, high_eqbound=d, plots=F, desc=F) #for Hab, equivalence hypothesis => TOST test (see https://rpsychologist.com/d3/equivalence/)

ratings.subject %>% filter(phase == "Acq") %>% t.test(rating ~ threat, ., paired=T) %>% apa::t_apa(es_ci=T)
ratings.subject %>% filter(phase == "Acq") %>% group_by(threat) %>% summarise(rating.m = mean(rating, na.rm=T), rating.sd = sd(rating, na.rm=T), rating.se = se(rating, na.rm=T))

# ANOVA Generalization Phase (SPAI)
ratings.subject.gen.diagnostic = ratings.gen %>% 
  group_by(subject, SPAI, SPAI.z, STAI, STAI.z, threat, diagnostic, pairs) %>% 
  summarise(rating.m=mean(rating, na.rm=T), rating.se=se(rating, na.rm=T))
ez::ezANOVA(data=ratings.subject.gen.diagnostic, 
            dv=.(rating.m), wid=.(subject), 
            within=.(threat, diagnostic), 
            #between=.(pairs),
            between=.(SPAI.z), observed=SPAI.z,
            detailed=T, type=2) %>% apa::anova_apa(force_sph_corr=T)

ez::ezANOVA(data=ratings.subject.gen.diagnostic, 
            dv=.(rating.m), wid=.(subject), 
            within=.(threat, diagnostic), 
            #between=.(pairs),
            between=.(STAI.z), observed=STAI.z,
            detailed=T, type=2) %>% apa::anova_apa(force_sph_corr=T)


# threat main effect: t-tests of adjacent threat-levels
ratings.subject.gen = ratings.gen %>% group_by(subject, threat, threat_num) %>% summarise(rating.m=mean(rating, na.rm=T), rating.se=se(rating, na.rm=T))
ratings.subject.first.gen = ratings.first.gen %>% group_by(subject, threat, threat_num) %>% summarise(rating.m=mean(rating, na.rm=T), rating.se=se(rating, na.rm=T))  # for first half
for (i in (min(ratings.subject.gen$threat_num)+1):max(ratings.subject.gen$threat_num)) {
  levels = c(i-1, i)
  cat(paste0("\n\nComparing levels: ", paste(levels, collapse=" vs. "), "\n"))
  ratings.subject.gen %>% filter(threat_num %in% levels) %>%
    t.test(rating.m ~ threat, ., paired=T) %>% apa::t_apa(es_ci=T)
}

# marginal diagnostic main effect (descriptives because dichotomous)
ratings.gen %>% group_by(diagnostic, subject) %>% summarise(rating = mean(rating, na.rm=T)) %>% summarise(rating.m = mean(rating, na.rm=T), rating.sd = sd(rating, na.rm=T), rating.se = se(rating, na.rm=T)) %>% arrange(desc(rating.m))

# t-tests of diagnostic regions within threat-levels
# ratings.subject.gen.diagnosticXthreat = ratings.gen %>% group_by(subject, SPAI, STAI, threat, threat_num, diagnostic) %>% summarise(rating.m=mean(rating, na.rm=T), rating.se=se(rating, na.rm=T))
# for (i in 1:6) {
#   cat(paste0("\n\nDiagnostic regions in level: ", i, "\n"))
#   ratings.subject.gen.diagnosticXthreat %>% filter(threat_num == i) %>% 
#     t.test(rating.m ~ diagnostic, ., paired=T) %>% apa::t_apa(es_ci=T)
# }

#SPAI main effect (i.e., correlation)
ratings.subject.gen.lvl = ratings.subject.gen.diagnostic %>% group_by(subject, SPAI, STAI) %>% summarise(rating.se = se(rating.m, na.rm=T), rating.m = mean(rating.m, na.rm=T))
ratings.subject.gen.lvl %>% with(cor.test(rating.m, SPAI, alternative="greater")) %>% apa::cor_apa(r_ci=T) #%>% correlation_out()
ratings.subject.gen.lvl %>% ggplot(aes(x=SPAI, y=rating.m, color=SPAI, fill=SPAI)) +
  geom_errorbar(aes(ymin=rating.m-rating.se*1.96, ymax=rating.m+rating.se*1.96), width=spai.width) +
  stat_smooth(method="lm", color = "black") +
  #geom_point(size=4, shape=21, color="black") +
  geom_point(size=4) + 
  ylab("Average Rating (1-5)") +
  scale_color_viridis_c() + scale_fill_viridis_c() + myGgTheme + theme(legend.position = "none")


#STAI main effect (i.e., correlation)
ratings.subject.gen.lvl = ratings.subject.gen.diagnostic %>% group_by(subject, SPAI, STAI) %>% summarise(rating.se = se(rating.m, na.rm=T), rating.m = mean(rating.m, na.rm=T))
ratings.subject.gen.lvl %>% with(cor.test(rating.m, STAI, alternative="greater")) %>% apa::cor_apa(r_ci=T) #%>% correlation_out()
ratings.subject.gen.lvl %>% ggplot(aes(x=STAI, y=rating.m, color=STAI, fill=STAI)) +
  geom_errorbar(aes(ymin=rating.m-rating.se*1.96, ymax=rating.m+rating.se*1.96), width=.4) +
  stat_smooth(method="lm", color = "black") +
  #geom_point(size=4, shape=21, color="black") +
  geom_point(size=4) + 
  ylab("Average Rating (1-5)") +
  scale_color_viridis_c() + scale_fill_viridis_c() + myGgTheme + theme(legend.position = "none")

#SPAI * threat * diagnostic
# ratings.subject.gen.diagnostic %>% ggplot(aes(x=SPAI, y=rating.m, color=diagnostic, fill=diagnostic, shape=diagnostic)) +
#   facet_wrap(vars(threat)) +
#   geom_smooth(method="lm", color = "black") + geom_smooth(method="lm", color = "black", fill=NA) + #print regression lines on top of confidence bands
#   geom_point(size=4) +
#   ylab("Average Rating (1-5)") +
#   myGgTheme + myGgTheme + theme(legend.position = "none")
  
# ratings.subject.gen.diagnostic %>% ggplot(aes(x=SPAI, y=rating.m, color=threat, fill=threat, shape=threat)) +
#   facet_wrap(vars(diagnostic)) +
#   geom_smooth(method="lm", color = "black") + geom_smooth(method="lm", color = "black", fill=NA) + #print regression lines on top of confidence bands
#   #geom_point(size=4) +
#   ylab("Average Rating (1-5)") +
#   scale_color_manual(values=colors) + scale_fill_manual(values=colors) + 
#   guides(colour=guide_legend(reverse=T), shape=guide_legend(reverse=T), fill=guide_legend(reverse=T)) + 
#   myGgTheme + myGgTheme #+ theme(legend.position = "none")
#diagnostic eyes: low SPAI <=> more differentiation
#diagnostic m/n:  low SPAI <=> less fear generalization (but similar differentiation)

# ANOVA Generalization Phase per Block 
ratings.subject.gen.block = ratings.gen %>% 
  mutate(block = as_factor(block)) %>%
  group_by(subject, block, SPAI, SPAI.z, STAI, STAI.z, threat, diagnostic, pairs) %>% 
  summarise(rating.m=mean(rating, na.rm=T), rating.se=se(rating, na.rm=T))
#ratings.subject.gen.block %>% filter(rating.m %>% is.na()) #subject 57: all ratings of CS+ NA in block 3
ratings.subject.gen.block = ratings.subject.gen.block %>% filter(subject!=57)
ez::ezANOVA(data=ratings.subject.gen.block, 
            dv=.(rating.m), wid=.(subject), 
            within=.(threat, diagnostic, block), 
            #between=.(pairs),
            #between=.(SPAI.z), observed=SPAI.z,
            between=.(STAI.z), observed=STAI.z,
            detailed=T, type=2) %>% apa::anova_apa(force_sph_corr=T)

# Block Main Effect
ratings.subject.gen.block %>% 
  group_by(block)%>%
  summarise(mean = mean(rating.m),
            sd = sd(rating.m))

# Block x Threat Interaction
ratings.gen.threat.block <- ratings.subject.gen.block %>% 
  group_by(block, threat)%>%
  summarise(mean = mean(rating.m),
            se = se(rating.m))
ratings.gen.threat.block

print(ratings.gradient.plot.blocks <- ratings.gen.threat.block %>% 
        ggplot(aes(x=threat, y=mean, color=block, group=block)) + 
        #try ggbeeswarm ? or position_jitter ? or geom_quasirandom ? or violinplot?
        #geom_path(data=ratings.ga.gen %>% filter(threat_num %in% c(1, threatLevels.n)), color = "black", size=1.5) + #generalization line
        geom_line(size=1) + geom_point(size=4.5) + 
        geom_errorbar(aes(ymin=mean-se*1.96, ymax=mean+se*1.96), size=1.5) +
        scale_colour_discrete(name = "Generalization", labels = c("first half", "second half"))+
        ylab("Threat Rating (1-5)") + xlab("Threat") + #labs(color="Generalization", fill="Generalization") +
        myGgTheme + theme(legend.position="right"))
#CS+ goes up, all other stimuli down

# Block x Threat x Diagnostic Interaction
ratings.gen.diagnostic.block <- ratings.subject.gen.block %>% 
  group_by(block, threat, diagnostic)%>%
  summarise(mean = mean(rating.m),
            se = se(rating.m)) %>%
  mutate(
    group = case_when(block == 2 & diagnostic == "Eyes" ~ 1,
                      block == 2 & diagnostic == "Mouth/Nose" ~ 2,
                      block == 3 & diagnostic == "Eyes" ~ 3,
                      block == 3 & diagnostic == "Mouth/Nose" ~ 4) %>% as_factor())
ratings.gen.diagnostic.block

print(ratings.gradient.plot.blocks <- ratings.gen.diagnostic.block %>% 
        #ggplot(aes(x=threat, y=mean, color=group , group=group)) + 
        ggplot(aes(x=threat, y=mean, color=block, group=block)) + facet_wrap(vars(diagnostic)) +
        #try ggbeeswarm ? or position_jitter ? or geom_quasirandom ? or violinplot?
        #geom_path(data=ratings.ga.gen %>% filter(threat_num %in% c(1, threatLevels.n)), color = "black", size=1.5) + #generalization line
        geom_line(size=1) + geom_point(size=4.5) + 
        geom_errorbar(aes(ymin=mean-se*1.96, ymax=mean+se*1.96), size=1.5) +
        #scale_colour_discrete(name = "Block & Diagnostic Region", labels = c("first eyes", "first mouth/nose", "second eyes", "second mouth/nose"))+
        scale_colour_discrete(name = "Generalization", labels = c("first half", "second half"))+
        ylab("Threat Rating (1-5)") + xlab("Threat") + #labs(color="Block", fill="Block") +
        myGgTheme + theme(legend.position="right"))
#biggest changes between blocks: Eyes: GS4, Mouth/Nose: GS2 & 3
#TODO pair-wise t-tests first vs. second half across all diagnostic x threat

# ANOVA First Generalization Phase
ratings.subject.first.gen.diagnostic = ratings.first.gen %>% group_by(subject, SPAI, SPAI.z, STAI, STAI.z, threat, diagnostic, pairs) %>%
  summarise(rating.m=mean(rating, na.rm=T), rating.se=se(rating, na.rm=T))
# ez::ezANOVA(data=ratings.subject.first.gen.diagnostic,
#             dv=.(rating.m), wid=.(subject),
#             within=.(threat, diagnostic),
#             #between=.(pairs),
#             between=.(SPAI.z), observed=SPAI.z,
#             #between=.(STAI.z), observed=STAI.z,
#             detailed=T, type=2) %>% apa::anova_apa(force_sph_corr=T)

# Gradient Analysis -----------------------------------------------------------
individualPlots = F

#subject-level means generalization phase
gradients.simple = data.frame(subject=numeric(), lds=integer(), diff=numeric(), level=numeric(), stringsAsFactors=F)
for (s in ratings.subject.gen$subject %>% unique()) {
  ratings.temp = ratings.subject.gen %>% filter(subject==s)
  
  gradient = gradient.analysis(ratings.temp$rating.m)
  gradients.simple = bind_rows(gradients.simple, data.frame(subject=s, lds=gradient["lds"], diff=gradient["diff"], level=gradient["level"], row.names = NULL, stringsAsFactors=F))
  
  if (individualPlots) {
    gradient = gradient %>% signif(3)
    print(ratings.temp %>% ggplot(aes(x=threat_num, y=rating.m, color=threat, group=NA)) + #generalization line
            geom_path(data=ratings.temp %>% filter(threat_num %in% c(1, 6)), aes(group=NA), color = "black", size=1.5) +
            geom_errorbar(aes(ymin=rating.m-rating.se*1.96, ymax=rating.m+rating.se*1.96)) +
            geom_line() + geom_point() + 
            scale_color_manual(values=colors) +
            ggtitle(paste0(s, ": lds = ", gradient[1], ", diff = ", gradient[2], ", level = ", gradient[3])) + myGgTheme)
    invisible(readline(prompt="Press [enter] to continue")) #TODO save to file instead
  }
}

names(gradients.simple) = c("subject", "Gen_all_lds", "Gen_all_diff", "Gen_all_level")

#subject-level means generalization phase by diagnostic
gradients.diagnostic = data.frame()
for (s in ratings.subject.gen.diagnostic$subject %>% unique()) {
  ratings.temp = ratings.subject.gen.diagnostic %>% filter(subject==s)
  
  gradient.eyes = "Eyes" %>% {filter(ratings.temp, diagnostic==.)} %>% .$rating.m %>% gradient.analysis()
  gradient.m_n = "Mouth/Nose" %>% {filter(ratings.temp, diagnostic==.)} %>% .$rating.m %>% gradient.analysis()
  gradients.diagnostic = rbind(gradients.diagnostic, c(s, gradient.eyes, gradient.m_n))
  
  if (individualPlots) {
    print(ratings.temp %>% ggplot(aes(x=threat, y=rating.m, color=diagnostic, shape=diagnostic, group=diagnostic)) + 
            #geom_path(data=ratings.ga.m.diagnostic %>% filter(threat_num %in% c(1, 6)), aes(group=NA), color = "black", size=1.5) + #generalization line
            geom_errorbar(aes(ymin=rating.m-rating.se*1.96, ymax=rating.m+rating.se*1.96), position=dodge, width=dodge.width) + 
            geom_line(position=dodge) + geom_point(position=dodge) + 
            ggtitle(s) + myGgTheme)
    invisible(readline(prompt="Press [enter] to continue")) #TODO save to file instead
  }
}
names(gradients.diagnostic) = c("lds", "diff", "level") %>% {c("subject", paste0("Gen_eyes_", .), paste0("Gen_mn_", .))}

gradients = full_join(gradients.simple, gradients.diagnostic, by="subject")

# #subject-level means first half generalization phase
# gradients.simple.first.half = data.frame(subject=character(), lds=numeric(), diff=numeric(), level=numeric(), stringsAsFactors=F)
# for (s in ratings.subject.first.gen$subject %>% unique()) {
#   ratings.temp = ratings.subject.first.gen %>% filter(subject==s)
# 
#   gradient = gradient.analysis(ratings.temp$rating.m)
#   gradients.simple.first.half = bind_rows(gradients.simple.first.half, data.frame(subject=s %>% as.character(), lds=gradient["lds"], diff=gradient["diff"], level=gradient["level"], row.names = NULL, stringsAsFactors=F))
# 
#   if (individualPlots) {
#     gradient = gradient %>% signif(3)
#     print(ratings.temp %>% ggplot(aes(x=threat_num, y=rating.m, color=threat, group=NA)) + #generalization line
#             geom_path(data=ratings.temp %>% filter(threat_num %in% c(1, 6)), aes(group=NA), color = "black", size=1.5) +
#             geom_errorbar(aes(ymin=rating.m-rating.se*1.96, ymax=rating.m+rating.se*1.96)) +
#             geom_line() + geom_point() +
#             scale_color_manual(values=colors) +
#             ggtitle(paste0(s, ": lds = ", gradient[1], ", diff = ", gradient[2], ", level = ", gradient[3])) + myGgTheme)
#     invisible(readline(prompt="Press [enter] to continue")) #TODO save to file instead
#   }
# }
# names(gradients.simple.first.half) = c("subject", "Gen_all_lds", "Gen_all_diff", "Gen_all_level")
# 
# #subject-level means first half generalization phase by diagnostic
# gradients.diagnostic.first.half = data.frame()
# for (s in ratings.subject.first.gen.diagnostic$subject %>% unique()) {
#   ratings.temp = ratings.subject.first.gen.diagnostic %>% filter(subject==s)
# 
#   gradient.eyes = "Eyes" %>% {filter(ratings.temp, diagnostic==.)} %>% .$rating.m %>% gradient.analysis()
#   gradient.m_n = "Mouth/Nose" %>% {filter(ratings.temp, diagnostic==.)} %>% .$rating.m %>% gradient.analysis()
#   gradients.diagnostic.first.half = rbind(gradients.diagnostic.first.half, c(s, gradient.eyes, gradient.m_n))
# 
#   if (individualPlots) {
#     print(ratings.temp %>% ggplot(aes(x=threat, y=rating.m, color=diagnostic, shape=diagnostic, group=diagnostic)) +
#             #geom_path(data=ratings.ga.m.diagnostic %>% filter(threat_num %in% c(1, 6)), aes(group=NA), color = "black", size=1.5) + #generalization line
#             geom_errorbar(aes(ymin=rating.m-rating.se*1.96, ymax=rating.m+rating.se*1.96), position=dodge, width=dodge.width) +
#             geom_line(position=dodge) + geom_point(position=dodge) +
#             ggtitle(s) + myGgTheme)
#     invisible(readline(prompt="Press [enter] to continue")) #TODO save to file instead
#   }
# }
# names(gradients.diagnostic.first.half) = c("lds", "diff", "level") %>% {c("subject", paste0("Gen_eyes_", .), paste0("Gen_mn_", .))}
# 
# gradients.diagnostic.first.half$subject <- as.character(gradients.diagnostic.first.half$subject)
# gradients.first.half = full_join(gradients.simple.first.half, gradients.diagnostic.first.half, by="subject")
# gradients.first.half$subject <- as.numeric(gradients.first.half$subject)

# Wide format for correlations ----------------------------------------------------
#all phases except generalization
ratings.wide = ratings %>% filter(phase != "Gen", !(subject %in% exclusions.rating)) %>% 
  group_by(subject, phase, threat) %>% summarise(rating.m = mean(rating, na.rm=T))
ratings.wide$threat = recode(ratings.wide$threat, `1` = "CS-", `6` = "CS+")
ratings.wide = ratings.wide %>% unite(temp, phase, threat) %>% spread(temp, rating.m)

#all ratings in generalization
ratings.wide.threat = ratings.subject.gen %>% select(-rating.se, -threat_num) %>% spread(threat, rating.m)
names(ratings.wide.threat) = c("subject", "Gen_all_CS-", paste0("Gen_all_GS", 1:4), "Gen_all_CS+")
ratings.wide = full_join(ratings.wide, ratings.wide.threat, by="subject")

# #ratings in first half of generalization
# ratings.wide.first.threat = ratings.subject.first.gen %>% select(-rating.se, -threat_num) %>% spread(threat, rating.m)
# names(ratings.wide.first.threat) = c("subject", "Gen_all_CS-", paste0("Gen_all_GS", 1:4), "Gen_all_CS+")
# ratings.first.wide = full_join(ratings.wide, ratings.wide.first.threat, by="subject")
# 
# #gradients first half
# ratings.first.wide = full_join(ratings.wide, gradients.first.half, by="subject") %>% 
#   select("subject", contains("Hab"), contains("Acq"), contains("Gen"), everything()) %>% tibble()

#gradients
ratings.wide = full_join(ratings.wide, gradients, by="subject") %>% 
  select("subject", contains("Hab"), contains("Acq"), contains("Gen"), everything()) %>% tibble()

# #gradients first half
# ratings.first.wide = full_join(ratings.first.wide, gradients.first.half, by="subject") %>%
#   select("subject", contains("Hab"), contains("Acq"), contains("Gen"), everything()) %>% tibble()

#all(ratings.wide == read_rds("ratings.wide.rds" %>% paste0(path.rds, .)), na.rm=T) #check equivalence of processing
#ratings.wide %>% write_rds("ratings.wide.rds" %>% paste0(path.rds, .))
