if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)

requirePackage = function(name, load=T) {
  package = as.character(name)
  if (package %in% rownames(installed.packages()) == FALSE) install.packages(package)
  if (load) library(package, character.only=TRUE)
}

{ #install packages needed
  requirePackage("scales", load=F)
  requirePackage("cowplot", load=F)
  requirePackage("readxl", load=F)
  requirePackage("apaTables", load=F)
  requirePackage("schoRsch", load=F)
  requirePackage("apa", load=F)
  requirePackage("ez", load=F)
  requirePackage("TOSTER", load=F)
  requirePackage("psych", load=F)
}


{ # Variables ---------------------------------------------------------------
  # a priori exclusions for all variables
  exclusions = c(
     5, #pain rating < 4
    40, #different faces between discrimination training and main experiment
    60, #test subject for new experimenter (paradigm not performed)
    88, #experiment aborted after discrimination task
    91  #pain rating < 4
  ) %>% unique() %>% sort()
  
  exclusions.onlyGen = c(
    44, 50 #restart during acquisition => exclude habituation & acquisition (but keep generalization)
  ) %>% unique() %>% sort()
  
  preAcqEnd = 20.5 #no pain during first 20 trials
  acqEnd = preAcqEnd + 40 #acquisition for 40 trials
  gen1End = acqEnd + 70
  trials.n = 200 #number of trials that shall be analyzed (if more trials, last ones will be taken)
  
  recalibration = c(acqEnd, gen1End) #border(s) separating trials
  
  preStim = 2000
  #expositionStart = 0
  ratingStart = 4000
  shockTime = ratingStart + 1850
  shockEnd = shockTime + 150
  trialEnd = ratingStart + 2000
  itiEnd = trialEnd + c(3150, 9500) #not theoretical but empirical bounds
  
  threatLevels.n = 6
  
  sample.rate = 500 #samples/second
  trial.duration = sample.rate * trialEnd / 1000 #seconds
  
  #startID = "CONDITION" #identifier for trial start messages
  expoID = "Stimulus " #identifier for exposition start messages
  stimExt = ".jpg"
  
  screen.height = 1080 #height of screen in pix
  screen.width  = 1920 # width of screen in pix
}

# Paths -------------------------------------------------------------------
path = "C:/Users/jat41gk/Documents/Projekte/Visual Exploration Social Anxiety/" #Janna
#path = "C:/Data/C10.3.6 Fear Generalization Anxiety/" #Mario @work
#path = sub("C:/Data", "D:/Arbeit", path, fixed=TRUE) #Mario @home

#load behavioral data (logs)
path.ratings = "Data/log/" %>% paste0(path, .)
files.rating.prefix = "vp"
files.rating.extension = "-LookAtMe.log"

path.eye = "Data/EyeLink/Output/" %>% paste0(path, .) #eye tracking data
path.pupil  = "Data/Pupil/" %>% paste0(path, .)

path2 = "C:/Users/mar84qk/Dropbox/Arbeit/C10 - Gamer Fear/3 Diagnostic Generalization/" #@work
#path2 = sub("C:/Users/mar84qk", "D:", path2, fixed=TRUE) #@home

path.rois = "2 Material/Look At Me 3 Anx/Stimuli/mask/" %>% paste0(path2, .)

#path.rds = "Analyse/" %>% paste0(path, .) #Janna
#path.rds = "3 Analysis/main Anx/" %>% paste0(path2, .) #Mario
path.rds = getwd() %>% paste0("/") #sync to git

path.phys = "Data/physio/" %>% paste0(path, .)

path.rpeaks = "rpeaks/" %>% paste0(path.phys, .)
path.rpeaks.postfix = "_rpeaks.csv"

{ # Functions ---------------------------------------------------------------
  se = function(x, na.rm = FALSE) {
    sd(x, na.rm) / sqrt(if(!na.rm) length(x) else sum(!is.na(x)))
  }
  
  correlation_out = function(coroutput) {
    names = coroutput$data.name %>% strsplit(" and ") %>% unlist()
    cat(paste0("r(", names[1], ", ", names[2], "): ", coroutput %>% apa::cor_apa(print=F)), "\n")
  }
  
  F_apa = function(x) {
    cat("F(", paste(x$parameter, collapse=", "), ") = ", 
        round(x$statistic, 2), 
        ", p ",
        ifelse(x$p.value < .001, "< .001",
               paste0("= ", round(x$p.value, 2))),
        "\n", sep="")
  }
  
  ez.ci = function(ez, conf.level = .95, sph.corr=T) {
    for (effect in ez$ANOVA$Effect) {
      index.sphericity = which(ez$`Sphericity Corrections`$Effect == effect)
      GGe = ifelse(sph.corr==F || length(index.sphericity)==0, 1, ez$`Sphericity Corrections`$GGe[index.sphericity])
      
      index.effect = which(ez$ANOVA$Effect == effect)
      ez$ANOVA %>% with(apaTables::get.ci.partial.eta.squared(F[index.effect], DFn[index.effect]*GGe, DFd[index.effect]*GGe, conf.level = conf.level)) %>% 
        sapply(round, digits=2) %>% 
        paste0(collapse=", ") %>% paste0(effect, ": ", round(conf.level*100), "% CI [", ., "]\n") %>% cat()
    }
  }
  
  lmer.ci = function(lmer, conf.level = .95, twotailed=T) {
    values = lmer %>% summary() %>% .$coefficients %>% .[, c("Estimate", "df")] %>% data.frame()
    effects = values %>% rownames()
    
    for (i in seq(effects)) {
      psych::r.con(values$Estimate[i], values$df[i], p = conf.level, twotailed=twotailed) %>% 
        round(digits=2) %>% 
        paste0(collapse=", ") %>% paste0(effects[i], ": ", round(conf.level*100), "% CI [", ., "]\n") %>% cat()
    }
  }
  
  color.gradient.discrete = function(color.low, color.high, n) {
    scales::seq_gradient_pal(low=color.low, high=color.high)(seq(0, 1, length.out = n)) 
    #ggplot2::scale_color_gradient ? doesn't work well :/ also not ggplot2::scale_color_gradientn
  }
  
  color.gradient.discrete.single = function(color, n, color.low="white", includeEnd=F) {
    if (includeEnd) return(color.gradient.discrete(color.low, color, n))
    else return(color.gradient.discrete(color.low, color, n+1)[-1])
  }
  
  color.gradient.discrete.divergent = function(color.low, color.high, n, color.middle = "white") {
    if (n %% 2 == 0) { #even number => don't include middle color
      c(color.gradient.discrete.single(color.low, n/2, color.middle, includeEnd=F) %>% rev(),
        color.gradient.discrete.single(color.high, n/2, color.middle, includeEnd=F)
      ) %>% return()
    } else { #odd number => round up but include middle color only once
      c(color.gradient.discrete.single(color.low, ceiling(n/2), color.middle, includeEnd=T) %>% rev(),
        color.gradient.discrete.single(color.high, ceiling(n/2), color.middle, includeEnd=T)
      ) %>% unique() %>% return()
    }
  }
  
  gradient.analysis = function(x) {
    n = length(x)
    level = mean(x)
    diff = x[n] - x[1]
    lds = mean(x[c(1, n)]) - mean(x[2:(n-1)])
    result = c(lds, diff, level); names(result) = c("lds", "diff", "level")
    return(result)
  }
  
  lds = function(x, standardize=F) {
    n = length(x)
    lds = mean(x[c(1, n)]) - mean(x[2:(n-1)])
    
    if (standardize) lds = lds / abs(x[n] - x[1])
    return(lds)
  }
  
  pathToCode = function(path, path.sep="/", file.ext="\\.") {
    first = path %>% gregexpr(path.sep, .) %>% sapply(max) %>% {. + 1}
    last = path %>% gregexpr(file.ext, .) %>% sapply(max) %>% {. - 1}
    return(path %>% substring(first, last))
  }
  
  codeToNum = function(code) code %>% gsub("\\D+", "", .) %>% as.integer()
  
  read.phys = function(path) {
    #Problem: subjects until 36 have first column "min", from 37 this column is not present
    
    #read.delim(path, na.strings="", skip=9) %>%
    #read_table(path, na="", skip=9, progress=F, show_col_types=F) %>%
    read_delim(path, delim="\t", na="", skip=9, progress=F, show_col_types=F) %>% suppressMessages() %>% 
      rename(EDA = "CH1", ECG = "CH2", Trigger = "CH28") %>%
      #read_delim(path, delim="\t", na="", skip=11, col_names=c("min", "EDA", "ECG", "Trigger", "blank"), progress=F, show_col_types=F) %>% 
      
      filter(Trigger <= 2^8)
  }
  
  get.trigger.onsets = function(x) {
    return(x %>% diff() %>% {. > 0} %>% which() %>% {. + 1})
    
    # triggers.all = which(x != 0) #all trigger indices (including trains)
    # triggers = c(first(triggers.all), #take very first instance
    #              triggers.all[triggers.all != lag(triggers.all)+1]) #and every instance that is not the same as its precursor index + 1
    # return(triggers[!is.na(triggers)])
  }
  
  recode.triggers = function(x, triggercode=1) {
    recode = rep.int(0, times=length(x))
    recode[get.trigger.onsets(x)] = triggercode
    return(recode)
  }
  
  expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
}

# plotting ----------------------------------------------------------------
colors = color.gradient.discrete("blue", "red", n=threatLevels.n)
#scales::show_col(colors)
dodge.width = .6 #good for factors on x-axis
dodge = position_dodge(width=dodge.width)
dodge_half = position_dodge(width=dodge.width/2)

dodge_stai = position_dodge(width=1)
spai.width = .1

#ggplot general theme
#theme_set( #don't use theme_set because it has to be executed every session but myGgTheme can be saved in global environment
myGgTheme <- theme_bw() + theme(
  #aspect.ratio = 1,
  plot.title = element_text(hjust = 0.5),
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
#)

# Files -------------------------------------------------------------------
files.rating = list.files(path.ratings, pattern=paste0("^", files.rating.prefix, ".*", files.rating.extension, "$"), full.names=TRUE)
files.phys = list.files(path.phys, pattern=".*.txt") #load physiology

questionnaires = suppressMessages( #avoid name repair message (affected columns will be deselected anyway)
  readxl::read_excel("Daten_Fragebögen.xlsx" %>% paste0(path, "Data/", .),
                     sheet = "Daten_FB_gesamt", skip=2, col_names=T)) %>% select(1:5) %>% rename(subject = `VP-Nr`) %>% 
  mutate(subject = subject %>% gsub("vp", "", ., fixed=T) %>% as.integer(),
         SPAI = SPAI / 22, #German SPAI version: mean item scores, not sum scores (for comparison with cut-off values)
         discr = Screening_nachher - Screening_vorher,
         problem = abs(discr) > 1.5, problem = ifelse(problem %>% is.na(), F, problem)) %>% #preregistration
  filter(subject %>% is.na() == F, subject > 0, #get rid of test subject & empty rows
         subject %in% exclusions == F) #apply a priori exclusions

spai.cutoff = 2.79 #Glischinski et al. (2018) for German version of SPAI (N = 359 patients)
questionnaires %>% summarize(mean = mean(SPAI), SD = sd(SPAI))
paste0("subjects meeting the cut-off: ", {mean(questionnaires$SPAI >= spai.cutoff)*100} %>% round(digits=2), "% ", 
       "(N = ", sum(questionnaires$SPAI > spai.cutoff), "; ",
       "z = ", with(questionnaires, (spai.cutoff - mean(SPAI)) / sd(SPAI)) %>% round(digits=2), ")")

questionnaires %>%
  ggplot(aes(x=SPAI)) + geom_histogram(binwidth=.25, boundary=spai.cutoff, color="black", fill="grey") + 
  geom_vline(xintercept=spai.cutoff, color="red") + myGgTheme + scale_y_continuous(breaks=scales::breaks_pretty())

#questionnaires %>% filter(problem==T)
#exclusions = exclusions %>% c(questionnaires %>% filter(problem==T) %>% .$subject) #rather don't exclude (deviation from preregistration)

#questionnaires %>% select(SPAI, STAI) %>% summarise(across(.fns = function(x) {max(x) - min(x)})) %>% mutate(width = SPAI/STAI) #=> round spai.width to .1 (see above)

