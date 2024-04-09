if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
#source("0 General.R" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/C10 - Gamer Fear/3 Diagnostic Generalization/3 Analysis/main Anx/", .))

requirePackage("signal", load=F)

exclusions.eda.num = c(exclusions) %>% 
  c() #exclusion by visual inspection is obsolete by ucr scoring
exclusions.eda = exclusions.eda.num %>% toPhysString()

#ratings = readRDS("ratings.rds" %>% paste0(path.rds, .))

#output plots to check if everything works as expected
checkPlots = F #check plot after downsampling, smoothing, trimming?
checkPlotsTrial = T #check plots for trial scoring of SCRs

#downsampling
downsampling = T
sample.rate.new = 20 #sample down to (Hz)
if (sample.rate %% sample.rate.new != 0) warning(paste0("sample.rate must be a multiple of sample.rate.new. Rounding to sample.rate.new = ", sample.rate / round(sample.rate / sample.rate.new), " Hz"))

#smoothing
smoothing = T
lowPassFreq = 2 #low pass filter (Hz)

#ucr analysis
ucrPlots = T
ucrMinWindow = c(.5, 3) #within how many seconds after UCS must the detected minimum appear?
#ucrLatestPeak = 10 #how many seconds after UCS must the detected peak have been appeared at the latest? (not used anymore)
ucrMin = .1 #minimum amplitude of an individual unconditioned response (in µS)

ucrMeanValid = .5
ucrMeanAmp = .2

#cr analysis
crPlots = T
checkPlotsTrialGenOnly = T
crMinWindow = ucrMinWindow #within how many seconds after CS must the detected minimum appear?
crMin = .02 #minimum amplitude of an individual conditioned response (in µS)

#save waveform into png
createSclPlots = F
saveSclPlotsToFile = T && createSclPlots
eda.middleplot.time = sample.rate * 60 #in seconds, used to create middle plots
ecg.middleplot.time = sample.rate * 10 #ecg plot  not needed due to detect beats script giving better impression about data quality

#save files
saveLedalab = F #shall downsampled, smoothened, and trimmed signal be saved to file (e.g. for ledalab)?
analysis.border = ifelse(downsampling, sample.rate.new, sample.rate) * 30 #trim around triggers (plus border in seconds)

{ # Functions ---------------------------------------------------------------
  sample.down = function(signal, conversionRate) {
    suppressWarnings(matrix(signal, ncol=conversionRate, byrow=T)) %>% apply(1, mean) %>% 
      head(-1) #discard last sample because it gets distorted by zero padding
  }
  
  closestRepresentative = function(x, values, returnIndices=F) {
    indices = x %>% sapply(function(x, data) which.min(abs(data-x)), data=values)
    if (returnIndices) return(indices)
    return(values[indices])
  }
  
  localMaxima = function(signal, integerMode=F) {
    result = signal %>% c(-ifelse(integerMode, .Machine$integer.max, Inf), .) %>% #prepend SMALLEST value to allow that first sample point is also local maximum
      diff() %>% {. > 0} %>% #which sample points are rising compared to their predecessor?
      rle() %>% .$lengths %>% #for how long are segments rising consecutively?
      cumsum() %>% #cumulative sum to find indices of local extrema
      {.[seq.int(from=1, to=length(.), by=2)]} #only take every other point
    if (signal[[1]] == signal[[2]]) result = result[-1]
    
    return(result)
  }
  
  localMinima = function(signal, integerMode=F) {
    result = signal %>% c(ifelse(integerMode, .Machine$integer.max, Inf), .) %>% #prepend BIGGEST value to allow that first sample point is also local maximum
      diff() %>% {. > 0} %>% #which sample points are rising compared to their predecessor?
      rle() %>% .$lengths %>% #for how long are segments rising consecutively?
      cumsum() %>% #cumulative sum to find indices of local extrema
      {.[seq.int(from=1, to=length(.), by=2)]} #only take every other point
    if (signal[[1]] == signal[[2]]) result = result[-1]
    
    return(result)
  }
  
  inflectionPoints = function(signal, integerMode=F, upOnly=T) {
    extrema = signal %>% c(ifelse(integerMode, .Machine$integer.max, Inf), .) %>% diff()
    curvature = extrema %>% c(ifelse(integerMode, .Machine$integer.max, Inf)) %>% #for curvature, append biggest value (first and last sample point have no defined curvature)
      diff()
    result = curvature %>% #{. == 0} %>% which() #not working for discrete graphs
      {. < 0} %>% #which sample points are curved downwards (compared to adjacent)?
      rle() #run length of consecutive curvature
    
    if (upOnly) keep = result$values
    
    result = result %>% .$lengths %>% #for how long are segments curved upwards consecutively?
      cumsum() #cumulative sum to find indices of inflection points
      
    if (upOnly) result = result[keep]
    
    return(result)
  }
  
  create.SCL.plots = function(files, saveToFile=T, ecgPlot=F) {
    subdirEda = "plots eda/"; dir.create(paste0(path.rds, subdirEda))
    #subdirEdaMiddlePlots = "plots eda/middle plots/"; dir.create(paste0(path.rds, subdirEdaMiddlePlots))
    if (ecgPlot) {subdirEcg = "plots ecg/"; dir.create(paste0(path.rds, subdirEcg))}
    for (file in files) {
      cat(file, " ... ")
      
      filename = file %>% sub("\\..*", "", .) #get rid of file extension
      eda = read.phys(paste0(path.phys, file))
      triggers = get.trigger.onsets(eda$Trigger)
      if (length(triggers) < trials.n) { print(paste0("Skipping file due to too few trials: ", file)); next }
      triggers.first = min(triggers %>% tail(trials.n)) #if physio was recorded during timing training, triggers will be in there
      triggers.last = max(triggers) + trial.duration
      edaCut = eda[triggers.first:triggers.last,]
      
      #EDA plot cut to experiment range
      with(edaCut, plot.ts(EDA)); title(filename)
      if (saveToFile) { dev.copy(png, file=paste0(path.rds, subdirEda, filename, ".png"), width=800, height=600); dev.off() }
      
      #full EDA plot
      # with(eda, plot.ts(EDA)); title(filename)
      # abline(v=c(triggers.first, triggers.last), lty=3, col="red") #mark start and end of experiment trials
      # if (saveToFile) { dev.copy(png, file=paste0(path.rds, "plots eda/", filename, ".png"), width=800, height=600); dev.off() }
      
      #middle EDA plot
      # middle = mean(c(triggers.first, triggers.last))
      # subplot = (middle-eda.middleplot.time/2):(middle+eda.middleplot.time/2)
      # with(eda, plot.ts(EDA[subplot])); title(filename)
      # abline(v=triggers-middle-eda.middleplot.time/2, lty=3) #mark trigger onsets
      # if (saveToFile) { dev.copy(png, file=paste0(path.rds, subdirEdaMiddlePlots, filename, " 1min.png"), width=800, height=600); dev.off() }
      
      
      #middle ECG plot (full plot not useful)
      if (ecgPlot) {
        subplot = (middle-ecg.middleplot.time/2):(middle+ecg.middleplot.time/2)
        with(eda, plot.ts(ECG[subplot])); title(filename)
        abline(v=triggers-middle-ecg.middleplot.time/2, lty=3) #mark trigger onsets
        if (saveToFile) { dev.copy(png, file=paste0(path.rds, subdirEcg, filename, " 10sec.png"), width=800, height=600); dev.off() }
      }
    }
  }
}

# Preprocessing -----------------------------------------------------------
if (createSclPlots) create.SCL.plots(files.phys, saveToFile=saveSclPlotsToFile)

#read and tidy data
files.phys.included = files.phys %>% setdiff(exclusions.eda %>% toPhysFileName())

edas.list = list()
edas.df.list = list()
eda.maxima.list = list()
eda.minima.list = list()
eda.inflection.list = list()
eda.ucr = tibble(subject = files.phys.included %>% sub("\\..*", "", .), 
                 ucr = 0, valid = 0, lat = 0, rise = 0)
for (file in files.phys.included) {
  #file = files.phys.included %>% sample(1) #for testing
  cat("... ", file, "\n", sep="")
  
  filename = file %>% pathToCode()
  code = filename %>% codeToNum()
  eda = read.phys(paste0(path.phys, file)) %>% #select(-2) %>% 
    mutate(Trigger = Trigger %>% recode.triggers())
  
  if (exclusions.phys.trials[[filename]] %>% is.null() == F) {
    toExclude = exclusions.phys.trials[[filename]]
    triggers = eda$Trigger[eda$Trigger != 0]
    triggers[toExclude] = 0
    eda$Trigger[eda$Trigger!=0] = triggers
  }
  
  triggers.n = sum(eda$Trigger != 0)
  
  ratingfile = ratings %>% filter(subject==code) #TODO what if excluded from ratings?
  
  if (triggers.n < trials.n) { warning(paste0("Warning: Too few triggers in ", file, ". Using tail of conditions.")) }
  if (triggers.n > trials.n) { warning(paste0("Warning: Too many triggers in ", file, ". Using tail of triggers.")) }
  
  shocks = ratingfile$shock %>% tail(n=triggers.n) #idea: if trials are missing, assume the first ones are => if 158 of 160 trials, take all but first two
  conditions = ratingfile %>% unite("condition", pair, threat_num, sep="") %>% .$condition %>% as.integer() %>% 
    tail(n=triggers.n) %>% #idea: if trials are missing, assume the first ones are => if 158 of 160 trials, take all but first two
    c(rep(0, times=max(c(triggers.n - length(.), 0))), .) #prepend 0s until length of triggers.n is reached
  
  
  eda = eda %>% mutate(sample = 1:n(), time = 1:nrow(eda) / sample.rate - (1 / sample.rate),
                       Trigger = conditions %>% replace(Trigger, Trigger != 0, .) %>% as.numeric()) #inject conditions as triggers of onsets
  #if (checkPlots) {eda %>% ggplot(aes(x=time, y=EDA)) + geom_line() + geom_vline(xintercept = eda$time[eda$Trigger!=0]) + myGgTheme} %>% print()
  
  if (downsampling) {
    conversion = round(sample.rate / sample.rate.new)
    eda.down = data.frame(time = sample.down(eda$time, conversion),
                          EDA = sample.down(eda$EDA, conversion), 
                          Trigger = 0) %>% 
      mutate(sample=1:n()) %>% select(sample, everything())
    
    #calculate closest position for trigger onsets in downsampled time
    triggers.time.old = eda$time[eda$Trigger != 0]
    triggers.indices.new = triggers.time.old %>% closestRepresentative(eda.down$time, returnIndices = T) #for each old trigger time, find index of closest existing downsampled time
    eda.down$Trigger[triggers.indices.new] = eda$Trigger[eda$Trigger != 0] #inject conditions as triggers of onsets
    
    eda = eda.down; rm(eda.down)
    if (checkPlots) {eda %>% ggplot(aes(x=time, y=EDA)) + geom_line() + geom_vline(xintercept = eda$time[eda$Trigger!=0]) + myGgTheme} %>% print()
  }
  
  if (smoothing) {
    reps = 100 #filters can produce artifacts on the edges => fill edges with first and last values
    eda.filt = c(rep(first(eda$EDA), reps), eda$EDA, rep(last(eda$EDA), reps)) %>% 
      signal::filtfilt(signal::butter(round(12/6), lowPassFreq/(ifelse(downsampling, sample.rate.new, sample.rate)/2)), .) #low-pass filter
    eda$EDA = eda.filt[(reps+1):(length(eda.filt)-reps)]; rm(eda.filt)
    
    if (checkPlots) {
      eda %>% ggplot(aes(x=time, y=EDA)) + geom_line() + geom_vline(xintercept = eda$time[eda$Trigger!=0]) + myGgTheme} %>% print()
  }
  
  if (saveLedalab) {
    output = eda
    
    #trim range according to triggers +/- analysis.border
    triggers.index = which(output$Trigger!=0)
    output = output[max(0, first(triggers.index)-analysis.border):min(nrow(output), last(triggers.index)+analysis.border),]
    
    dir.create(paste0(path.phys, "downsampled/"), showWarnings = F)
    write.table(output, paste0(path.phys, "downsampled/", file),row.names=F, col.names=F, quote=F)
    rm(output)
  }
  
  #peek at result
  #eda %>% filter(Trigger != 0) %>% glimpse()
  #eda %>% filter(Trigger != 0) %>% View()
  #eda %>% ggplot(aes(x=time, y=EDA)) + geom_line() + geom_vline(xintercept = eda$time[eda$Trigger!=0]) + myGgTheme
  
  
  #one overview table for each subject with one row per trial
  eda.vp = eda %>% filter(Trigger != 0) %>% select(-EDA) %>% 
    mutate(trial = 1:n()+(trials.n-triggers.n), condition = Trigger, threat = condition %% 10, pair = condition %/% 10,
           phase = case_when(trial < preAcqEnd ~ "Hab",
                             trial <    acqEnd ~ "Acq",
                             TRUE ~              "Gen"),
           shock = shocks, shockPrior = c(F, lag(shocks)[-1]),
           time.start = time, time.end = time + min(itiEnd)/1000,
           sample.start = sample, 
           sample.end = {time.end * ifelse(downsampling, sample.rate.new, sample.rate)} %>% round()
    ) %>% select(-Trigger, -time, -sample) %>% select(trial, condition, everything())
  
  #nest EDA time series data into overview table
  # eda.vp$EDA = eda.vp$sample.start %>% lapply(function(start) 
  #   eda %>% filter(sample >= start, 
  #                  sample < start + min(itiEnd)*2/1000*ifelse(downsampling, sample.rate.new, sample.rate)))
  #eda.vp$EDA[[1]] %>% glimpse() #access EDA time series data for trial 1
  
  edas.df.list[[filename]] = eda.vp
  edas.list[[filename]] = eda
  
  # Quantify URs ------------------------------------------------------------
  eda.vp.ucr = eda.vp %>% filter(shock == T) %>% select(-contains("end")) %>% 
    mutate(index = {time.start + shockTime/1000} %>% closestRepresentative(eda$time, returnIndices=T),
           time.start = eda$time[index], sample.start = eda$sample[index])

  eda.maxima = data.frame(sample = localMaxima(eda$EDA)) %>% 
    mutate(n = 1:n(), time=eda$time[sample], EDA = eda$EDA[sample]) %>% select(n, everything())
  eda.minima = data.frame(sample = localMinima(eda$EDA)) %>% 
    mutate(n = 1:n(), time=eda$time[sample], EDA = eda$EDA[sample]) %>% select(n, everything())
  eda.inflection = data.frame(sample = inflectionPoints(eda$EDA)) %>% 
    mutate(n = 1:n(), time=eda$time[sample], EDA = eda$EDA[sample]) %>% select(n, everything())
  
  eda.maxima.list[[filename]] = eda.maxima
  eda.minima.list[[filename]] = eda.minima
  eda.inflection.list[[filename]] = eda.inflection
  
  if (checkPlots || saveSclPlotsToFile) {
    sclPlot = eda %>% ggplot(aes(x=time, y=EDA)) + 
      geom_vline(xintercept = eda$time[eda$Trigger!=0], color=ifelse(eda.vp$shock, "orange", "black")) +
      geom_line() + 
      #geom_point(data=eda.inflection, color="purple") + 
      geom_point(data=eda.maxima, color="red") + 
      geom_point(data=eda.minima, color="blue") + 
      ggtitle(filename) + myGgTheme
    if (checkPlots) sclPlot %>% print()
    if (saveSclPlotsToFile) sclPlot %>% ggsave(paste0(path.rds, "plots eda/", filename, ".png"), plot=., device="png", width=1920/300, height=1080/300, dpi=300)
  }
  
  for (t in 1:nrow(eda.vp.ucr)) {
    start = eda.vp.ucr$time.start[t]
    eda.minima.trial = eda.minima %>% filter(time >= start + min(ucrMinWindow),
                                             time <= start + max(ucrMinWindow))
    eda.maxima.trial = eda.maxima %>% filter(time > suppressWarnings(min(eda.minima.trial$time))) %>% 
      head(nrow(eda.minima.trial)) #same amount of maxima as minima
    
    #old (obsolete by adding inflection points; see below)
    #check if segment onset needs to be added as minimum (in case of overlapping SCRs without local minimum)
    # minTooLate = eda.minima %>% filter(time >= start + max(ucrMinWindow)) %>% head(1) %>% .$time #take first true minimum that is too late
    # if (is_empty(minTooLate) || minTooLate > eda.maxima.trial %>% head(1) %>% .$time) #if first minimum that is too late occurs after first maximum in segment (i.e., true minimum is before segment onset => overlapping SCRs without minimum)
    #   eda.minima.trial = eda.minima.trial %>% bind_rows(eda %>% filter(time == start) %>% select(-Trigger) %>% mutate(n = 0), .)
    
    #if no minimum found within time range use deflection point as minimum
    if (nrow(eda.minima.trial)==0) { 
      warning(paste0("No minimum found. Using inflection point for subject ", filename, ", trial ", t))
      eda.minima.trial = eda.inflection %>% filter(time >= start + min(ucrMinWindow),
                                                   time <= start + max(ucrMinWindow))
      eda.maxima.trial = eda.maxima %>% filter(time > suppressWarnings(min(eda.minima.trial$time))) %>% 
        head(nrow(eda.minima.trial)) #same amount of maxima as minima
      
      #check if mins and maxs alternate (in case inflection points were used)
      while (nrow(eda.maxima.trial)>1) {
        latestMaxTimes = eda.maxima.trial %>% tail(2) %>% .$time
        if (eda.minima.trial %>% filter(time < max(latestMaxTimes), #check if any inflection point is between last two maxima (may be false if two inflection points before earlier max)
                                        time > min(latestMaxTimes)) %>% 
            nrow() %>% {. > 0}) break
        else { #if no inflection point between last two maxima
          eda.minima.trial = eda.minima.trial %>% head(-1) #delete last inflection point
          eda.maxima.trial = eda.maxima.trial %>% head(nrow(eda.minima.trial)) #same amount of maxima as minima
        }
      }
      
      #check if true min occurs between inflection point and corresponding max (with true min being too late)
      toKeep = c()
      if (nrow(eda.minima.trial) > 0) { 
        for (i in 1:nrow(eda.minima.trial)) {
          if (eda.minima %>% filter(time > eda.minima.trial[i, "time"], time < eda.maxima.trial[i, "time"]) %>% 
              nrow() %>% {. == 0}) toKeep = toKeep %>% c(i)
        }
        eda.minima.trial = eda.minima.trial[toKeep,]; eda.maxima.trial = eda.maxima.trial[toKeep,] #only keep inflection points that don't have a true minimum between themself and their max
      }
    }
    
    #calculate SCRs
    eda.minima.trial.renamed = eda.minima.trial; names(eda.minima.trial.renamed) = names(eda.minima.trial.renamed) %>% paste0(".min")
    eda.maxima.trial.renamed = eda.maxima.trial; names(eda.maxima.trial.renamed) = names(eda.maxima.trial.renamed) %>% paste0(".max")
    scr = bind_cols(eda.minima.trial.renamed, eda.maxima.trial.renamed)
    #names(scr) = names(scr) %>% {ifelse(endsWith(., "1"), gsub("1", ".max", ., fixed=T), paste(., "min", sep="."))}
    scr = scr %>% mutate(SCR = EDA.max - EDA.min,
                         Lat = time.min - start,
                         Rise = time.max - start) %>% select(SCR, Lat, Rise, everything()) %>% 
      filter(SCR == max(SCR)) #take maximum SCR
    
    if (nrow(scr) == 0) {scr = scr %>% bind_rows(tibble(SCR = 0, Lat = NA, Rise = NA))
    } else if (scr$SCR <= 0) {scr = scr %>% mutate(SCR = 0, Lat = NA, Rise = NA)}
    
    scr.out = scr %>% select(1:3)
    eda.vp.ucr[t, names(scr.out)] = scr.out
    
    if (checkPlotsTrial) {
      {eda %>% filter(time >= eda.vp.ucr$time.start[t], time <= max(itiEnd)/1000 - shockTime/1000 + eda.vp.ucr$time.start[t]) %>% 
              ggplot(aes(x=time, y=EDA)) + 
              geom_vline(xintercept = ucrMinWindow + eda.vp.ucr$time.start[t], linetype="dashed", color="blue") + 
              #geom_vline(xintercept = itiEnd/1000 + eda.vp.ucr$time.start[t]) +
              geom_line() + #EDA
              geom_point(data=bind_rows(eda.minima.trial %>% filter(time == scr$time.min),
                                        eda.maxima.trial %>% filter(time == scr$time.max)),
                         pch=21, size=3) + #selected extreme values first as highlighting circles
              geom_point(data=eda.maxima.trial, color="red") + #all maxima
              geom_point(data=eda.minima.trial, color="blue") + #all minima
              myGgTheme + ggtitle(paste0(filename, ": UR ", t, "/", nrow(eda.vp.ucr))) + xlab("Time (s)")} %>% 
        ggsave(paste0("plots eda/UCS/trials/", filename, " UR ", t, ".png"), ., width=1920, height=1080, units="px")
        
      #   print()
      # invisible(readline(prompt="Press [enter] to continue"))
    }
  }
  
  # #old: just as max - min (but max could be before min and other problems)
  # eda.vp.ucr = eda.vp %>% filter(shock == T) %>%
  #   mutate(UCR = EDA %>% sapply(function(df) {
  #     temp = df %>% mutate(sample = sample - min(sample) + 1, time = time - min(time)) %>%  #relative to trial onset
  #       filter(time > shockEnd/1000, time < shockEnd/1000 + latestPeak)
  #     return(temp$EDA %>% {max(.) - min(.)})
  #   }),
  #   valid = UCR %>% unlist() %>% {. >= ucrMin},
  #   UCR = UCR %>% unlist() %>% ifelse(. < ucrMin, 0, .))
  # 
  # #test = eda.vp.ucr %>% {sapply(.$EDA, function(df, i) df %>% mutate(trial = i), i = .$trial)}
  
  eda.ucr[eda.ucr$subject==filename, -1] = eda.vp.ucr %>% mutate(valid = SCR > ucrMin) %>% 
    summarise(scr = mean(SCR, na.rm=T), 
              valid = mean(valid),
              lat = mean(Lat, na.rm=T),
              rise = mean(Rise, na.rm=T)
    )
  
  if (ucrPlots) {
    path = paste0(path.rds, "plots eda/UCS/")
    dir.create(path, showWarnings=F)
    
    unified = eda.vp.ucr %>% mutate(EDA = sample.start %>% lapply(function(start) 
      eda %>% filter(sample >= start,
                     sample <= start + max(itiEnd)/1000*ifelse(downsampling, sample.rate.new, sample.rate)) %>% select(-Trigger))
    )
    
    for (t in 1:nrow(unified)) {
      unified$EDA[[t]] = unified$EDA[[t]] %>% 
        mutate(trial = t, sample = sample - min(sample) + 1, time = time - min(time)) #unify starting time to allow overlap
    }
    
    path.ucs = paste0(path.rds, "plots eda/UCS/")
    if (dir.exists(path.ucs)==F) dir.create(path.ucs, recursive=T)
    unified$EDA %>% bind_rows() %>% 
      mutate(trial = as.factor(trial)) %>% #filter(time > shockEnd/1000) %>% 
    {ggplot(., aes(x=time, y=EDA, color=trial)) + 
        #geom_rect(xmin=0, xmax=max(itiEnd)/1000, ymin=-Inf, ymax=Inf, color="grey", alpha=.1) +
        #geom_vline(xintercept=c(0, max(itiEnd)/1000), color="grey") + #borders of min/max scoring
        geom_vline(xintercept = ucrMinWindow, linetype="dashed", color="blue") + 
        geom_vline(xintercept=(itiEnd-shockTime)/1000) + #borders of ITI
        geom_path() + scale_color_viridis_d() +
        ggtitle(filename) + myGgTheme + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))} %>% 
      #print()
      ggsave(paste0(path.ucs, filename, ".png"), plot=., device="png", width=1920/300, height=1080/300, dpi=300)
  }
  
  if (checkPlots && !checkPlotsTrial) invisible(readline(prompt="Press [enter] to continue"))
}
rm(eda, eda.vp)

eda.ucr = eda.ucr %>% mutate(subject = subject %>% gsub("vp", "", ., fixed=T) %>% as.integer(),
                             include = ucr > ucrMeanAmp & valid > ucrMeanValid, 
                             ln_ucr = log(ucr + 1))
#all(eda.ucr == read_rds("eda.ucr.rds" %>% paste0(path.rds, .)), na.rm=T) #check equivalence of processing
#eda.ucr %>% write_rds("eda.ucr.rds" %>% paste0(path.rds, .))


# Exclusions: UCS non-responder -------------------------------------------
#eda.ucr = read_rds("eda.ucr.rds" %>% paste0(path.rds, .))
#eda.ucr = eda.ucr %>% mutate(include = ifelse(subject %in% c(4, 28), F, include)) #discard manually according to UCS plot

sum(eda.ucr$include == F) %>% paste0(" subjects excluded (", round(100 - mean(eda.ucr$include)*100, digits = 2), "%)")
#eda.ucr %>% filter(include==FALSE) %>% arrange(desc(valid))

with(eda.ucr, hist(valid, breaks=seq(0, 1, length.out=20+1), main="Valid URs")); abline(v=ucrMeanValid, col="red", lwd=3, lty=2)
with(eda.ucr, hist(ucr, breaks=seq(0, ceiling(max(ucr)), by=ucrMeanAmp), main="Size of URs")); abline(v=ucrMeanAmp, col="red", lwd=3, lty=2)

#too few reactions
#sum(eda.ucr$valid < ucrMeanValid) / nrow(eda.ucr)
#eda.ucr[eda.ucr$valid < ucrMeanValid,]

#too little average reaction
#sum(eda.ucr$ucr < ucrMeanAmp) / nrow(eda.ucr)
#eda.ucr[eda.ucr$ucr < ucrMeanAmp,]


# Quantify CRs ------------------------------------------------------------
#eda.ucr = read_rds("eda.ucr.rds" %>% paste0(path.rds, .))
subjects = eda.ucr %>% .$subject #no exclusions by UCR
#subjects = eda.ucr %>% filter(include) %>% .$subject #with exclusions by UCR
for (s in seq(subjects)) {
  subject = subjects[s]
  cat("... ", subject, "\n", sep="")
  
  eda.vp = edas.df.list[[subject %>% toCode()]]
  eda = edas.list[[subject %>% toCode()]]
  eda.maxima = eda.maxima.list[[subject %>% toCode()]]
  eda.minima = eda.minima.list[[subject %>% toCode()]]
  eda.inflection = eda.inflection.list[[subject %>% toCode()]]
  
  for (t in 1:nrow(eda.vp)) {
    trial = eda.vp$trial[t]
    start = eda.vp$time.start[t]
    eda.minima.trial = eda.minima %>% filter(time >= start + min(crMinWindow),
                                             time <= start + max(crMinWindow))
    eda.maxima.trial = eda.maxima %>% filter(time > suppressWarnings(min(eda.minima.trial$time))) %>% 
      head(nrow(eda.minima.trial)) #same amount of maxima as minima
    if (nrow(eda.maxima.trial) < nrow(eda.minima.trial)) { #if recording ended before next max
      eda.maxima.trial = eda.maxima.trial %>% bind_rows(
        eda %>% filter(time >= max(eda.minima.trial$time)) %>% filter(EDA==max(EDA)) %>% mutate(n = 0) %>% select(n, everything()) %>% select(-Trigger)
      )
    }
    
    #old (obsolete by adding inflection points; see below)
    #check if segment onset needs to be added as minimum (in case of overlapping SCRs without local minimum)
    # minTooLate = eda.minima %>% filter(time >= start + max(crMinWindow)) %>% head(1) %>% .$time #take first true minimum that is too late
    # if (is_empty(minTooLate) || minTooLate > eda.maxima.trial %>% head(1) %>% .$time) #if first minimum that is too late occurs after first maximum in segment (i.e., true minimum is before segment onset => overlapping SCRs without minimum)
    #   eda.minima.trial = eda.minima.trial %>% bind_rows(eda %>% filter(time == start) %>% select(-Trigger) %>% mutate(n = 0), .)
    
    #if no minimum found within time range use deflection point as minimum
    if (nrow(eda.minima.trial)==0) { 
      warning(paste0("No minimum found. Using inflection point for subject ", filename, ", trial ", trial))
      eda.minima.trial = eda.inflection %>% filter(time >= start + min(crMinWindow),
                                                   time <= start + max(crMinWindow))
      eda.maxima.trial = eda.maxima %>% filter(time > suppressWarnings(min(eda.minima.trial$time))) %>% 
        head(nrow(eda.minima.trial)) #same amount of maxima as minima
      if (nrow(eda.maxima.trial) < nrow(eda.minima.trial)) { #if recording ended before next max
        eda.maxima.trial = eda.maxima.trial %>% bind_rows(
          eda %>% filter(time >= max(eda.minima.trial$time)) %>% filter(EDA==max(EDA)) %>% mutate(n = 0) %>% select(n, everything()) %>% select(-Trigger)
        )
      }
      
      #check if mins and maxs alternate (in case inflection points were used)
      while (nrow(eda.maxima.trial)>1) {
        latestMaxTimes = eda.maxima.trial %>% tail(2) %>% .$time
        if (eda.minima.trial %>% filter(time < max(latestMaxTimes), #check if any inflection point is between last two maxima (may be false if two inflection points before earlier max)
                                        time > min(latestMaxTimes)) %>% 
            nrow() %>% {. > 0}) break
        else { #if no inflection point between last two maxima
          eda.minima.trial = eda.minima.trial %>% head(-1) #delete last inflection point
          eda.maxima.trial = eda.maxima.trial %>% head(nrow(eda.minima.trial)) #same amount of maxima as minima
        }
      }
      
      #check if true min occurs between inflection point and corresponding max (with true min being too late)
      toKeep = c()
      if (nrow(eda.minima.trial) > 0) { 
        for (i in 1:nrow(eda.minima.trial)) {
          if (eda.minima %>% filter(time > eda.minima.trial[i, "time"], time < eda.maxima.trial[i, "time"]) %>% 
              nrow() %>% {. == 0}) toKeep = toKeep %>% c(i)
        }
        eda.minima.trial = eda.minima.trial[toKeep,]; eda.maxima.trial = eda.maxima.trial[toKeep,] #only keep inflection points that don't have a true minimum between themself and their max
      }
    }
    
    #calculate SCRs
    eda.minima.trial.renamed = eda.minima.trial; names(eda.minima.trial.renamed) = names(eda.minima.trial.renamed) %>% paste0(".min")
    eda.maxima.trial.renamed = eda.maxima.trial; names(eda.maxima.trial.renamed) = names(eda.maxima.trial.renamed) %>% paste0(".max")
    scr = bind_cols(eda.minima.trial.renamed, eda.maxima.trial.renamed)
    #names(scr) = names(scr) %>% {ifelse(endsWith(., "1"), gsub("1", ".max", ., fixed=T), paste(., "min", sep="."))}
    scr = scr %>% mutate(SCR = EDA.max - EDA.min,
                         Lat = time.min - start,
                         Rise = time.max - start) %>% select(SCR, Lat, Rise, everything()) %>% 
      filter(SCR == max(SCR)) #take maximum SCR
    
    if (nrow(scr) == 0) {scr = scr %>% bind_rows(tibble(SCR = 0, Lat = NA, Rise = NA))
    } else if (scr$SCR <= 0) {scr = scr %>% mutate(SCR = 0, Lat = NA, Rise = NA)}
    
    scr.out = scr %>% select(1:3)
    eda.vp[t, names(scr.out)] = scr.out
    
    if (checkPlotsTrial && (checkPlotsTrialGenOnly==F || (checkPlotsTrialGenOnly==T && eda.vp$phase[t]=="Gen"))) {
      {eda %>% filter(time >= eda.vp$time.start[t], time <= max(itiEnd)/1000 - shockTime/1000 + eda.vp$time.start[t]) %>% 
              ggplot(aes(x=time, y=EDA)) + 
              geom_vline(xintercept = crMinWindow + eda.vp$time.start[t], linetype="dashed", color="blue") + 
              #geom_vline(xintercept = itiEnd/1000 + eda.vp$time.start[t]) +
              geom_line() + #EDA
              geom_point(data=bind_rows(eda.minima.trial %>% filter(time == scr$time.min),
                                        eda.maxima.trial %>% filter(time == scr$time.max)),
                         pch=21, size=3) + #selected extreme values first as highlighting circles
              geom_point(data=eda.maxima.trial, color="red") + #all maxima
              geom_point(data=eda.minima.trial, color="blue") + #all minima
                myGgTheme + ggtitle(paste0(subject %>% toCode(), ": CR ", trial, "/", trials.n)) + xlab("Time (s)")} %>% 
      ggsave(paste0("plots eda/CS/trials/", subject %>% toCode(), " CR ", trial, ".png"), ., width=1920, height=1080, units="px")
    
    #   print()
    # invisible(readline(prompt="Press [enter] to continue"))
    }
  }
  
  # #old: just as max - min (but max could be before min and other problems)
  # eda.vp = eda.vp %>% mutate(
  #   CR = EDA %>% sapply(function(df) {
  #     temp = df %>% mutate(sample = sample - min(sample) + 1, time = time - min(time)) %>%  #relative to trial onset
  #       filter(time < min(itiEnd))
  #     return(temp$EDA %>% {max(.) - min(.)})
  #   }),
  #   #valid = CR %>% {. >= crMin}, #not needed! invalid CRs are coded as 0 and used for mean calculation!
  #   CR = CR %>% ifelse(. < crMin, 0, .),
  #   subject = subject)
  
  edas.df.list[[subject]] = eda.vp
  
  if (crPlots) {
    unified = eda.vp %>% mutate(EDA = sample.start %>% lapply(function(start) 
      eda %>% filter(sample >= start,
                     sample <= start + max(itiEnd)/1000*ifelse(downsampling, sample.rate.new, sample.rate)) %>% select(-Trigger))
    ) %>% filter(phase == "Gen")
    
    for (t in 1:nrow(unified)) {
      unified$EDA[[t]] = unified$EDA[[t]] %>% 
        mutate(trial = t, threat = eda.vp$threat[t], sample = sample - min(sample) + 1, time = time - min(time)) #unify starting time to allow overlap
    }
    
    path.cs = paste0(path.rds, "plots eda/CS/")
    if (dir.exists(path.cs)==F) dir.create(path.cs, recursive=T)
    unified$EDA %>% bind_rows() %>% 
      mutate(trial = as.factor(trial)) %>% #filter(time > shockEnd/1000) %>% 
      {ggplot(., aes(x=time, y=EDA, color=trial)) + facet_wrap(vars(threat)) + 
          #geom_rect(xmin=0, xmax=max(itiEnd)/1000, ymin=-Inf, ymax=Inf, color="grey", alpha=.1) +
          #geom_vline(xintercept=c(0, max(itiEnd)/1000), color="grey") + #borders of min/max scoring
          geom_vline(xintercept = crMinWindow, linetype="dashed", color="blue") + 
          geom_vline(xintercept=(itiEnd-shockTime)/1000) + #borders of ITI
          geom_path() + scale_color_viridis_d() +
          ggtitle(subject %>% toCode()) + myGgTheme + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))} %>% 
      #print()
      ggsave(paste0(path.cs, subject %>% toCode(), ".png"), plot=., device="png", width=1920/300, height=1080/300, dpi=300)
  }
  
  if (checkPlots && !checkPlotsTrial) invisible(readline(prompt="Press [enter] to continue"))
}
rm(eda, eda.vp)

#eda dataframe without time series data
# eda.df = edas.df.list[[subjects[1]]] %>% select(-EDA)
# for (i in 2:length(subjects)) eda.df = edas.df.list[[subjects[i]]] %>% select(-EDA) %>% bind_rows(eda.df, .)
eda.df = edas.df.list %>% bind_rows(.id="subject") %>% 
  mutate(diagnostic = {pair %in% c(1, 3)} %>% ifelse("Eyes", "Mouth/Nose") %>% as.factor(),
         pairs = {pair %in% 2:3} %>% ifelse(2, 1) %>% as.factor(),
         ln_cr = log(SCR + 1)) %>% select(subject, everything())
rm(edas.df.list, edas.list)
#all(eda.df == read_rds("eda.rds" %>% paste0(path.rds, .)), na.rm=T) #check equivalence of processing
#eda.df %>% write_rds("eda.rds" %>% paste0(path.rds, .))


# Inference Tests ---------------------------------------------------------
#eda.ucr = read_rds("eda.ucr.rds" %>% paste0(path.rds, .))
#eda.df = read_rds("eda.rds" %>% paste0(path.rds, .))
eda.df = eda.df %>% filter(subject %in% {eda.ucr %>% filter(include) %>% .$subject}) #only analyze US responders

#habituation
eda.df.hab = eda.df %>% filter(phase == "Hab") %>% group_by(subject, threat) %>% 
  summarise(ln_cr = mean(ln_cr, na.rm=T)) %>% mutate(threat = threat %>% recode(`1` = "CS-", `6` = "CS+"))
# d = .5 #effect size that constitutes equivalence bounds
# eda.df.hab %>% spread(threat, ln_cr) %>% 
#   TOSTER::dataTOSTpaired(pairs=list(c(i1="CS+", i2="CS-")), low_eqbound=-d, high_eqbound=d, plots=F, desc=F) #for Habituation, equivalence hypothesis => TOST test (see https://rpsychologist.com/d3/equivalence/)
eda.df.hab %>% spread(threat, ln_cr) %>% 
  t.test(x=.$"CS-", y=.$"CS+", alternative="two.sided", mu=0, paired=T) %>% apa::t_apa(es_ci=T)
eda.df.hab %>% group_by(threat) %>% summarise(ln_cr.se = se(ln_cr, na.rm=T), ln_cr = mean(ln_cr, na.rm=T)) %>%
  ggplot(aes(x=threat, y=ln_cr, group=NA)) +
  geom_point(size=5) + #geom_path(size=2) +
  geom_errorbar(aes(ymin=ln_cr-ln_cr.se*1.96, ymax=ln_cr+ln_cr.se*1.96), size=2) +
  geom_line(data=eda.df.hab, mapping=aes(x=threat, y=ln_cr, group=subject)) + myGgTheme

#acquisition

#TODO apply exclusions.onlyGen to habituation & acquisition

eda.df.acq = eda.df %>% filter(phase == "Acq", shock == F, shockPrior == F) %>% 
  group_by(subject, threat) %>% summarise(ln_cr = mean(ln_cr, na.rm=T)) %>% mutate(threat = threat %>% recode(`1` = "CS-", `6` = "CS+"))
eda.df.acq %>% spread(threat, ln_cr) %>% 
  t.test(x=.$"CS-", y=.$"CS+", alternative="less", mu=0, paired=T) %>% apa::t_apa(es_ci=T)
eda.df.acq %>% group_by(threat) %>% summarise(ln_cr.se = se(ln_cr, na.rm=T), ln_cr = mean(ln_cr, na.rm=T)) %>% 
  ggplot(aes(x=threat, y=ln_cr, group=NA)) + 
  geom_point(size=5) + #geom_path(size=2) +
  geom_errorbar(aes(ymin=ln_cr-ln_cr.se*1.96, ymax=ln_cr+ln_cr.se*1.96), size=2) +
  geom_line(data=eda.df.acq, mapping=aes(x=threat, y=ln_cr, group=subject)) + myGgTheme

#generalization
eda.df.gen = eda.df %>% filter(phase == "Gen", shock == F, shockPrior == F) %>% 
  group_by(subject, threat, diagnostic, pairs) %>% summarise(ln_cr = mean(ln_cr, na.rm=T)) %>% 
  ungroup() %>% mutate(threat = as.factor(threat)) %>% #don't recode yet! messes up the ggplot
  #mutate(threat = threat %>% recode(`1` = "CS-", `2` = "GS1", `3` = "GS2", `4` = "GS3", `5` = "GS4", `6` = "CS+"))
merge(questionnaires, ., by="subject") %>% tibble()

eda.df.gen %>% #mutate(SPAI = scale(SPAI), STAI = scale(STAI)) %>% #no effect - done implicitly?
  ez::ezANOVA(dv=.(ln_cr), wid=.(subject), 
              within=.(threat, diagnostic), 
              #between=.(pairs),
              between=.(SPAI), observed=SPAI,
              #between=.(STAI), observed=STAI,
              detailed=T, type=2) %>% apa::anova_apa(force_sph_corr=T)

#SPAI main effect
# eda.df.gen %>% group_by(subject) %>% summarise(ln_cr.se = se(ln_cr, na.rm=T), ln_cr = mean(ln_cr, na.rm=T)) %>% merge(questionnaires, by="subject") %>% 
#   with(cor.test(ln_cr, SPAI)) %>% correlation_out()
# eda.df.gen %>% group_by(subject) %>% summarise(ln_cr.se = se(ln_cr, na.rm=T), ln_cr = mean(ln_cr, na.rm=T)) %>% merge(questionnaires, by="subject") %>% 
#   ggplot(aes(x=SPAI, y=ln_cr, color=SPAI, fill=SPAI)) +
#   geom_errorbar(aes(ymin=ln_cr-ln_cr.se*1.96, ymax=ln_cr+ln_cr.se*1.96), width=spai.width) +
#   stat_smooth(method="lm", color = "black") +
#   #geom_point(size=4, shape=21, color="black") +
#   geom_point(size=4) + 
#   ylab("ln(1 + SCR)") +
#   scale_color_viridis_c() + myGgTheme + theme(legend.position = "none")

#STAI main effect
# eda.df.gen %>% group_by(subject) %>% summarise(ln_cr.se = se(ln_cr, na.rm=T), ln_cr = mean(ln_cr, na.rm=T)) %>% merge(questionnaires, by="subject") %>% 
#   with(cor.test(ln_cr, STAI)) %>% correlation_out()
# eda.df.gen %>% group_by(subject) %>% summarise(ln_cr.se = se(ln_cr, na.rm=T), ln_cr = mean(ln_cr, na.rm=T)) %>% merge(questionnaires, by="subject") %>% 
#   ggplot(aes(x=STAI, y=ln_cr, color=STAI, fill=STAI)) +
#   geom_errorbar(aes(ymin=ln_cr-ln_cr.se*1.96, ymax=ln_cr+ln_cr.se*1.96), width=1) +
#   stat_smooth(method="lm", color = "black") +
#   #geom_point(size=4, shape=21, color="black") +
#   geom_point(size=4) + 
#   ylab("ln(1 + SCR)") +
#   scale_color_viridis_c() + myGgTheme + theme(legend.position = "none")

#threat main effect
eda.df.gen.subj = eda.df.gen %>% group_by(threat, subject) %>% summarise(ln_cr.se = se(ln_cr, na.rm=T), ln_cr = mean(ln_cr, na.rm=T))
eda.df.gen %>% group_by(threat, subject) %>% summarise(ln_cr = mean(ln_cr, na.rm=T)) %>% summarise(ln_cr.se = se(ln_cr, na.rm=T), ln_cr = mean(ln_cr, na.rm=T)) %>% 
  ggplot(aes(x=threat, y=ln_cr, color=threat, group=NA)) + 
  geom_dotplot(data=eda.df.gen.subj, mapping=aes(group=threat, fill=threat), binaxis="y", alpha=.25, color="black", stackratio=1, stackdir="centerwhole", dotsize=.5) +
  #geom_path(data=eda.df.gen.ga.threat %>% filter(threat %in% c(1, 6)), color = "black", size=1.5) + #generalization line
  geom_errorbar(aes(ymin=ln_cr-ln_cr.se*1.96, ymax=ln_cr+ln_cr.se*1.96), size=1.5) +
  geom_line(size=1) + 
  geom_point(size=4.5) +
  #geom_line(data=eda.df.gen, mapping=aes(x=threat, y=ln_cr, group=subject), alpha=.1) + #individual gradients
  scale_color_manual(values=colors, guide=guide_legend(reverse=T)) +
  #scale_fill_manual(values=colors, guide=guide_legend(reverse=T)) +
  scale_fill_manual(values=rep("grey", 6), guide=guide_legend(reverse=T)) +
  scale_x_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) +
  ylab("LN(1 + SCR)") + xlab("Threat") + labs(color="Threat", fill="Threat") + myGgTheme +
  theme(legend.position = "none")

for (i in 2:6) { #CS- vs. rest
  levels = c(1, i)
  cat(paste0("\n\nComparing levels: ", paste(levels, collapse=" vs. "), "\n"))
  eda.df.gen.subj %>% filter(threat %in% levels) %>%
    t.test(ln_cr ~ threat, ., paired=T) %>% apa::t_apa(es_ci=T)
}

for (i in 1:5) { #CS+ vs. rest
  levels = c(i, 6)
  cat(paste0("\n\nComparing levels: ", paste(levels, collapse=" vs. "), "\n"))
  eda.df.gen.subj %>% filter(threat %in% levels) %>%
    t.test(ln_cr ~ threat, ., paired=T) %>% apa::t_apa(es_ci=T)
}

for (i in 1:6) { #GS3 vs. rest
  levels = c(i, 4)
  if (min(levels) == max(levels)) next
  cat(paste0("\n\nComparing levels: ", paste(levels, collapse=" vs. "), "\n"))
  eda.df.gen.subj %>% filter(threat %in% levels) %>%
    t.test(ln_cr ~ threat, ., paired=T) %>% apa::t_apa(es_ci=T)
}


#threat * diagnostic
eda.df.gen %>% group_by(threat, diagnostic, subject) %>% summarise(ln_cr = mean(ln_cr, na.rm=T)) %>% summarise(ln_cr.se = se(ln_cr, na.rm=T), ln_cr = mean(ln_cr, na.rm=T)) %>%
  ggplot(aes(x=threat, y=ln_cr, color=diagnostic, group=diagnostic)) +
  geom_point(size=3, position=dodge) + geom_line(size=1.5, , position=dodge) +
  geom_errorbar(aes(ymin=ln_cr-ln_cr.se*1.96, ymax=ln_cr+ln_cr.se*1.96), position=dodge, width=dodge.width) +
  #geom_line(data=eda.df.gen, mapping=aes(x=threat, y=ln_cr, group=subject), alpha=.1) + #individual gradients
  #scale_color_manual(values=colors, guide=guide_legend(reverse=T)) +
  scale_x_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) +
  ylab("LN(1 + SCR)") + xlab("Threat") + labs(color="Threat") + myGgTheme

# eda.df.gen %>% group_by(threat, diagnostic, pairs) %>% summarise(ln_cr.se = se(ln_cr, na.rm=T), ln_cr = mean(ln_cr, na.rm=T)) %>%
#   #ggplot(aes(x=threat, y=ln_cr, color=diagnostic, shape=pairs, group=interaction(diagnostic, pairs))) +
#   ggplot(aes(x=threat, y=ln_cr, color=diagnostic, group=diagnostic)) + facet_wrap(vars(pairs), labeller = "label_both") +
#   geom_point(size=3, position=dodge) + geom_line(size=1.5, , position=dodge) +
#   geom_errorbar(aes(ymin=ln_cr-ln_cr.se*1.96, ymax=ln_cr+ln_cr.se*1.96), position=dodge, width=dodge.width) +
#   #geom_line(data=eda.df.gen, mapping=aes(x=threat, y=ln_cr, group=subject), alpha=.1) + #individual gradients
#   #scale_color_manual(values=colors, guide=guide_legend(reverse=T)) +
#   scale_x_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) +
#   ylab("LN(1 + SCR)") + xlab("Threat") + labs(color="Threat") + myGgTheme

#SPAI * threat
eda.df.gen %>% group_by(SPAI, threat, subject) %>% summarise(ln_cr = mean(ln_cr, na.rm=T)) %>% summarise(ln_cr.se = se(ln_cr, na.rm=T), ln_cr = mean(ln_cr, na.rm=T)) %>%
  #ggplot(aes(x=SPAI, y=ln_cr, color=threat, fill=threat, group=threat)) +
  ggplot(aes(x=SPAI, y=ln_cr, color=threat, fill=threat, group=threat)) + facet_wrap(vars(threat), labeller="label_both") +
  geom_point(size=3) + 
  geom_smooth(method="lm", size=1.5, alpha = .2) +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors) + myGgTheme
#no positive slope for GS1 & GS3


# Gradient Analysis -------------------------------------------------------
eda.df.gen.simple = eda.df.gen %>% group_by(subject, threat, pairs) %>% summarise(ln_cr = mean(ln_cr, na.rm=T)) %>% mutate(diagnostic="all")
# eda.gradients.simple = tibble(subject=numeric(), lds=numeric(), diff=numeric(), level=numeric())
# for (s in eda.df.gen.simple$subject %>% unique()) {
#   eda.gradients.simple = eda.gradients.simple %>% bind_rows(
#     eda.df.gen.simple %>% filter(subject==s) %>% .$ln_cr %>% gradient.analysis() %>% #gradient analysis
#       tibble(subject=as.numeric(s), values=., names=names(.)) %>% spread(names, values)) #awkward stuff to append to tibble [c(s, .) converts everything to string]
# }
# names(eda.gradients.simple)[2:4] = "EDA_Gen_all_" %>% paste0(names(eda.gradients.simple)[2:4])

eda.gradients = bind_rows(eda.df.gen %>% mutate(diagnostic=as.character(diagnostic)), #factor to char, to merge with simple gradients
                          eda.df.gen.simple) %>% 
  mutate(diagnostic = recode_factor(diagnostic, "Eyes"="eyes", "Mouth/Nose"="mn")) %>% 
  #arrange(subject, threat, diagnostic) %>% 
  group_by(subject, diagnostic) %>% 
  #do(temp = gradient.analysis(.$ln_cr)) %>% unnest() %>% group_by(subject, diagnostic) %>% mutate(n=1:n()) %>% spread(n, temp) %>% rename("lds" = `1`, "diff" = `2`, "level" = `3`)
  summarise(lds = gradient.analysis(ln_cr)[1], 
            diff = gradient.analysis(ln_cr)[2], 
            level = gradient.analysis(ln_cr)[3]) %>% 
  gather("measure", "value", lds:level) %>% unite("temp", diagnostic, measure) %>% spread(temp, value) %>% select(subject, contains("lds"), contains("diff"), contains("level"), everything())
names(eda.gradients)[-1] = "EDA_Gen_" %>% paste0(names(eda.gradients)[-1])

# Wide format for correlations --------------------------------------------
eda.wide = merge(eda.gradients %>% ungroup() %>% mutate(subject = as.numeric(subject)), #eda.gradients.simple,
                 eda.ucr %>% filter(include) %>% transmute(subject = as.numeric(subject), EDA_ln_ucr = ln_ucr), 
                 by="subject", all=T)

#all(eda.wide == read_rds("eda.wide.rds" %>% paste0(path.rds, .)), na.rm=T) #check equivalence of processing
#eda.wide %>% write_rds("eda.wide.rds" %>% paste0(path.rds, .))
