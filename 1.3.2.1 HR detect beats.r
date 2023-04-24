##################################################################
# Project: SFB-TRR 58, C10, Experiment 3
#
# Tested with RStudio 2022.07.1 Build 554 running R 4.2.1
#
# Analyse ECG raw data:
# - high-pass filter (2 Hz default)
# - Look for r-peaks
# - enable manual correction if necessary
# - users may have to switch the read.table statement below to conform to their data format
#
# Extensions:
# showAll: If TRUE, all segments are plotted, else only implausible ones are shown (exceeding the range from minhr to maxhr)
#          Press Enter to proceed to next segment, 
#          Enter a number to change the threshold for amplitude detection (e.g., "0.5")
#          "m" for manual correction mode ("manual")
#          "i" to interpolate rpeaks (based on minHrChange in %)
#          "l" to change segment length (buggy, so rather end script, change segment.length, and restart :/ )
#          "s" to save a screenshot of the current segment
#          "b" for previous segment ("back"),
#          "j", enter, and a number to jump to the specified segment ("jump")
#          "f" to save R-peaks and finalize current file ("finish")
#          "t" to toggle markers on or off
#          "-" revert y-axis for current data file
commandList = c("", "m", "i", "l", "s", "b", "j", "f", "t", "-")
#
# unzipFiles: If TRUE, zip files are searched, unpacked, loaded, and deleted again. If FALSE, txt data will be searched.
#
# If an rpeaks file already exists, rpeaks are loaded and shown within segments.
# Any changes will override loaded files! (Exception: Abort script with "Esc" prior to finishing current subject)
#
# Several minor changes to enhance usability.
#
# TODOs:
# Bug: Interpolated R-peaks cannot be removed in manual mode 
# => type a threshold to override r peaks and then enter manual mode
#
# Bug: switch of segment length during script leads to misaligned r peaks
# => end script, change segment.length, and restart

requirePackage = function(name, load=T) {
  package = as.character(name)
  if (package %in% rownames(installed.packages()) == FALSE) install.packages(package)
  if (load) library(package, character.only=TRUE)
}

requirePackage("tidyverse")
requirePackage("signal", load=F)

#input
if (exists("path.phys")==F)
  path.phys = "TODO enter input path here"
if (path.phys %>% endsWith("/") == F) path.phys = "/" %>% paste0(path.phys, .)

if (exists("exclusions.phys.trials")==F) exclusions.phys.trials = list()
if (exists("sample.rate")==F) sample.rate = 500 #sample rate in Hz
if (exists("trials.n")==F) trials.n = 200 #number of trials that shall be analyzed (if more trials, last ones will be taken); Turn off by assigning: .Machine$integer.max

#output
if (exists("path.rpeaks")==F)
  path.rpeaks = "TODO enter output path here"
if (path.rpeaks %>% endsWith("/") == F) path.rpeaks = "/" %>% paste0(path.rpeaks, .)
path.rpeaks.postfix = "_rpeaks.csv"
dir.create(path.rpeaks, showWarnings=F) #make sure that path.rpeaks exists

path.screenshots = "plots HR/" %>% paste0(path.rpeaks, "../", .) #paste0(ifelse(exists("path.rds"), path.rds, path.phys), .)
dir.create(path.screenshots, showWarnings=F)

#graphics.off()

# Settings ----------------------------------------------------------------
analysis.range = c(-10, 20) #time in seconds before first and after last marker that shall be analyzed

showAll = T #whether all segments shall be shown or only those with problems (see parameters below)
minhr =  50 #minimum plausible heart rate
maxhr = 120 #maximum plausible heart rate

minHrChange = .7 #smallest percentage to which the bpm may drop within a single beat (e.g., 100 bpm to 70 bpm)

showMarkers = F #whether markers should be shown (Recommendation: FALSE for blind scoring)

segment.length = 30 # r-peak detection segments in seconds
high.pass = 2 # high-pass filter (Hz)
ecg.filter = signal::butter(2, high.pass / (sample.rate / 2), "high") #2nd order filter, i.e. 2*6 = 12 dB/oct roll-off

unzipFiles = F #switch to T if files are not in txt but zip format (packed individually)

# Remove outlier before scoring?
limitrange = FALSE
plimit = 0.003


# Functions ---------------------------------------------------------------
findpeaks = function (series, amplitude, dcrit=20) {
  above = which(series > amplitude)
  if (sum(diff(above) > 1) > 0) {
    ast = c(1,(2:length(above))[diff(above)>1])
    if (length(ast)==1) {
      aen = length(above)		
    } else {
      aen = c(ast[2:length(ast)]-1,length(above))		
    }
    
    result = numeric()
    for (i in seq(ast)) {
      result = c(result,above[ast[i]:aen[i]][series[above[ast[i]:aen[i]]]==max(series[above[ast[i]:aen[i]]])][1])
    }
    
    # Remove duplicates
    if (length(result)>1) {
      fresult = c(result[1],result[2:length(result)][diff(result)>dcrit])
    }
  } else {
    fresult = numeric()
  }
  
  fresult
}

# Test if x/y fall into a rectangle
inrect = function(x, y, bx1, by1, bx2, by2) {
  x.low = min(c(bx1, bx2))
  x.high = max(c(bx1, bx2))
  y.low = min(c(by1, by2))
  y.high = max(c(by1, by2))
  
  return((x > x.low) & (x < x.high) & (y > y.low) & (y < y.high))
}

manualcheck = function(t_trial, ecg_trial, rpeak_trial) {
  ready = FALSE # only true if trial was checked completely
  while (!ready) {
    # Check manually (delete systolic peaks within a selected time range)
    ok = FALSE
    while (!ok) {	
      # Plot raw data
      plot(t_trial,ecg_trial,type="l",xlab="time (s)",ylab="ECG")
      box(col="red",lwd=2)
      abline(v=rpeak_trial,col="red")
      if (showMarkers) abline(v=markerwind, col="purple")
      title("Manual editing")
      
      # Plot Heart Rate
      if (length(rpeak_trial)>0) {
        hr = 60/diff(rpeak_trial)
        problem = ifelse((min(hr)<minhr) | (max(hr)>maxhr),1,0)
        plot(t_trial,ecg_trial,type="n",xlab="time (s)",ylab="HR (bpm)", ylim=range(hr))
        lines(rpeak_trial[2:length(rpeak_trial)],hr,col=c("black","red")[problem+1],lwd=2)
      } else {
        plot(t_trial,ecg_trial,type="n",xlab="time (s)",ylab="HR (bpm)",ylim=c(-10,10))
      }
      
      location = locator(n=2)
      
      if (is.null(location) | (length(location$x)==1)) {
        ok = TRUE
      } else {
        # critsys = t_trial[(t_trial>(min(location$x))) & (t_trial<(max(location$x)))]
        # rpeak_trial = rpeak_trial[!(rpeak_trial %in% critsys)]
        rpeak_trial = rpeak_trial[!(rpeak_trial>min(location$x) & rpeak_trial<max(location$x))]
      }
    }
    
    # Check manually (mark r-peaks)
    ok = FALSE
    while (!ok) {	
      # Plot ECG
      plot(t_trial,ecg_trial,type="l",xlab="time (s)",ylab="ECG")
      box(col="green",lwd=2)
      abline(v=rpeak_trial,col="red")
      if (showMarkers) abline(v=markerwind, col="purple")
      title("Manual editing")
      
      # Plot Heart Rate
      if (length(rpeak_trial)>0) {
        hr = 60/diff(rpeak_trial)
        problem = ifelse((min(hr)<minhr) | (max(hr)>maxhr),1,0)
        plot(t_trial,ecg_trial,type="n",xlab="time (s)",ylab="HR (bpm)", ylim=range(hr))
        lines(rpeak_trial[2:length(rpeak_trial)],hr,col=c("black","red")[problem+1],lwd=2)
      } else {
        plot(t_trial,ecg_trial,type="n",xlab="time (s)",ylab="HR (bpm)",ylim=c(-10,10))
      }
      
      # Button
      rect(min(t_trial)-2, min(hr), min(t_trial), max(hr), border=NA, col="lightblue")
      text(min(t_trial)-1, (max(hr)-min(hr))/2+min(hr), "Delete", adj=c(0.5,0.5), srt=90)
      box(lwd=2)    
      
      location = locator(n=1)
      if (is.null(location)) {
        ok = TRUE
        ready = TRUE
      } else {
        if (inrect(location$x, location$y, min(t_trial)-2, min(hr), min(t_trial), max(hr))) {
          ok = TRUE
        } else {
          #critsys = (1:length(t_trial))[(t_trial>(location$x-0.15)) & (t_trial<(location$x+0.15))]
          critsys = (1:length(t_trial))[(t_trial>(location$x-0.02)) & (t_trial<(location$x+0.02))] #use which-function instead?
          rpeak_trial = c(rpeak_trial,t_trial[critsys[ecg_trial[critsys]==max(ecg_trial[critsys])][1]])
          rpeak_trial = sort(rpeak_trial)
          # Delete duplicates
          difs = diff(rpeak_trial)
          rpeak_trial = rpeak_trial[c(TRUE,difs!=0)]
        }
      }
    }
  }
  
  return(rpeak_trial)
}

pathToCode = function(path, path.sep="/", file.ext="\\.") {
  first = path %>% gregexpr(path.sep, .) %>% lapply(max) %>% unlist() %>% {. + 1}
  last = path %>% gregexpr(file.ext, .) %>% lapply(max) %>% unlist() %>% {. - 1}
  return(path %>% substring(first, last))
}

promptIndex = function(codes) {
  index = NA
  while(is.na(index) || index < 1) {
    cat("\n "); cat(paste0(seq(codes), ": ", codes, "\n"))
    prompt = readline("Proceed with subject number: ")
    index = tryCatch({as.integer(prompt)}, error=NA) %>% suppressWarnings()
  }
  
  return(index)
}

implausibility = function(hrtrial, minimum=minhr, maximum=maxhr) {
  #sum(hrtrial < minimum | hrtrial > maximum) #sum is too unspecific
  
  tooLow = hrtrial[hrtrial < minimum]
  tooLow = sum(minimum - tooLow)
  
  tooHigh = hrtrial[hrtrial > maximum]
  tooHigh = sum(tooHigh - maximum)
  
  return(tooLow + tooHigh)
}


subjects.ecg = list.files(path.phys, pattern=ifelse(unzipFiles, ".zip", ".txt"), full.names=TRUE)
codes.ecg = pathToCode(subjects.ecg)

# MAIN --------------------------------------------------------------------
index = 0
while(T) {
  prompt = readline(paste0("Proceed with subject number: ", index+1, "? (Type \"n\" to select subject from list) "))
  index = ifelse(prompt %>% substring(1, 1) == "n", promptIndex(codes.ecg), index+1)
  if (index > length(codes.ecg)) {
    print("Reached end of subject list. Exiting script.")
    break
  }

#for (index in seq(subjects.ecg)) {
  subject = subjects.ecg[index]
  code = codes.ecg[index]
  
  code %>% paste0("Subject: ", .) %>% print()
  print("Loading data...")
  
  #external window (needed for script!)
  x11()
  par(mar=c(2,4,1,1))
  layout(matrix(c(1,2),nrow=2), heights=c(2,1))
  mainWindow = dev.cur()
  
  if (exists("r.threshold")) rm("r.threshold") #clear subject threshold
  
  # Load & Prepare Data -----------------------------------------------------
  # Relevante Variablen:
  # - ecg
  # - timeline
  # - marker (markertime)
  if (unzipFiles) code %>% paste0(".zip") %>% unzip(exdir=path.phys)
  data = read.phys(subject)
  #data = read.table(subject, sep="\t", col.names=c("min", "EDA", "ECG", "Trigger", "empty"), na.strings="", skip=15) %>% select(2:4)
  if (unzipFiles) subject %>% file.remove()
  
  ecg = data$ECG
  ecg[is.na(ecg)] = 0  # Remove NA
  timeline = seq(ecg) / sample.rate      # in Sekunden
  
  mst = data$Trigger %>% diff() %>% {. > 0} %>% which() %>% {. + 1}
  if (length(mst)==0) mst = c(1, nrow(data)) #if no triggers present, score all data
  if (exclusions.phys.trials[[code]] %>% is.null() == F) {
    print(paste0("Manually excluding trials: ", paste0(exclusions.phys.trials[[code]], collapse = ", ")))
    mst = mst[-exclusions.phys.trials[[code]]]
  }
  mst = mst %>% tail(trials.n)
  marker = timeline[mst]
  
  # Restrict range to majority of data
  if (limitrange) {
    ecg[ecg < quantile(ecg,plimit)] = quantile(ecg, plimit)
    ecg[ecg > quantile(ecg,1-plimit)] = quantile(ecg, 1-plimit)
  }
  
  data$ECG = signal::filtfilt(ecg.filter, ecg) #filtering
  rm(ecg)
  
  
  
  # R Scoring -------------------------------------------------------------
  
  # Restrict scoring range
  scorest = {head(marker, 1) + min(analysis.range)} %>% c(0) %>% max()
  scoreen = {tail(marker, 1) + max(analysis.range)} %>% c(tail(timeline, 1)) %>% min()
  analysis.sequence = seq(scorest, scoreen, segment.length)
  
  if (scorest == 0) warning(paste0(code, ": Too little time at the START of subject."))
  if (scoreen == tail(timeline, 1)) warning(paste0(code, ": Too little time at the END of subject."))
  
  # Determine peak detection threshold and find rpeaks adaptively for successive time windows
  print(paste0("Analyzing data: (", round(scorest/60, 2), " - ", round(scoreen/60, 2), " min)"))
  
  allrpeak.list = list() #very important if there is an old list (previous subject) that is longer than the current analysis segment length
  if (code %>% paste0(path.rpeaks, ., path.rpeaks.postfix) %>% file.exists()) { #r peaks exist. Load them from file.
    print("R peaks file detected.")
    allrpeak = read.csv2(code %>% paste0(path.rpeaks, ., path.rpeaks.postfix))[, 1]
    allrpeak.list = list()
    for (i in seq(analysis.sequence)) {
      st = analysis.sequence[i]
      en = {st + segment.length} %>% c(max(timeline)) %>% min()
      allrpeak.list[[toString(i)]] = allrpeak[allrpeak >= st & allrpeak < en]
    }
  } else if (exists("allrpeak")) rm(allrpeak)
  
  manualok = FALSE
  reassign = F
  
  i = 1
  finishFlag = F
  while (i <= length(analysis.sequence) && !finishFlag) {
    windowtxt = paste0(i, " of ", length(analysis.sequence))
    if (exists("r.threshold_new")) rm("r.threshold_new")
    
    st = analysis.sequence[i]
    en = {st + segment.length} %>% c(max(timeline)) %>% min()
    
    timewind = timeline[timeline>=st & timeline<en]
    ecgwind  = data$ECG[timeline>=st & timeline<en]
    markerwind = marker  [marker>=st &   marker<en]
    
    #allrpeak.list needs to be filled a priori since finish command can be prompted, i.e. not all segments are visited and previous file would be shrunken down to visited segments
    # if (exists("allrpeak")) {
    #   rpeak.temp = allrpeak[allrpeak[, 1] >= st & allrpeak[, 1] < en, 1] #take rpeaks of current segment from allrpeak
    #   if (length(rpeak.temp) > 0) allrpeak.list[[toString(i)]] = rpeak.temp
    # }
    
    if (allrpeak.list[[toString(i)]] %>% is.null() %>% !. && !reassign) { #this can also be reached from a "back" prompt even though allrpeak does not exist
      rpeak = allrpeak.list[[toString(i)]]
      hrtrial = 60/diff(rpeak)
      # print("rpeaks read 1")
    } else { #no r peaks in this segment, detect them
      
      #does threshold exist? if not, create one
      if (!exists("r.threshold")) {
        r.threshold.initial = 1
        while (r.threshold.initial==1 || min(hrtrial) < minhr) { # Successively lower threshold until minimum plausible heart rate is found ##TODO: this could go in a separate function with while loop below
          r.threshold.old = as.numeric(quantile(ecgwind, r.threshold.initial)) #determine threshold as upper percentil of time series
          
          r.threshold.initial = r.threshold.initial - 0.01
          r.threshold = as.numeric(quantile(ecgwind, r.threshold.initial)) #determine threshold as upper percentil of time series
          
          rpeak = timewind[findpeaks(ecgwind, r.threshold)]
          # print("rpeaks created 1")
          hrtrial = 60/diff(rpeak)
          
          
          rpeak.old = timewind[findpeaks(ecgwind, r.threshold.old)]
          hrtrial.old = 60/diff(rpeak.old)
          
          #compare number of implausible values
          if(length(rpeak.old) > 1 && implausibility(hrtrial) > implausibility(hrtrial.old)) {
            r.threshold.initial = r.threshold.initial + 0.01
            r.threshold = r.threshold.old #keep old threshold (fewer mistakes)
            break
          }
        }      
      }
      
      #apply threshold
      rpeak = timewind[findpeaks(ecgwind, r.threshold)]
      # print("rpeaks created 2")
      hrtrial = 60/diff(rpeak)
      
      #check for plausibility
      if (min(hrtrial) < minhr || max(hrtrial) > maxhr) { #note: if the while loop would start here, you might get stuck in endless loop (e.g. 99% threshold < minHR but 98% threshold > maxHR)
        if (max(hrtrial) > maxhr) r.threshold.initial = 1 #reset threshold to 1 if maxhr is exceeded
        
        while (r.threshold.initial==1 || min(hrtrial) < minhr) { # Successively lower threshold until minimum plausible heart rate is found ##TODO: this could go in a separate function with while loop above
          r.threshold.old = as.numeric(quantile(ecgwind, r.threshold.initial)) #determine threshold as upper percentil of time series
          
          r.threshold.initial = r.threshold.initial - 0.01
          r.threshold = as.numeric(quantile(ecgwind, r.threshold.initial)) #determine threshold as upper percentil of time series
          
          rpeak = timewind[findpeaks(ecgwind, r.threshold)]
          # print("rpeaks created 3")
          hrtrial = 60/diff(rpeak)
          
          rpeak.old = timewind[findpeaks(ecgwind, r.threshold.old)]
          hrtrial.old = 60/diff(rpeak.old)
          
          #compare number of implausible values
          if(length(rpeak.old) > 1 && implausibility(hrtrial) > implausibility(hrtrial.old)) {
            r.threshold.initial = r.threshold.initial + 0.01
            r.threshold = r.threshold.old #keep old threshold (fewer mistakes)
            rpeak = rpeak.old
            hrtrial = hrtrial.old
            break
          }
        }      
      }
    }
    
    problem = 1 #just to implement do-while loop
    showThis = showAll
    while (showThis || problem) {
      # if (allrpeak.list[[toString(i)]] %>% is.null() %>% !. && !reassign) rpeak = allrpeak.list[[toString(i)]]
      # else  rpeak = timewind[findpeaks(ecgwind, r.threshold)]
      # hrtrial = 60/diff(rpeak)
      
      problem = 0
      if (length(hrtrial)==0 || 
          min(hrtrial) < minhr || max(hrtrial) > maxhr || 
          min(tail(hrtrial / lag(hrtrial), -1)) < minHrChange) {
        problem = 1
      }
      
      if (showThis || problem) {   # Plot trial and ask for different threshold
        # Plot raw data
        plot(timewind, ecgwind, type="l", xlab="time (s)", ylab="ECG")
        if (problem) box(col="orange",lwd=2)
        if (exists("r.threshold")) abline(h=r.threshold, col="blue")
        abline(v=rpeak, col="red")
        if (showMarkers) abline(v=markerwind, col="purple")
        currentTitle = paste0("Subject: ", code, " / Window: ", windowtxt, " (", st," - ", en, " sec)")
        title(currentTitle)
        # Plot heart rate
        bpmRange = range(hrtrial)
        if (Inf %in% bpmRange || -Inf %in% bpmRange) { #special case: no rpeak detected by threshold
          plot(timewind, ecgwind, type="n", xlab="time (s)", ylab="HR (bpm)", ylim=c(-1, 1))
          abline(h=0, col="red", lwd=2)
        } else {
          plot(timewind, ecgwind, type="n", xlab="time (s)", ylab="HR (bpm)", ylim=bpmRange)
          lines(rpeak[2:length(rpeak)], hrtrial, col=c("black","red")[problem+1], lwd=2)
          abline(h=c(minhr, maxhr), lty=2)
        }
        
        reassign = F #needs to be here since it is used to determine if old rpeaks shall be loaded or peak detection with new threshold
        
        promptText = "New threshold"
        if (exists("r.threshold")) promptText = paste0(promptText, " (", round(r.threshold, 2), ")")
        promptText = paste0(promptText, ": ")
        r.threshold_new = readline(promptText)
        invalid = T
        while (invalid) {
          #if (r.threshold_new=="") {
          if (substring(r.threshold_new, 1, 1)=="m") {   # Manual editing
            rpeak = manualcheck(timewind, ecgwind, rpeak)
            allrpeak.list[[toString(i)]] = rpeak
            hrtrial = 60/diff(rpeak)
            title(currentTitle)
            invalid = F
            showThis = T
            # print("rpeaks manual")
            #} else if (substring(r.threshold_new, 1, 1)=="o") {
          } else if (substring(r.threshold_new, 1, 1)=="i") {
            hrMod = tail(hrtrial / lag(hrtrial), -1)
            hrModProblems = which(hrMod < minHrChange) + 1 #rpeaks AFTER these indices are missing
            #hrModProblems = which(hrMod < minHrChange & lead(hrMod) > 1/minHrChange)  + 1 #rpeaks AFTER these indices are missing
            print(paste0("Found ", length(hrModProblems), " problems in heart rate modulation:"))
            
            for (p in hrModProblems) {
              if (p-1 <= 0 || p+1 > length(hrtrial)) 
                print(paste0("Cannot interpolate on the border of a segment. Change segment length via \"l\" to proceed."))
              else { #interpolation algorithm
                neighbors.rpeak = rpeak[p+c(0, 1)] #adjacent valid heart beats
                neighbors.hr = hrtrial[p+c(-1, 1)] #adjacent valid heart rates not affected by missing beat
                missing = min(neighbors.rpeak) + abs(diff(neighbors.rpeak)) * last(neighbors.hr)/sum(neighbors.hr)
                
                rpeak = rpeak %>% c(missing)
              }
            } 
            rpeak = rpeak %>% sort()
            hrtrial = 60/diff(rpeak)
            allrpeak.list[[toString(i)]] = rpeak
            #invalid = F; showThis = F; i = i - 1 #workaround to reload
            invalid = F
            showThis = T
          } else if (r.threshold_new=="") {
            manualok = TRUE
            problem = 0   # Manual acceptance
            invalid = F
            showThis = F
            allrpeak.list[[toString(i)]] = rpeak
          } else if (substring(r.threshold_new, 1, 1)=="b") {
            problem = 0
            invalid = F
            showThis = F
            allrpeak.list[[toString(i)]] = rpeak
          } else if (substring(r.threshold_new, 1, 1)=="s") {
            name = paste0(path.screenshots, code, " ", i, " @", segment.length, ".png")
            print(paste("Saving Screenshot to:", name))
            dev.copy(png, file=name, width=1920, height=1080); dev.off(); dev.set(mainWindow)
          } else if (substring(r.threshold_new, 1, 1)=="f") {
            problem = 0
            invalid = F
            showThis = F
            finishFlag = T
            allrpeak.list[[toString(i)]] = rpeak
          } else if (substring(r.threshold_new, 1, 1)=="t") {
            showMarkers = !showMarkers
            print(paste0("showing markers: ", showMarkers))
            invalid = F
            showThis = T
          } else if (substring(r.threshold_new, 1, 1)=="-") {
            data$ECG = -data$ECG
            ecgwind = -ecgwind
            if (exists("r.threshold")) rm("r.threshold")
            invalid = F
            showThis = T
          } else if (substring(r.threshold_new, 1, 1)=="j") {
            jumpTo = readline("Jump to segment: ") %>% as.integer()
            while (jumpTo %>% is.na()) {
              print(paste0("Error. Cannot read new segment: ", jumpTo))
              jumpTo = readline("Jump to segment: ") %>% as.integer()
            }
            
            problem = 0
            invalid = F
            showThis = F
          } else if (substring(r.threshold_new, 1, 1)=="l") {
            input = readline("Change segment length (seconds): ")
            while (input %>% as.integer() %>% is.na()) {
              print(paste0("Error. Cannot read new segment length: ", jumpTo))
              input = readline("Change segment length (seconds): ")
            }
            
            segment.length = input %>% as.integer()
            allrpeak.list[[toString(i)]] = rpeak
            
            analysis.sequence = seq(scorest, scoreen, segment.length)
            i = which(analysis.sequence <= st) %>% last()
            
            allrpeak = allrpeak.list %>% unlist() %>% unname() %>% unique() %>% sort()
            allrpeak.list = list()
            for (ii in seq(analysis.sequence)) {
              st.ii = analysis.sequence[ii]
              en.ii = {st.ii + segment.length} %>% c(max(timeline)) %>% min()
              allrpeak.list[[toString(ii)]] = allrpeak[allrpeak >= st.ii & allrpeak < en.ii]
            }
            rpeak = allrpeak.list[[toString(i)]]
            hrtrial = 60/diff(rpeak)
            
            st = analysis.sequence[i]
            en = {st + segment.length} %>% c(max(timeline)) %>% min()
            timewind = timeline[timeline>=st & timeline<en]
            ecgwind  = data$ECG[timeline>=st & timeline<en]
            
            problem = 0
            invalid = F
            showThis = T
            #showThis = F; i = i - 1 #workaround to reload
          } else { #no keyword, so new threshold provided as numeric
            r.threshold_new = tryCatch(as.numeric(r.threshold_new), error=NA)
            invalid = is.na(r.threshold_new)
            reassign = !invalid
          }
          
          if (invalid) {
            if (r.threshold_new %in% commandList == F)
              print(paste0("Error. Cannot read new threshold: ", r.threshold_new))
            
            r.threshold_new = readline(paste0("New threshold (", 
                                              ifelse(exists("r.threshold"), round(r.threshold, 2), ""), 
                                              "): "))
          } else if (reassign) {
            r.threshold = as.numeric(r.threshold_new)
            rpeak = timewind[findpeaks(ecgwind, r.threshold)]
            # print("rpeaks created 4")
            hrtrial = 60/diff(rpeak)
            allrpeak.list[[toString(i)]] = rpeak
          }
        }
      }
      if (showThis==F) allrpeak.list[[toString(i)]] = rpeak
    }
    
    #allrpeak.list[[toString(i)]] = rpeak
    if (exists("r.threshold_new")==F || substring(r.threshold_new, 1, 1) != "j") jumpTo = i + 1 #jumpTo has to exist even if not evaluated :/
    i = case_when(
      showAll && substring(r.threshold_new, 1, 1) == "b" ~ {i - 1} %>% c(1) %>% max(),
      showAll && substring(r.threshold_new, 1, 1) == "j" ~ jumpTo %>% c(1) %>% max(),
      TRUE ~ i + 1)
    if (showAll==F) cat(".")
  } 
  
  cat("\n")
  
  allrpeak = allrpeak.list %>% unlist() %>% unname() %>% unique() %>% sort()
  # allrpeak = numeric()
  # for (i in seq(allrpeak.list)) allrpeak = c(allrpeak, allrpeak.list[[toString(i)]])
  
  # Delete duplicates
  keep = c(TRUE, diff(allrpeak) > 60/maxhr/2)
  if (sum(!keep)>0) print(paste0(sum(!keep), " duplicate(s) found and bigger voltage kept (heartrate > ", 2 * maxhr, ")"))
  doubleR = which(!keep)
  keepLater = data$ECG[which(timeline %in% allrpeak[doubleR])] > data$ECG[which(timeline %in% allrpeak[doubleR - 1])]
  keep[doubleR] = keepLater; keep[doubleR-1] = !keepLater
  
  allrpeak = allrpeak[keep]
  
  # Interpolate extrasystoles
  allHr = 60/diff(allrpeak)
  allHrMod = tail(allHr / lag(allHr), -1)
  
  hrModProblems = which(allHrMod < minHrChange & lead(allHrMod) > 1/minHrChange) + 1 #rpeaks AFTER these indices are missing
  interpolate = T
  while (length(hrModProblems) > 0 && interpolate) {
    print(paste0("Found ", length(hrModProblems), " problems in heart rate modulation:"))
    allrpeak[hrModProblems] %>% sapply(function(x) return({x > analysis.sequence} %>% which() %>% tail(1))) %>% 
      paste0("Segment ", ., " (of ", length(analysis.sequence), ")")
    
    interpol = readline("Shall problems be solved by interpolation? (\"y\" to interpolate) ")
    if (substring(interpol, 1, 1)!="y") {
      interpolate = F
      } else {
      for (i in hrModProblems) {
        #allHr[i] #this value is too low because beat is missing between
        #allrpeak[i+0:1]
        
        #missing = mean(allrpeak[i+0:1]) #interpolate by mean time (i.e., allHr = constant; hrMod == 1)
        
        #interpolate by mean heart rate
        missingHr = mean(allHr[i+c(-1, 1)])
        missing = allrpeak[i] + 60/missingHr
        
        #weighted average according to bpm #doesn't work well
        # weights = {allHr[i-1] / (allHr[i+c(-1, 1)] %>% sum())} %>% c(., {1 - .})
        # #sum(weights) #must be 1
        # missing = weighted.mean(allrpeak[i+0:1], weights)
        
        allrpeak = allrpeak %>% c(missing)
      } 
      allrpeak = allrpeak %>% sort()
      allHr = 60/diff(allrpeak)
      allHrMod = tail(allHr / lag(allHr), -1)
      hrModProblems = which(allHrMod < minHrChange & lead(allHrMod) > 1/minHrChange) + 1 #rpeaks AFTER these indices are missing
    }
  }
  
  # Save positions of R-peaks
  out = data.frame(rpeaks = allrpeak)
  write.csv2(out, code %>% paste0(path.rpeaks, ., path.rpeaks.postfix), row.names=FALSE, quote=FALSE)
  
  dev.off()
}


# Error detection ---------------------------------------------------------
#range of heart rate changes (in percent) across all subjects
ranges = list.files(path.rpeaks, full.names = T) %>% 
  .[c(-88)] %>% #exclude due to extrasystoles
  #.[-c(63, 86)] %>% #exclude subjects from range analysis if they have already been checked manually
  sapply(function(x) { allrpeak = read.csv2(x); allHr = 60/diff(allrpeak$rpeaks); allHrMod = tail(allHr / lag(allHr), -1); return(range(allHrMod))})
ranges %>% as.vector() %>% range() %>% {(. - 1) * 100} %>% round(2) %>% paste0("%", collapse=", ") %>% cat("\nRange of all heart range changes (beat-by-beat, in percent): ", ., "\n")
ranges %>% colnames() %>% .[ranges %>% which.min() %>% {. / 2} %>% ceiling()] %>% paste("Min:", .) %>% cat("\n")
ranges %>% colnames() %>% .[ranges %>% which.max() %>% {. / 2} %>% ceiling()] %>% paste("Max:", .) %>% cat("\n")

#find minimum in percental heart range change for last subject
#time.min = allrpeak[which.min(allHrMod)]; {time.min > analysis.sequence} %>% which() %>% tail(1) %>% paste0("Segment ", ., " (of ", length(analysis.sequence), ")")
