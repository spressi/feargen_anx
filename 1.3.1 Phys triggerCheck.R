if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
#source("0 General.R" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/C10 - Gamer Fear/3 Diagnostic Generalization/3 Analysis/main Anx/", .))

toPhysString = function(x, width=2) formatC(x, width=width, format="d", flag="0") #add leading zeros (and convert to string)
toPhysFileName = function(x, width=2, prefix="vp", postfix=".txt") x %>% toPhysString() %>% paste0(prefix, ., postfix)
toCode = function(x, prefix="vp", postfix="") x %>% toPhysString() %>% paste0(prefix, ., postfix)

exclusions.phys.trials = list()
exclusions.phys.trials[[toCode(10)]] =  1:3 + gen1End-.5  #vp10: Gen2 restarted after  3 trials: (skip triggers 131-133)
exclusions.phys.trials[[toCode(31)]] =  1:5               #vp31: Hab  restarted after  5 trials: (skip triggers 1-5)
exclusions.phys.trials[[toCode(41)]] =  1                 #vp41: Hab  restarted after  1 trial
exclusions.phys.trials[[toCode(44)]] =  1:21              #vp44: Hab  restarted after  21 trials
exclusions.phys.trials[[toCode(45)]] =  c(1:3 + acqEnd-.5, 1:4 + 3+gen1End-.5) #vp45: Gen1 restarted after 3 trials & Gen2 restarted after 4 #careful: have to add the 3 excessive Gen1 trials to Gen2, too
exclusions.phys.trials[[toCode(46)]] =  1 + gen1End-.5    #vp46: Gen2 restarted after  1 trial
exclusions.phys.trials[[toCode(50)]] =  1:28              #vp50: Hab  restarted after 28 trials
exclusions.phys.trials[[toCode(82)]] =  1:7 + acqEnd-.5   #vp82: Gen1 restarted after  7 trials
exclusions.phys.trials[[toCode(85)]] =  1:2 + acqEnd-.5   #vp85: Gen1 restarted after  2 trials
exclusions.phys.trials[[toCode(95)]] =  1:5               #vp95: Hab  restarted after  5 trials


breakPositions.theory = floor(c(acqEnd, gen1End))
breaks.theory = length(breakPositions.theory)
maxDistBlock = max(itiEnd)/1000 * sample.rate * 1.2 #max trial time (in seconds) * sampling rate * 20% buffer
for (file in files.phys) {
#for (file in files.phys[c(38)]) {
  data = read.phys(paste0(path.phys, file)) %>% select(Trigger)
  mst = data$Trigger %>% diff() %>% {. > 0} %>% which() %>% {. + 1}
  
  endFlag = ""
  if (exclusions.phys.trials[[file %>% pathToCode()]] %>% is.null() == F) {
    endFlag = " (after manual trial exclusion)" 
    mst = mst[-exclusions.phys.trials[[file %>% pathToCode()]]]
  }
  
  # markerHist = mst %>% diff() %>% hist() %>% .$counts
  # #print(markerHist)
  # separator = {markerHist == 0} %>% which() %>% min() #find first bin with zero frequency
  # breaks.detected = markerHist[separator:length(markerHist)] %>% sum()
  
  markerDist = mst %>% diff()
  breakPositions.detected = which(markerDist > maxDistBlock) #markerDist %>% hist() %>% .$counts
  #markerDist[breakPositions.detected]/sample.rate/60
  breaks.detected = length(breakPositions.detected)
  
  cat(paste0(file %>% pathToCode(), ": ", length(mst), " trials, ", 
               breaks.detected, " break(s) at ", paste(breakPositions.detected, collapse=", "), endFlag), "\n")
  
  if (length(mst)!=trials.n || breaks.theory!=breaks.detected || any(breakPositions.theory!=breakPositions.detected)) 
    warning(file)
}

