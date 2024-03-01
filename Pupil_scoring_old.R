###############################################################################
# Project: Fear Generalization
# Mario Reutter
#
# Extract pupil responses from EDF-files:
# - Trial segmentation
# - Downsampling
# - Low-pass filtering
# - Conversion to mm
# - Save single-trial data (+ grand average)

# INFO: This script requires the edfR package, which is only available for Mac
# Update: eyelinkReader to the rescue!

#require("spatstat")
require("eyelinkReader")
#require("edfR")
require("png")
#require("signal") #don't load because it will mask dplyr::filter
if ("signal" %in% rownames(installed.packages()) == FALSE) install.packages("signal")

options(warn=2)

#rm(list=ls())

etpath <- "C:/Users/jat41gk/Documents/Projekte/Visual Exploration Social Anxiety/Data/EyeLink/" #"/Users/gamer/Documents/UKE/Projekte/SFB-TRR58/Generalization_Mario/Daten/Eyelink/"
path <- "C:/Users/jat41gk/Documents/Projekte/Visual Exploration Social Anxiety/Data/" #"/Users/gamer/Documents/UKE/Projekte/SFB-TRR58/Generalization_Mario/Daten/"
protpath <- "C:/Users/jat41gk/Documents/Projekte/Visual Exploration Social Anxiety/Data/prot/" #"/Users/gamer/Documents/UKE/Projekte/SFB-TRR58/Generalization_Mario/Daten/prot/"
# Save Grand Average under:
savepath <- "C:/Users/jat41gk/Documents/Projekte/Visual Exploration Social Anxiety/Data/Pupil/" #"/Users/gamer/Documents/UKE/Projekte/SFB-TRR58/Generalization_Mario/Daten/prot/"

#vpdat <- read.csv2(paste0(path,"inclusions.eye.csv"))
vpn.eye = fixations$subject %>% unique() %>% setdiff(exclusions.eye.num) %>% sort()
#exclusions.test <- c(41,42,43) # testweise ausschließen
vpn <- vpn.eye #%>% setdiff(exclusions.test) #vpdat$filename

#vpn <- c(48)

hz <- 1000
newhz <- 100

smoothing <- TRUE    # Smooth time series (2 Hz)?
conversion <- TRUE   # Convert to mm according to Hayes & Petrov (2016)
grandaverage <- TRUE # Save grand average across all trials?

# Distance Eye - Screen Center (mm)
distance <- 500 #TODO: check if distance was the same for tower mount and chin rest

# Scoring-Window
st <- -1000; en <- 6000  # Prestimulus-Baseline start and end

# Loop over subjects
ga <- numeric()  # Grand Averages per participant

for (vp in vpn) {
  #code <- substr(as.character(vpn[i]),1,7)
  #vp <- 41
  code <- paste0("vp",ifelse(vp<10,"0",""), vp)
  print(code)
  prot_all <- data.frame()
  
  for (block in 1:3) {  
    print(block)
    #block <- 1
    # Load raw EDF data
    #trials <- edf.trial(paste(etpath,code,".edf",sep=""),samples=T,eventmask=T)
    trials <- read_edf(paste(etpath,code,"_",block,".edf",sep=""),import_samples=TRUE)

    # Determine trial onsets
    msg <- trials$events[grep("Stimulus",trials$events$message),]
    ntrial <- nrow(msg)
    
    # Throw out error when ntrial!=160
    if(block ==1 & ntrial!=60) { print(paste("Incorrect number of trials:",ntrial)) }
    if(block == 2 & ntrial!=70) { print(paste("Incorrect number of trials:",ntrial)) }
    if(block == 3 & ntrial!=70) { print(paste("Incorrect number of trials:",ntrial)) }
    #if (ntrial!=60) { print(paste("Incorrect number of trials:",ntrial)) }
      
    
    # Loop over trials to determine trial-by-trial pupil width
    rawpd <- numeric()
    for (trialnr in 1:ntrial) {
      # Determine onset (in ms)
      #trialnr <- 1
      #trialnr <- 41
      onset <- msg$sttime[trialnr]+1000
      
      # Get pupil data
      trialdat <- trials$samples[(trials$samples$time>(onset+st)) & (trials$samples$time<=(onset+en)),]
      if(length(trialdat$trial) == 0) {
         trialdat.fill<- matrix(NA, ncol = length(trialdat), nrow = 1)
         trialdat.fill <- as.data.frame(trialdat.fill)
         names(trialdat.fill) <- names(trialdat)
         trialdat <- trialdat.fill
      }
      
      #include blinks
      blinksdat <- trials$blinks
      blinksdat <- dplyr::filter(blinksdat, blinksdat$trial == trialnr)
      trialdat$blink <- 0
    
      # solve with loop? 
      for(i in 1:length(blinksdat$trial)) {
        # i <- 1
        for (j in 1:length(trialdat$trial)) {
          #j <- 3479
          #trialdat$paR[j]
          ifelse((trialdat$time[j] >= blinksdat$sttime[i] & trialdat$time[j] <= blinksdat$entime[i]), trialdat$blink[j] <- 1, trialdat$blink[j] <- trialdat$blink[j])
        }
      }
  
      # Use pupil data of the right eye / Exception: Look023 -> left eye was erroneously tracked
      if (length(trialdat$trial) != 7000) {
        f <- matrix(NA, ncol = length(trialdat), nrow = 7000-length(trialdat$trial))
        filling <- data.frame(f)
        names(filling) <- names(trialdat)
        trialdat <- rbind(trialdat, filling) # mit bind_rows() versuchen
      }
      
      if (code=="vp48") {   #left eye was tracked for this subject 
        pd <- trialdat$paL
      } else {
        pd <- trialdat$paR
      }
    
      # Interpolate blinks
      pd[trialdat$blink==1] <- NA
      if (sum(!is.na(pd))>=2) {
        pdi <- approx(trialdat$time,pd,trialdat$time)
      } else {
        pdi <- list()
        pdi$x <- trialdat$time
        pdi$y <- pd
      }
      
      # Downsampling
      pdids       <- suppressWarnings(apply(matrix(pdi$y,ncol=hz/newhz,byrow=TRUE),1,mean))
      timelinenew <- suppressWarnings(apply(matrix(pdi$x,ncol=hz/newhz,byrow=TRUE),1,mean))
      
      # Smooting?
      if (smoothing) {
        # Multiply first and last value to get rid of smoothing artefacts
        pdidsf <- c(rep(pdids[1],100),pdids,rep(pdids[length(pdids)],100))
        pdidsf[is.na(pdidsf)] <- mean(pdidsf,na.rm=TRUE)
        bf <- signal::butter(2, 1/(newhz/2)*2)     # 2 Hz low-pass filter
        pdidssmooth <- signal::filtfilt(bf, pdidsf) # apply filter
        pdidssmooth <- pdidssmooth[101:(length(pdidssmooth)-100)]
      } else {
        pdidssmooth <- pdids
      }
      
      # Conversion to mm?
      # Based on Hayes & Petrov (2016)
      if (conversion) {
        pdidssmooth[pdidssmooth<0] <- 0
        pdidssmooth <- 1.70*10^(-4)*distance*sqrt(pdidssmooth)
      }
      
      # Store data
      rawpd <- rbind(rawpd,pdidssmooth)
    }
    
    # Calculate grand average across trials
    ga <- rbind(ga, apply(rawpd,2,mean,na.rm=TRUE))
    
    # Save proprocessed data
    prot <- data.frame(trial=1:ntrial,rawpd)
    names(prot) <- c("trial",paste("pd",1:(ncol(rawpd)),sep=""))
    write.csv2(prot,paste0(savepath,"Blockwise/",code,"_",block,"_pupil.csv"),row.names=FALSE,quote=FALSE)
    
    zip(paste0(savepath,code,"_",block,"_pupil.zip"),paste0(protpath,code,"_pupil.csv"),flags="-j9X")
    #file.remove(paste0(protpath,"Blockwise/",code,"_pupil.csv"))
    prot_all <- rbind(prot_all, prot)
  }
  prot_all$trial <- 1:(nrow(prot_all))
  write.csv2(prot_all,paste0(savepath,code,"_pupil.csv"),row.names=FALSE,quote=FALSE)
}

options(warn=0)

# Plot single participants
require(ggplot2)

plotdat <- data.frame(subj=rep(1:nrow(ga),each=ncol(ga)),
                      sec=seq(st/1000,en/1000-1/newhz,1/newhz),
                      y=as.numeric(t(ga)))

ggplot(plotdat, aes(y=y, x=sec, group=subj)) +
  geom_line(aes(color=subj)) +
  labs(x = "Time (s)") +
  labs(y = "Pupil change (mm)")

ggsave(paste(savepath,"Pupil-Changes_SingleSubj.png",sep=""),width=12,height=9,units="cm",scale=1.2)

# Plot grand averages
y    <- apply(ga,2,mean,na.rm=TRUE)
se   <- apply(ga,2,sd,na.rm=TRUE)/sqrt(nrow(ga))
sec  <- seq(st/1000,en/1000-1/newhz,1/newhz)
plotdat <- data.frame(sec,y,se)

ggplot(plotdat, aes(y=y, x=sec)) +
  geom_line() +
  geom_ribbon(data=plotdat,aes(ymin=y-se,ymax=y+se),alpha=0.3) +
  labs(x = "Time (s)") +
  labs(y = "Pupil change (mm)")

ggsave(paste(savepath,"Pupil-Changes_GAst0.png",sep=""),width=12,height=9,units="cm",scale=1.2)
