if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
#source("0 General.R" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/C10 - Gamer Fear/3 Diagnostic Generalization/3 Analysis/main Anx/", .))

requirePackage("ggforce", load=F) #for drawing circles
requirePackage("png") #for reading png
#requirePackage("raster", load=F) #used to be needed for checkPlot=T

{ # Global Parameters -------------------------------------------------------
  exclusions.eye.num = c() %>% c(exclusions) #a priori exclusions, e.g. calibration not successful
  numToEye = function(nums) return(nums %>% formatC(width=2, format="d", flag="0") %>% paste0("vp", .)) #add leading zeros and prefix "vp"
  exclusions.eye = exclusions.eye.num %>% numToEye()
  
  #baseline validation
  baseline = c(-300, 0) #Baseline in ms relative to stimulus onset; min(baseline) = start; max(baseline) = end
  useAllBaselines = list("83" = 2, "84" = 3) #manually allow all baselines of vp83 phase 2 (Gen1) & vp84 phase 3 (Gen2)
  saveBaselinePlots = T
  driftPlots = T #c("vp30", "vp33")
  maxDeviation_rel = 3 #|z| score
  outlierLimit.eye = .5 #maximum percentage of invalid baselines per subject
  maxSpread = 150 #maximum spread of valid baselines (diameter / edge length)
  usePointDistance = T #use point (vector) distance instead of coordinates independently?
  
  #ROI analysis
  validFixTime.trial = .5 #percentage of valid fixation time within trial in order to analyze trial
  validFixTime.subj = .5 #percentage of trials with sufficient valid fixation time in order to analyze subject
  diagnosticDwell = .5 #percentage of trials per subject that need at least one fixation towards the diagnostic ROI
  showTrialPlots = F
  unifyRois = T
  binResolution = 500 #resolution of temporal bins in ms
  bins = seq(0, trialEnd, binResolution) #bins for temporal dynamic analysis (in ms)
  
  wsx <- 1920 ; wsy <- 1080 #screen resolution (24" ASUS VG248QE)
  pixsize <- 0.27675  # Pixel size in mm
  distance <- 560 # Distance Camera - Eye 560 mm (chin rest with remote tracking)
}

{ # Functions ---------------------------------------------------------------
  outlier_remove <- function(x, z=3) {
    # Remove outlier iteratively
    insample = 1 - is.na(x) #insample <- rep(1,length(x)); insample[is.na(x)] <- 0
    ok <- FALSE
    while (!ok && sum(insample) > 3) {
      xminpos <- (1:length(x))[x==min(x[insample==1], na.rm=T)] #can contain NAs
      xminpos <- xminpos[xminpos %in% (1:length(insample))[insample==1]][1] #eliminates NAs (and takes only first minimum if several are present)
      xmaxpos <- (1:length(x))[x==max(x[insample==1], na.rm=T)] #can contain NAs
      xmaxpos <- xmaxpos[xmaxpos %in% (1:length(insample))[insample==1]][1] #eliminates NAs (and takes only first maximum if several are present)
      tempinsample <- insample; tempinsample[c(xminpos,xmaxpos)] <- 0
      subx <- x[tempinsample==1]
      
      if (x[xminpos] < (mean(subx) - z*sd(subx))) {
        insample[xminpos] <- 0
        out1 <- TRUE
      } else {
        out1 <- FALSE
      }
      
      if (x[xmaxpos] > (mean(subx) + z*sd(subx))) {
        insample[xmaxpos] <- 0
        out2 <- TRUE
      } else {
        out2 <- FALSE
      }
      
      if (!out1 & !out2) { ok <- TRUE }
    }
    return(insample==1)
  }
  
  #removes values until (absolute) deviation is satisfied
  outlier_remove_spread = function(x, deviation, insample=NULL) {
    if (is.null(insample)) insample = T %>% rep(length(x))
    insample[is.na(x)] = F
    
    spread = x[insample] %>% range() %>% diff()
    
    while (spread > deviation) {
      m = x[insample] %>% mean()
      x[!insample] = m #set outliers to m to ignore them in next step
      out.next = {x - m} %>% abs() %>% which.max()
      
      insample[out.next] = F #crucial: use x, not x[insample] to match vector length of insample
      spread = x[insample] %>% range() %>% diff()
    }
    return(insample)
  }
    
  
  #TODO after outlier removal is finished, look once (?) more if "outliers" can be INCLUDED (can happen for biased distributions)
  outlier_remove_abs <- function(x, deviation) {
    # Remove outlier iteratively
    insample = 1 - is.na(x) #insample <- rep(1,length(x)); insample[is.na(x)] <- 0
    ok <- FALSE
    while (!ok && sum(insample) > 2) {
      xminpos <- (1:length(x))[x==min(x[insample==1], na.rm=T)]
      xminpos <- xminpos[xminpos %in% (1:length(insample))[insample==1]][1]
      xmaxpos <- (1:length(x))[x==max(x[insample==1], na.rm=T)]
      xmaxpos <- xmaxpos[xmaxpos %in% (1:length(insample))[insample==1]][1]
      tempinsample <- insample; tempinsample[c(xminpos,xmaxpos)] <- 0
      subx <- x[tempinsample==1]
      
      if (x[xminpos] < (mean(subx) - deviation)) {
        insample[xminpos] <- 0
        out1 <- TRUE
      } else {
        out1 <- FALSE
      }
      
      if (x[xmaxpos] > (mean(subx) + deviation)) {
        insample[xmaxpos] <- 0
        out2 <- TRUE
      } else {
        out2 <- FALSE
      }
      
      if (!out1 & !out2) { ok <- TRUE }
    }
    return(insample==1)
  }
  
  pointDistance = function(x1, y1, x2, y2) {
    return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
  }
  
  outlier_remove_point = function(x, y, z=3, insample=NULL) {
    if (length(x) != length(y)) warning("outlier_remove_point: x- and y-coordinates of unequal length.")
    if (is.null(insample)) insample = T %>% rep(min(c(length(x), length(y))))
    insample[is.na(x) | is.na(y)] = F
    
    ok = F
    while (!ok && sum(insample) > 3) {
      xmean = mean(x[insample]); ymean = mean(y[insample]) #calculate centroid
      distances = pointDistance(x, y, xmean, ymean) #have to be updated every iteration because centroid shifts
      distances[!insample] = 0 #ignore previous outliers for this iteration
      out.next = distances %>% which.max() #find most extreme point (that has not yet been removed)
      insample[out.next] = F #temporarily remove point
      
      if (distances[out.next] <= mean(distances[insample]) + z*sd(distances[insample])) {
        insample[out.next] = T
        ok = T
      }
    }
    
    return(insample)
  }
  
  outlier_remove_point_spread = function(x, y, deviation, insample=NULL) {
    if (length(x) != length(y)) warning("outlier_remove_point_spread: x- and y-coordinates of unequal length.")
    if (is.null(insample)) insample = T %>% rep(min(c(length(x), length(y))))
    insample[is.na(x) | is.na(y)] = F
    
    ok = F
    while (!ok && sum(insample) > 2) {
      xmean = mean(x[insample]); ymean = mean(y[insample]) #calculate centroid
      distances = pointDistance(x, y, xmean, ymean) #have to be updated every iteration because centroid shifts
      distances[!insample] = 0 #ignore previous outliers for this iteration
      out.next = distances %>% which.max() #find most extreme point (that has not yet been removed)
      insample[out.next] = F #temporarily remove point
      
      if (distances[out.next] <= deviation/2) { #deviation = diameter but with distances only radius can be checked => deviation/2
        insample[out.next] = T
        ok = T
      }
    }
    
    return(insample)
  }
  
  pointInBounds = function(point.x, point.y, x.bounds, y.bounds) {
    if (is.na(point.x) || is.na(point.y) || length(point.x)==0 || length(point.y)==0) return(FALSE)
    return(point.x >= min(x.bounds) & point.x <= max(x.bounds) & 
             point.y >= min(y.bounds) & point.y <= max(y.bounds))
  }
  
  degToCm = function(angle, distance) {
    return(distance*tan(angle*pi/180))
  }
  
  degToPix = function(angle, distance, resolution, screenSize) {
    return(degToCm(angle, distance) * resolution / screenSize)
  }
  
  cmToPix = function(cm, resolution, screenSize) {
    return(cm * resolution / screenSize)
  }
  
  cmToDeg = function(cm, distance) {
    return(tan(cm/distance)*180/pi)
  }
  
  # Calculate distance to monitor for deflected fixations
  distfix <- function(x,y,distance) {
    distcalc <- sqrt(x^2+y^2+distance^2)
    return(distcalc)
  }
  
  # Calculate angle between 2 vectors (to compute scanpath)
  angle <- function(x,y){
    dot.prod <- x%*%y 
    norm.x <- norm(x,type="2")
    norm.y <- norm(y,type="2")
    theta <- acos(dot.prod / (norm.x * norm.y))
    theta.deg <- theta/pi*180
    as.numeric(theta.deg)
  }
  
  loadFixations = function(filePath, screen.height) {
    fixations = read_delim(filePath, delim="\t", col_names=F, skip=1, locale=locale(decimal_mark=","), na=".", show_col_types=F) #read.table(filePath, skip=1, dec=",", sep="\t", na.strings=".")
    #fixations = fixations %>% select(-2) #drop 2nd column containing only "Trial: " #only works if sep-parameter is set to default
    names(fixations) = c("subject","trial","start","end","x","y")
    fixations = fixations %>% separate("subject", c("subject", "block")) %>% 
      mutate(subject = subject %>% sub("vp", "", .) %>% as.integer(),
             block = block %>% as.numeric(),
             trial = trial %>% sub("Trial: ","", .) %>% as.numeric(),
             trial = trial + case_when(block==2 ~ acqEnd-.5,
                                       block==3 ~ gen1End-.5,
                                       TRUE ~ 0),
             y = screen.height - y)
    return(fixations)
  }
  
  loadMessages = function(filePath) {
    messages = read_delim(filePath, delim="\t", col_names=F, skip=1, locale=locale(decimal_mark=","), na=".", show_col_types=F) #read.table(filePath, skip=1, dec=",", sep="\t", na.strings=".")
    names(messages) = c("subject","trial","time","event")
    messages = messages %>% separate("subject", c("subject", "block")) %>% 
      mutate(subject = subject %>% sub("vp", "", .) %>% as.integer(),
             block = block %>% as.numeric(), 
             trial = trial %>% sub("Trial: ","", .) %>% as.numeric(),
             trial = trial + case_when(block==2 ~ acqEnd-.5,
                                       block==3 ~ gen1End-.5,
                                       TRUE ~ 0))
    return(messages)
  }
  
  plotRoi = function(roi, matrixRoi=T) {
    #plot(raster::raster(roi)) #error
    
    # raster = reshape2::melt()
    # raster$value = as.factor(raster$value)
    # with(raster, plot(Var2, Var1, col=ifelse(value==1, "green", "white"), asp=1))
    
    if (matrixRoi) #image function expects x then y but matrix is y then x => switch manually
      roi = roi %>% apply(2, rev) %>% t()
    roi %>% ifelse(1, 0) %>% image(col=c("white", "green"), useRaster=T, asp=1)
  }
  
  loadRois = function(folderPath, checkPlot=F) {
    files = list.files(folderPath, ".png")
    rois = vector("list", length(files))
    for (i in 1:length(files)) {
      png = png::readPNG(paste0(folderPath, files[i])) #y (big values are down), x, RGBA
      roi = matrix(png[, , 4]!=0, #select all pixels that are not completely transparent as ROI
                   ncol=dim(png)[2])
      
      if (checkPlot) plotRoi(roi, matrixRoi=T)
        
      rois[[i]] = roi
    }
    return(rois)
  }
  
  #creates new variables that are accessible in the global environment (see "<<-")
  validateBaselines = function(fixs, mess, exclusions, 
                               maxDeviation_rel, maxSpread,
                               saveBaselinePlots=FALSE, postfix="") {
    if (saveBaselinePlots) dir.create(paste0(path.rds, "BL plots/"), showWarnings=F)
    
    vpn = fixs$subject %>% unique() %>% sort() %>% as.character() #all subjects in fixations
    vpn = vpn[!(vpn %in% exclusions)] #minus a priori exclusions
    vpn.n = length(vpn)
    
    trials.n = mess %>% group_by(subject) %>% summarise(trials = max(trial) - min(trial) + 1) %>% .$trials %>% max()
    
    baselines.trial = vector("list", length(vpn.n)) #list of evaluations of baselines per trial for every subject
    baselines.summary = data.frame(subject=character(vpn.n), ntrials=numeric(vpn.n), nValid=numeric(vpn.n), 
                                   invalid=numeric(vpn.n), na=numeric(vpn.n), 
                                   sd_x=numeric(vpn.n), sd_y=numeric(vpn.n),
                                   range_x=numeric(vpn.n), range_y=numeric(vpn.n),
                                   stringsAsFactors=FALSE)
    baselines.summary[,] = NA
    
    print("Validating baselines")
    for (vpIndex in 1:vpn.n) { #TODO make inner of loop enclosed function and call it with apply?
      code = vpn[vpIndex]
      cat(code, " ... ")
      
      # Determine trial number
      fix.vp = fixs %>% filter(subject==code)
      msg.vp = mess %>% filter(subject==code)
      trials.min = min(fix.vp$trial)
      trials.max = max(fix.vp$trial)
      #trials.n = trials.max - trials.min + 1 #must be same for every subject
      
      # baseline
      baseline.vp = data.frame(x=numeric(trials.n), y=numeric(trials.n), condition=character(trials.n), stringsAsFactors=FALSE)
      baseline.vp[,] = NA
      
      # Loop over trials to determine trial-by-trial baselines
      for (trial in fix.vp$trial %>% unique()) {
        # Select trial data
        fix.vp.trial = fix.vp[fix.vp$trial==trial,] #fix.vp %>% filter(trial==trial) #doesn't work
        msg.vp.trial = msg.vp[msg.vp$trial==trial,] #msg.vp %>% filter(trial==trial)
        
        # Determine onset (in ms)
        include <- grep(expoID, msg.vp.trial$event)
        if (length(include) > 0) {
          MsgIndex = grep(expoID, msg.vp.trial$event)
          onset <- msg.vp.trial$time[MsgIndex] %>% as.numeric()
          condition = msg.vp.trial$event[MsgIndex] %>% as.character(); condition = condition %>% substring(nchar(expoID)+1, nchar(condition)-4)
          
          # Subtract onset from startamps
          fix.vp.trial$start  <- fix.vp.trial$start - onset
          fix.vp.trial$end <- fix.vp.trial$end - onset
          
          # Caluculate baseline as weighted average of fixations
          fix.vp.trial.bl <- fix.vp.trial[fix.vp.trial$end>min(baseline) & fix.vp.trial$start<max(baseline),]
          if (nrow(fix.vp.trial.bl) > 0) {
            # Restrict fixation data to baseline
            fix.vp.trial.bl$start[1] <- max(c(min(baseline), first(fix.vp.trial.bl$start)))
            fix.vp.trial.bl$end[nrow(fix.vp.trial.bl)] <- min(c(max(baseline), last(fix.vp.trial.bl$end))) #ifelse(tail(fix.vp.trial.bl$end,1)>blen,blen,tail(fix.vp.trial.bl$end,1)) # = min(blen, tail...) ?
            fix.vp.trial.bl$dur <- fix.vp.trial.bl$end - fix.vp.trial.bl$start
            
            # Calculate baseline coordinates
            xbl <- sum(fix.vp.trial.bl$x*fix.vp.trial.bl$dur)/sum(fix.vp.trial.bl$dur)
            ybl <- sum(fix.vp.trial.bl$y*fix.vp.trial.bl$dur)/sum(fix.vp.trial.bl$dur)
            
            # Store values
            baseline.vp[trial-trials.min+1,] = c(xbl,ybl, condition)
          } else {
            # When no valid fixations are available store NA as baseline for current trial
            baseline.vp[trial-trials.min+1, "condition"] = condition
          }
        }
      }
      
      baseline.vp = baseline.vp %>% mutate(x = as.numeric(x), y = as.numeric(y))
      
      # Determine outlier
      if (usePointDistance) {
        baseline.vp = baseline.vp %>% mutate(
          blok = outlier_remove_point(x, y, maxDeviation_rel) %>%
            outlier_remove_point_spread(x, y, maxSpread, insample=.))
      } else {
        blxok = outlier_remove(baseline.vp$x, maxDeviation_rel) %>% outlier_remove_spread(baseline.vp$x, maxSpread, .)
        blyok = outlier_remove(baseline.vp$y, maxDeviation_rel) %>% outlier_remove_spread(baseline.vp$y, maxSpread, .)
        baseline.vp$blok = blxok & blyok #Baseline is valid when x and y coordinates are ok (i.e. no outlier)
      }
      
      # Store number of valid baselines per subject
      nValid = sum(baseline.vp$blok)
      invalid = 1 - nValid / trials.n
      nas = sum(is.na(baseline.vp$x)) / trials.n
      x.coords = baseline.vp %>% filter(blok) %>% .$x; mean_x = mean(x.coords); sd_x = sd(x.coords); range_x = x.coords %>% range() %>% diff()
      y.coords = baseline.vp %>% filter(blok) %>% .$y; mean_y = mean(y.coords); sd_y = sd(y.coords); range_y = y.coords %>% range() %>% diff()
      baselines.summary[vpIndex, 2:9] = c(trials.n, nValid, invalid, nas, sd_x, sd_y, range_x, range_y); baselines.summary[vpIndex , 1] = code
      baselines.trial[[code %>% as.character()]] = baseline.vp
      
      #plot subject
      if (saveBaselinePlots==TRUE || code %in% saveBaselinePlots || as.character(code) %in% saveBaselinePlots) {
        filename = paste0(path.rds, "BL plots/", postfix, "_", code, ".png")
          
        borders.rel.x = c(mean_x - maxDeviation_rel*sd_x, mean_x + maxDeviation_rel*sd_x)
        borders.rel.y = c(mean_y - maxDeviation_rel*sd_y, mean_y + maxDeviation_rel*sd_y)
        borders.abs.x = c(x.coords %>% range() %>% mean() - maxSpread/2, 
                          x.coords %>% range() %>% mean() + maxSpread/2)
        borders.abs.y = c(y.coords %>% range() %>% mean() - maxSpread/2, 
                          y.coords %>% range() %>% mean() + maxSpread/2)
        
        if (driftPlots==T || code %in% driftPlots || as.character(code) %in% driftPlots) {
          blplot = baseline.vp %>% mutate(trial=1:n()) %>% 
            ggplot(aes(x=x, y=y, color=trial)) + 
            xlim(0, screen.width) + ylim(0, screen.height) + #restrict area to screen
            geom_rect(xmin=0, ymin=0, xmax=screen.width, ymax=screen.height, color="black", fill=NA) +
            geom_hline(yintercept=mean_y, linetype="longdash") + geom_vline(xintercept = mean_x, linetype="longdash") #centroid of valid baselines
          
          #add borders depending on whether point distance was used (circles) or not (rectangles)
          if (usePointDistance) {
            distances = pointDistance(x.coords, y.coords, mean_x, mean_y)
            r_abs = mean(distances) + maxDeviation_rel*sd(distances)
            
            blplot = blplot + #circles
              ggforce::geom_circle(aes(x0=mean_x, y0=mean_y, r=r_abs), color="red", fill=NA, inherit.aes=F) + #borders for validity (relative)
              ggforce::geom_circle(aes(x0=mean_x, y0=mean_y, r=maxSpread/2), color="orange", fill=NA, inherit.aes=F) #borders for validity (absolute)
          } else {
            borders.rel.x = c(mean_x - maxDeviation_rel*sd_x, mean_x + maxDeviation_rel*sd_x)
            borders.rel.y = c(mean_y - maxDeviation_rel*sd_y, mean_y + maxDeviation_rel*sd_y)
            borders.abs.x = c(x.coords %>% range() %>% mean() - maxSpread/2, 
                              x.coords %>% range() %>% mean() + maxSpread/2)
            borders.abs.y = c(y.coords %>% range() %>% mean() - maxSpread/2, 
                              y.coords %>% range() %>% mean() + maxSpread/2)
            
            blplot = blplot + #rectangles
              geom_rect(xmin=min(borders.rel.x), ymin=min(borders.rel.y), xmax=max(borders.rel.x), ymax=max(borders.rel.y), color="red", fill=NA) + #borders for validity (relative)
              geom_rect(xmin=min(borders.abs.x), ymin=min(borders.abs.y), xmax=max(borders.abs.x), ymax=max(borders.abs.y), color="orange", fill=NA) #borders for validity (absolute)
          }
          
          #add points
          blplot = blplot + 
            geom_point() + scale_color_continuous(low="blue", high="green") + #all points (color coded by trial)
            geom_point(data=baseline.vp %>% filter(blok==F), mapping=aes(x=x, y=y), color="red") + #invalid baselines red
            geom_point(x=screen.width/2, y=screen.height/2, shape="+", size=5, color="black") + #fixation cross
            #theme(panel.border = element_rect(color = "black", fill=NA, size=5)) +
            coord_fixed() +
            ggtitle(paste0(code, "_", postfix, " (", round(invalid, digits=2)*100, "% out, ", round(nas, digits=2)*100, "% NAs)"))
          blplot %>% ggsave(filename=filename, plot=., device="png", dpi=300, units="in", width=1920/300, height = 1080/300)
        } else {
          borders.rel.x = c(mean_x - maxDeviation_rel*sd_x, mean_x + maxDeviation_rel*sd_x)
          borders.rel.y = c(mean_y - maxDeviation_rel*sd_y, mean_y + maxDeviation_rel*sd_y)
          borders.abs.x = c(x.coords %>% range() %>% mean() - maxSpread/2, 
                            x.coords %>% range() %>% mean() + maxSpread/2)
          borders.abs.y = c(y.coords %>% range() %>% mean() - maxSpread/2, 
                            y.coords %>% range() %>% mean() + maxSpread/2)
          
          #x11() #plot in new windows (max of 63)
          png(filename, 
              width=1920, height=1080)
          plot(baseline.vp$x, baseline.vp$y, pch=16, col=ifelse(baseline.vp$blok==0, "red", "black"), xlab="x (px)", ylab="y (px)", xlim=c(0, screen.width), ylim=c(0, screen.height), asp=1)
          title(paste0(code, "_", postfix, " (", round(invalid, digits=2)*100, "% out, ", round(nas, digits=2)*100, "% NAs)"))
          abline(h=c(mean_y, screen.height), v=mean_x); abline(v=borders.rel.x, h=borders.rel.y, col="red"); points(x=mean(c(0, screen.width)), y=screen.height/2, pch=3, col="blue")
          
          dev.off()
        }
      }
      
      #invisible(readline(prompt="Baseline created! Press [enter] to continue"))
    }
    
    baselines.trial <<- baselines.trial
    
    print("Baseline validation finished")
    return(baselines.summary)
  }
  
  roiAnalysis = function(fixs, baselines, showTrialPlots=F, offset=0) {
    trials.vp = fixs %>% .$trial %>% unique() %>% length()
    
    eye.vp = data.frame(trial = numeric(trials.vp),
                        condition = numeric(trials.vp),
                        dwell = numeric(trials.vp),
                        dwell.non = numeric(trials.vp),
                        ms = numeric(trials.vp),
                        ms.non = numeric(trials.vp),
                        diagnosticFirst = numeric(trials.vp), 
                        fixN = numeric(trials.vp),
                        mFixTime = numeric(trials.vp),
                        roiSwitch = numeric(trials.vp),
                        scanPath = numeric(trials.vp))
    for (i in 1:(length(bins)-1)) eye.vp = eye.vp %>% mutate(!!paste0("dwell.bin", i) := numeric(trials.vp))
    for (i in 1:(length(bins)-1)) eye.vp = eye.vp %>% mutate(!!paste0("dwell.non.bin", i) := numeric(trials.vp))
    
    bl.mean.x = mean(baselines$x[baselines$blok]) #calculate mean valid baseline for trials that have invalid baseline
    bl.mean.y = mean(baselines$y[baselines$blok])
    
    #for (i in fixs %>% .$trial %>% unique() %>% sort() %>% {. - offset}) {
    for (i in fixs %>% pull(trial) %>% unique() %>% seq()) {
      
      #print(i)
      
      #baseline correction
      trial = baselines[i, ]
      bl.corr.x = ifelse(is.na(trial$blok)==F && trial$blok, trial$x, bl.mean.x) %>% round(digits=1) #round to accuracy of eye tracker
      bl.corr.y = ifelse(is.na(trial$blok)==F && trial$blok, trial$y, bl.mean.y) %>% round(digits=1)
      
      i.total = i + offset
      
      #fixations
      fixations.trial = fixs %>% filter(trial==i.total)
      if (fixations.trial %>% nrow() == 0) {
        eye.vp[i, 1] = i.total; eye.vp[i, -1] = NA
        next
      }
      
      #messages
      condition = fixations.trial %>% transmute(condition = pic %>% sub(".png", "", ., fixed=T)) %>% .$condition %>% head(1)
      
      #rois
      roi = rois[[condition %>% substring(1, 1) %>% as.integer()]]
      roi.non = rois.nondiag[[condition %>% substring(1, 1) %>% as.integer()]]
      
      fixations.trial = fixations.trial %>% 
        mutate(x.cent = x - bl.corr.x,
               y.cent = y - bl.corr.y,
               dist = sqrt(x.cent^2 + y.cent^2), #distance from center
               x = x - bl.corr.x + screen.width/2, 
               y = y - bl.corr.y + screen.height/2,
               x.roi = {x - screen.width/2 + 1 + dim(roi)[2]/2} %>% round(), #old coordinate minus half screen + 1 => middle of screen & roi = 1. this + half roi => 1 = most left column of roi
               y.roi = {screen.height - y - screen.height/2 + 1 + dim(roi)[1]/2} %>% round(), #same as with x but reverse y-direction first (currently up = bigger but down = bigger needed for matrix indexing)
               #inRoi = ifelse(x.roi > 0 & y.roi > 0 & x.roi <= dim(roi)[2] & y.roi <= dim(roi)[1], roi[y.roi, x.roi], FALSE), #if coordinates within bounds of roi matrix, look up in matrix, else FALSE #doesn't work because whole vector is evaluated and throws error even though result FALSE shall be used
               start = start - onset, #start of exposition = time_0
               end = end - onset, #start of exposition = time_0
        )
      
      #roi analysis
      fixations.trial = fixations.trial %>% mutate(
        inRoi = x.roi > 0 & y.roi > 0 & x.roi <= dim(roi)[2] & y.roi <= dim(roi)[1],
        inRoi.non = x.roi > 0 & y.roi > 0 & x.roi <= dim(roi.non)[2] & y.roi <= dim(roi.non)[1]) #fixation on stimulus? (enables indexing without out of bounds error)
      #are fixations on stimuli also within roi?
      fixations.trial$inRoi[fixations.trial$inRoi] = roi[cbind(fixations.trial$y.roi[fixations.trial$inRoi], fixations.trial$x.roi[fixations.trial$inRoi])]
      fixations.trial$inRoi.non[fixations.trial$inRoi.non] = roi.non[cbind(fixations.trial$y.roi[fixations.trial$inRoi.non], fixations.trial$x.roi[fixations.trial$inRoi.non])]
      fixations.trial = fixations.trial %>% mutate(roi = case_when(
        inRoi ~ "diagnostic", inRoi.non ~ "non-diagnostic", TRUE ~ "no ROI") %>% as.factor())
      
      fixations.trial.analysis = fixations.trial %>% #filter(end > 0, start < ratingStart) %>% #replaced by start = ifelse...
        mutate(start = ifelse(start < 0, 0, start),
               end = ifelse(end > ratingStart, ratingStart, end),
               dur = end - start,
               x.st.mm = lag(x.cent) * pixsize, y.st.mm = lag(y.cent) * pixsize, #previous fixation = saccade start
               x.en.mm = x.cent * pixsize, y.en.mm = y.cent * pixsize, #current fixation = saccade end
               x.st.mm = ifelse(x.st.mm %>% is.na(), 0, x.st.mm), x.en.mm = ifelse(x.en.mm %>% is.na(), 0, x.en.mm), #set first fixation to fixation cross
               y.st.mm = ifelse(y.st.mm %>% is.na(), 0, y.st.mm), y.en.mm = ifelse(y.en.mm %>% is.na(), 0, y.en.mm), #set first fixation to fixation cross
               distfix.st = distfix(x.st.mm, y.st.mm, distance),
               distfix.en = distfix(x.en.mm, y.en.mm, distance)
        ) %>% filter(dur > 0) %>% 
        rowwise() %>% mutate(angle = angle(c(x.st.mm, y.st.mm, distfix.st), c(x.en.mm, y.en.mm, distfix.en)))
      
      fixations.trial.analysis.bins = fixations.trial %>% #filter(end > min(bins), start < max(bins)) %>% #replaced by start = ifelse...
        mutate(start = ifelse(start < min(bins), min(bins), start),
               end = ifelse(end > max(bins), max(bins), end),
               dur = end - start) %>% filter(dur > 0)
      
      #output variables
      dwell = with(fixations.trial.analysis, sum(dur[inRoi] / sum(dur)))
      dwell.non = with(fixations.trial.analysis, sum(dur[inRoi.non] / sum(dur[inRoi | inRoi.non])))
      ms = with(fixations.trial.analysis, start[inRoi] %>% head(1)); if (is_empty(ms)) ms = ratingStart
      ms.non = with(fixations.trial.analysis, start[inRoi.non] %>% head(1)); if (is_empty(ms.non)) ms.non = ratingStart
      diagnosticFirst = ms < ms.non
      fixN = fixations.trial.analysis %>% nrow()
      mFixTime = fixations.trial.analysis %>% pull(dur) %>% mean()
      roiSwitch = with(fixations.trial.analysis, roi[roi != "no ROI"]) %>% as.numeric() %>% diff() %>% {. != 0} %>% sum()
      scanPath = with(fixations.trial.analysis, sum(angle, na.rm=T))
      
      dwell.bins = c()
      dwell.non.bins = c()
      for (b in 2:length(bins)) {
        bin = bins[(b-1):b]
        fix.bin = fixations.trial.analysis.bins %>% filter(end >= min(bin) & start <= max(bin)) %>% #only fixations that touch current bin
          mutate(start = ifelse(start <= min(bin), min(bin), start), 
                 end = ifelse(end >= max(bin), max(bin), end),
                 dur = end - start) #prune values to borders of bin to calculate duration correctly
        dwell.bins = with(fix.bin, sum(dur[inRoi] / sum(dur))) %>% c(dwell.bins, .)
        dwell.non.bins = with(fix.bin, sum(dur[inRoi.non] / sum(dur[inRoi | inRoi.non]))) %>% c(dwell.non.bins, .)
      }
      
      if (showTrialPlots) {
        picPath = list.files(paste0(path.rois, ".."), "png", full.names=T); picPath = picPath[picPath %>% grep(condition, .)]
        #picPath = list.files(paste0(path.rois, ".."), "jpg", full.names=T)[condition%/%10 * 2-1 + round((condition%%10-1)/6)]
        fixationCross = head(fixations.trial.analysis, 1) %>% mutate(start = 0, end = 0, x = screen.width/2, y = screen.height/2, roi="no ROI", dur=0)
        fixations.trial.plot = rbind(fixationCross, fixations.trial.analysis) %>% mutate(roi = factor(roi, levels=c("diagnostic", "non-diagnostic", "no ROI")))
        { fixations.trial.plot %>% ggplot(aes(x=x, y=y)) + 
            annotation_custom(grid::rasterGrob(
              png::readPNG(picPath),
              #jpeg::readJPEG(picPath),
              width=unit(1,"npc"), height=unit(1,"npc")), 
              xmin = screen.width/2-642/2, xmax = screen.width/2+642/2, 
              ymin = screen.height/2-676/2, ymax = screen.height/2+676/2) +
            geom_point(aes(x=screen.width/2, y=screen.height/2), color="blue", shape="+", size=5) + #fixation cross
            #geom_point(aes(x=x[2], y=y[2], size=dur[2]), shape=21) + #first fixation after fix cross accentuated
            geom_point(aes(size=dur, color=roi), alpha=.25) + geom_path() + #all fixations as transparent circles with fill
            #geom_point(aes(size=dur, color=inRoi), shape=21) + geom_path() + #all fixations as opaque circles without fill with path
            #geom_rect(mapping=aes(xmin=0, xmax=screen.width, ymin=0, ymax=screen.height), fill=NA, color="black") + #screen border
            #xlim(0, screen.width) + ylim(0, screen.height) + #full screen
            scale_x_continuous(limits=c(screen.width/2-642/2, screen.width/2+642/2), expand=c(0,0)) +
            scale_y_continuous(limits=c(screen.height/2-676/2, screen.height/2+676/2), expand=c(0,0)) +
            coord_fixed() + ggtitle(paste0(code, ": trial ", i.total, ", condition ", condition)) + 
            scale_size(range = c(3, 10)) +
            labs(size = "duration") +
            scale_color_manual(values=c("green3", "red", "grey")) +
            theme_bw() + theme(plot.title = element_text(hjust = 0.5)) } %>% print()
      }
      
      #assign output variables
      eye.vp[i, -2] = c(i.total, dwell, dwell.non, ms, ms.non, diagnosticFirst, fixN, mFixTime, roiSwitch, scanPath, dwell.bins, dwell.non.bins); eye.vp[i, 2] = condition
    }
    
    return(eye.vp)
  }
}


# baseline validation --------------------------------------------------------------------
#ratings = read_rds("ratings.rds" %>% paste0(path.rds, .))
messages = loadMessages(paste0(path.eye, "Messages.txt")) %>% mutate(time = ifelse(time==0, 0, time + 1000)) #correct for duplicated pre-stimulus-baseline
fixations = loadFixations(paste0(path.eye, "Fixations.txt"), screen.height) %>% 
  #filter(subject %in% exclusions.eye.num == F) %>% #exclusion is done in baseline validation below
  left_join(ratings, by=c("subject", "trial", "block")) %>% #get conditions & other variables
  left_join(messages %>% filter(grepl(expoID, event, fixed=T)) %>% 
              rename(onset = time, pic = event) %>% 
              mutate(pic = gsub(expoID, "", pic, fixed=T)) %>% select(-pic), 
            by=c("subject", "trial", "block")) %>% arrange(subject, trial)
vpn.eye = fixations$subject %>% unique() %>% setdiff(exclusions.eye.num) %>% sort()

#determine valid trials by fixation time (invalid trials need not be evaluated by their baseline)
eye.valid.trial = fixations %>% mutate(start = start - preStim, end = end - preStim, #realign such that 0 = stim start
                                       start = ifelse(start < 0, 0, start), #discard fraction of fixation before stimulus
                                       end = ifelse(end > ratingStart, ratingStart, end), #discard fraction of fixation after rating onset
                                       dur = end - start) %>% filter(dur > 0) %>% 
  #filter(block >= 2) %>% #only generalization
  group_by(subject, trial, block) %>% summarise(valid = sum(dur) / ratingStart)
  #group_by(subject) %>% summarise(valid = sum(dur) / ratingStart / length(unique(trial)))
with(eye.valid.trial %>% filter(block==2), hist(valid, breaks=seq(0, 1, length.out=20+1), main=paste0("Valid Fixation Time (Block ", unique(block), ")"))); abline(v=validFixTime.trial, col="red", lwd=3, lty=2)
with(eye.valid.trial %>% filter(block==3), hist(valid, breaks=seq(0, 1, length.out=20+1), main=paste0("Valid Fixation Time (Block ", unique(block), ")"))); abline(v=validFixTime.trial, col="red", lwd=3, lty=2)

#exclude trials with insufficient valid fixations (need not be validated for their baseline)
fixations.valid = eye.valid.trial %>% filter(valid > validFixTime.trial) %>% 
  left_join(fixations %>% mutate(trial = as.numeric(trial)), by=c("subject", "trial", "block")) %>% select(-valid)


baselines.summary1 = validateBaselines(fixations.valid %>% filter(block==1), messages %>% filter(block==1), exclusions.eye, maxDeviation_rel, maxSpread, saveBaselinePlots, postfix="(pre)acq") %>% mutate(block=1) %>% select(subject, block, everything())
baselines.trial1 = baselines.trial; rm(baselines.trial)
baselines.summary2 = validateBaselines(fixations.valid %>% filter(block==2), messages %>% filter(block==2), exclusions.eye, maxDeviation_rel, maxSpread, saveBaselinePlots, postfix="gen1") %>% mutate(block=2) %>% select(subject, block, everything())
baselines.trial2 = baselines.trial; rm(baselines.trial)
baselines.summary3 = validateBaselines(fixations.valid %>% filter(block==3), messages %>% filter(block==3), exclusions.eye, maxDeviation_rel, maxSpread, saveBaselinePlots, postfix="gen2") %>% mutate(block=3) %>% select(subject, block, everything())
baselines.trial3 = baselines.trial; rm(baselines.trial)

# baseline validation summary
baselines.summary = baselines.summary1 %>% bind_rows(baselines.summary2) %>% bind_rows(baselines.summary3) %>% 
  mutate(included = invalid <= outlierLimit.eye & range_x <= maxSpread & range_y <= maxSpread)
rm(baselines.summary1, baselines.summary2, baselines.summary3)

for (b in baselines.summary$block %>% unique() %>% sort()) {
  with(baselines.summary %>% filter(block==b), hist(invalid, breaks=max(ntrials), main=paste0("Valid Baselines (Block ", unique(block), ")"))); abline(v = outlierLimit.eye, col="red", lwd=2, lty=2)
}

baselines.summary %>% group_by(block) %>% summarise(totalN = n(), includedN = sum(included), includedP = mean(included))

baselines.summary %>% filter(block > 1) %>% group_by(subject) %>% summarise(invalid = mean(invalid)) %>% arrange(desc(invalid))
baselines.summary %>% filter(block > 1) %>% group_by(subject) %>% summarise(invalid = mean(invalid)) %>% with(hist(invalid, breaks=20, main="Valid baselines during Generalization")); abline(v = outlierLimit.eye, col="red", lwd=2, lty=2)

baselines.summary %>% filter(block > 1) %>% group_by(subject) %>% summarise(invalid = mean(invalid)) %>% 
  filter(invalid < outlierLimit.eye) %>% summarise(invalid.m = mean(invalid), invalid.sd = sd(invalid))

# Exclusions --------------------------------------------------------------
#baselines
print(eye.invalid.bl <- baselines.summary %>% filter(block > 1, included==F) %>% .$subject %>% unique())

#valid fixations during trial
eye.valid = eye.valid.trial %>% group_by(subject, block) %>% summarise(invalid.n = sum(valid < validFixTime.subj), invalid.p = invalid.n / n())
#eye.valid %>% arrange(desc(invalid.p)) %>% filter(invalid.p > outlierLimit.eye)
print(eye.invalid.fix <- eye.valid %>% filter(block >= 2, invalid.p > validFixTime.subj) %>% .$subject %>% unique())

#exclusion
vpn.eye = vpn.eye %>% setdiff(eye.invalid.bl) %>% setdiff(eye.invalid.fix)

# Create ROIs -------------------------------------------------------------
rois = loadRois(path.rois, checkPlot=F) #diagnostic rois from picture masks

if (unifyRois) { #unify ROIs per diagnostic region: use logical OR (i.e. sum of numeric values)
  rois[[1]] = (apply(rois[[1]], 2, as.numeric) + apply(rois[[3]], 2, as.numeric)) > 0
  rois[[3]] = rois[[1]]
  rois[[2]] = (apply(rois[[2]], 2, as.numeric) + apply(rois[[4]], 2, as.numeric)) > 0
  rois[[4]] = rois[[2]]
}

#non-diagnostic ROIs
rois.nondiag = vector("list", length(rois))
rois.nondiag[[1]] = rois[[length(rois)]]
for (i in seq(rois)[-1]) rois.nondiag[[i]] = rois[[i-1]]

#write_rds(rois, "rois.rds" %>% paste0(path.rds, .)); write_rds(rois.nondiag, "rois.nondiag.rds" %>% paste0(path.rds, .))

# ROI analysis ------------------------------------------------------------
#rois = read_rds("rois.rds" %>% paste0(path.rds, .)); rois.nondiag = read_rds("rois.nondiag.rds" %>% paste0(path.rds, .))
eye.list = vector("list", length(vpn.eye))
for (code in vpn.eye) {
  #code = sample(vpn.eye, 1) #for testing
  cat(paste(code, "... "))
  
  fixations.vp = fixations.valid %>% filter(subject==code)
  baselines.trial1.vp = baselines.trial1[[code %>% as.character()]]
  baselines.trial2.vp = baselines.trial2[[code %>% as.character()]]
  baselines.trial3.vp = baselines.trial3[[code %>% as.character()]]
  
  if(useAllBaselines[code %>% as.character()] %in% 2) baselines.trial2.vp$blok=T
  if(useAllBaselines[code %>% as.character()] %in% 3) baselines.trial3.vp$blok=T
  
  eye.vp1 = tibble() #eye data of first block not used & subject 54: ratings of first block missing => conditions NA => error
  # eye.vp1 = roiAnalysis(fixations.vp %>% filter(block==1) %>% subset(subject %in% exclusions.onlyGen == F), #exclude habituation & acquisition for onlyGen subjects (see "0 General.R")
  #                      baselines.trial1.vp, showTrialPlots=showTrialPlots)
  eye.vp2 = roiAnalysis(fixations.vp %>% filter(block==2), baselines.trial2.vp, showTrialPlots=showTrialPlots, 
              offset=nrow(baselines.trial1.vp))
  eye.vp3 = roiAnalysis(fixations.vp %>% filter(block==3), baselines.trial3.vp, showTrialPlots=showTrialPlots, 
              offset=nrow(baselines.trial1.vp) + nrow(baselines.trial2.vp))
  
  eye.list[[code %>% as.character()]] = eye.vp1 %>% bind_rows(eye.vp2) %>% bind_rows(eye.vp3)
  if (showTrialPlots) invisible(readline(prompt="Press [enter] to continue")) #TODO save to file instead
}


# Tidy data frame(s) ------------------------------------------------------
eye = eye.list %>% bind_rows(.id="subject") %>% mutate(subject = subject %>% as.integer()) %>% 
  filter(condition != 0) #if no valid fixation during trial
rm(eye.list)

#conditions = read_rds("conditions.rds" %>% paste0(path.rds, .)); conditions.csp = read_rds("conditions.csp.rds" %>% paste0(path.rds, .))
eye = eye %>% separate(condition, c("pair", "morph"), convert=T) %>% 
  left_join(conditions, by=c("subject", "pair")) %>% #inject which end is CS+ (1 -> morph 100, 2 -> morph 0)
  #filter(subject %in% vpn.eye) %>% #apply exclusions #already done above
  mutate(morph = as.integer(morph),
         pair = pair %>% as.integer(),
         threat = morph %>% {. / 20 + 1}, #pretend 100 would always be CS+ (true if csp==1)
         threat = ifelse(csp==2, threatLevels.n+1 - threat, threat), #check if threat level has to be inverted (because csp==2)
         threat_num = threat, threat = as.factor(threat),
         diagnostic = as.factor(ifelse(pair %% 2 == 0, "Mouth/Nose", "Eyes")),
         sex = as.factor(ifelse(pair > 2, "female", "male")),
         phase = c("Hab", "Acq", "Gen") %>% cut(trial, breaks=c(0, preAcqEnd, acqEnd, Inf), labels=.),
         pair = as.factor(pair),
         dwell.rois = dwell / (dwell + dwell.non)) %>% select(1:5, dwell.rois, everything()) %>% arrange(subject, trial)

#exclusion based on diagnostic dwell
print(eye.invalid.dwell <- eye %>% group_by(subject) %>% summarise(nBad = sum(dwell == 0), nTotal = n(), valid = 1 - nBad/nTotal) %>% arrange(desc(nBad)) %>% 
  filter(valid < diagnosticDwell) %>% .$subject)
vpn.eye = vpn.eye %>% setdiff(eye.invalid.dwell)
eye = eye %>% filter(subject %in% vpn.eye) %>% merge(questionnaires, by="subject", all.x=T) %>% tibble()

#write_rds(eye, "eye.rds" %>% paste0(path.rds, .))

# Prepare data for analyses -----------------------------------------------
#eye = read_rds("eye.rds" %>% paste0(path.rds, .))
eye.gen = eye %>% filter(phase == "Gen") %>% group_by(subject) %>% 
  summarise(dwell.m=mean(dwell, na.rm=T), dwell.se=se(dwell, na.rm=T),
            dwell.non.m=mean(dwell.non, na.rm=T), dwell.non.se=se(dwell.non, na.rm=T),
            dwell.rois.m=mean(dwell.rois, na.rm=T), dwell.rois.se=se(dwell.rois, na.rm=T),
            ms.m = mean(ms, na.rm=T), ms.se = se(ms, na.rm=T),
            ms.non.m = mean(ms.non, na.rm=T), ms.non.se = se(ms.non, na.rm=T),
            diagnosticFirst.m = mean(diagnosticFirst, na.rm=T), diagnosticFirst.se = se(diagnosticFirst, na.rm=T),
            roiSwitch.m = mean(roiSwitch, na.rm=T), roiSwitch.se = se(roiSwitch, na.rm=T),
            SPAI = mean(SPAI)) %>% 
  select("subject", contains(".m"), everything())

eye.gen.diagnostic = eye %>% filter(phase == "Gen") %>% group_by(subject, diagnostic) %>% 
  summarise(dwell.m=mean(dwell, na.rm=T), dwell.se=se(dwell, na.rm=T),
            dwell.non.m=mean(dwell.non, na.rm=T), dwell.non.se=se(dwell.non, na.rm=T),
            dwell.rois.m=mean(dwell.rois, na.rm=T), dwell.rois.se=se(dwell.rois, na.rm=T),
            ms.m = mean(ms, na.rm=T), ms.se = se(ms, na.rm=T),
            ms.non.m = mean(ms.non, na.rm=T), ms.non.se = se(ms.non, na.rm=T),
            diagnosticFirst.m = mean(diagnosticFirst, na.rm=T), diagnosticFirst.se = se(diagnosticFirst, na.rm=T),
            roiSwitch.m = mean(roiSwitch, na.rm=T), roiSwitch.se = se(roiSwitch, na.rm=T)) %>% 
  select("subject", "diagnostic", contains(".m"), everything())


eye.diagnosticity = eye %>% filter(phase=="Gen") %>% select(-contains("bin")) %>% 
  gather(Diagnosticity, relDwell, c("dwell", "dwell.non")) %>% 
  mutate(ROI = as.factor(ifelse(Diagnosticity=="dwell", diagnostic, ifelse(diagnostic=="Eyes", 2, 1))) %>% recode_factor(`1` = "Eyes", `2` = "Mouth/Nose"),
         Diagnosticity = ifelse(as.numeric(ROI)==as.numeric(diagnostic), "Diagnostic", "Non-Diagnostic"))

eye.diagnosticity.bins = eye %>% filter(phase=="Gen") %>% 
  gather(Diagnosticity, relDwell, contains("dwell.bin"), contains("dwell.non.bin")) %>% 
  mutate(ROI = "dwell.non.bin" %>% grepl(Diagnosticity) %>% ifelse(ifelse(diagnostic=="Eyes", 2, 1), diagnostic) %>% as.factor() %>% recode_factor(`1` = "Eyes", `2` = "Mouth/Nose"),
         #bin = Diagnosticity %>% substr(., nchar(.), nchar(.)) %>% as.numeric() #doesn't work if >= 10 bins
         bin = Diagnosticity %>% gsub("dwell.bin", "", .) %>% gsub("dwell.non.bin", "", .) %>% as.numeric(),
         bin = bins[bin+1]/1000, #convert bin number to seconds
         Diagnosticity = ifelse(as.numeric(ROI)==as.numeric(diagnostic), "Diagnostic", "Non-Diagnostic"))
#eye.diagnosticity.bins[eye.diagnosticity.bins %>% .$relDwell %>% is.na() %>% which(), ]
eye.diagnosticity.bins[eye.diagnosticity.bins %>% .$relDwell %>% is.na() %>% which(), "relDwell"] = 0

eye.diagnosticity.bins.average = eye %>% filter(phase=="Gen") %>%
  gather(Diagnosticity, relDwell, "dwell", "dwell.non") %>% 
  mutate(ROI = "dwell.non" %>% grepl(Diagnosticity) %>% ifelse(ifelse(diagnostic=="Eyes", 2, 1), diagnostic) %>% as.factor() %>% recode_factor(`1` = "Eyes", `2` = "Mouth/Nose"),
         bin = (ratingStart + binResolution*1.5) / 1000,
         Diagnosticity = ifelse(as.numeric(ROI)==as.numeric(diagnostic), "Diagnostic", "Non-Diagnostic")) %>% 
  #select(subject, trial, diagnostic, relDwell, ROI, Diagnosticity, bin) %>%
  select(., which(names(.) %in% names(eye.diagnosticity.bins))) %>% #only keep columns that are in bins dataframe
  bind_rows(eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000), .)
#eye.diagnosticity.bins.average[eye.diagnosticity.bins %>% .$relDwell %>% is.na() %>% which(), ]

eye.diagnosticity.ms = eye %>% filter(phase=="Gen") %>% select(-contains("bin")) %>% 
  gather(Diagnosticity, ms, c("ms", "ms.non")) %>% 
  mutate(ROI = as.factor(ifelse(Diagnosticity=="ms", diagnostic, ifelse(diagnostic=="Eyes", 2, 1))) %>% recode_factor(`1` = "Eyes", `2` = "Mouth/Nose"),
         Diagnosticity = ifelse(as.numeric(ROI)==as.numeric(diagnostic), "Diagnostic", "Non-Diagnostic"))

# Hypotheses Dwell --------------------------------------------------------
#exclusions.eye.dwell = c(18, 42) #TODO don't exclude but prune subjects?
exclusions.eye.dwell = c() #results without exclusions

eye.diagnosticity %>% 
  filter(subject %in% exclusions.eye.dwell == F) %>% #manual exclusion because of extreme dwell
  #mutate(SPAI = scale(SPAI), STAI = scale(STAI)) %>% #no effect - done implicitly?
  ez::ezANOVA(dv=.(relDwell), wid=.(subject),
              within=.(ROI, Diagnosticity), within_full=.(threat),
              #within=.(threat, ROI, Diagnosticity),
              #between=.(SPAI), observed=SPAI,
              between=.(STAI), observed=STAI,
              detailed=T, type=2) %>% apa::anova_apa(force_sph_corr=T)

eye.diagnosticity.subj = eye.diagnosticity %>% group_by(diagnostic, Diagnosticity, ROI, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
  mutate(relDwell = relDwell*100) #improves readability of y-axis
eye.diagnosticity.subj.analysis = eye.diagnosticity.subj  %>% filter(subject %in% exclusions.eye.dwell == F) #manual exclusion because of extreme dwell
print(eye.main <- eye.diagnosticity.subj.analysis %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
        ggplot(aes(x=diagnostic, y=relDwell, color=ROI, shape=Diagnosticity)) + 
        #ggplot(aes(x=diagnostic, y=relDwell, color=as.numeric(ROI)==as.numeric(diagnostic))) + 
        geom_dotplot(data=eye.diagnosticity.subj.analysis %>% filter(ROI=="Eyes"), mapping=aes(group=interaction(ROI, diagnostic), fill=ROI), stackdir="down", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.75) +
        geom_dotplot(data=eye.diagnosticity.subj.analysis %>% filter(ROI!="Eyes"), mapping=aes(group=interaction(ROI, diagnostic), fill=ROI), stackdir="up", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.75) +
        geom_point(size=6, position=dodge) + geom_errorbar(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96), size=2, position=dodge) +
        #scale_x_discrete(labels=c("Eyes", "Mouth/Nose")) +
        scale_color_discrete(labels=c("Eyes", "Mouth/Nose")) +
        #scale_color_discrete(name="Diagnosticity of ROI", labels=c("Diagnostic", "Non-Diagnostic")) +
        #scale_color_viridis_d(labels=c("Eyes", "Mouth/Nose")) +
        ylab("Relative Dwell Time (%)") + xlab("Diagnostic Region") +
        ylim(c(0, 100)) + myGgTheme)

# Figure in German 
# eye.diagnosticity.subj.analysis <- eye.diagnosticity.subj.analysis %>%
#   mutate(ROI =
#            case_when(
#              ROI == "Eyes"~ "Augen",
#              ROI == "Mouth/Nose" ~ "Mund/Nase"
#            )) %>%
#   mutate(Diagnosticity =
#            case_when(
#              Diagnosticity == "Diagnostic" ~ "Diagnostisch",
#              Diagnosticity == "Non-Diagnostic" ~ "Nicht-Diagnostisch"
#            ))%>%
#   mutate(diagnostic =
#            case_when(
#              diagnostic == "Eyes" ~ "Augen",
#              diagnostic == "Mouth/Nose" ~ "Mund/Nase"
#            ))%>%
#   rename(Diagnostizität = Diagnosticity)
# print(eye.main <- eye.diagnosticity.subj.analysis %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
#         ggplot(aes(x=diagnostic, y=relDwell, color=ROI, shape=Diagnostizität)) + 
#         #ggplot(aes(x=diagnostic, y=relDwell, color=as.numeric(ROI)==as.numeric(diagnostic))) + 
#         geom_dotplot(data=eye.diagnosticity.subj.analysis %>% filter(ROI=="Augen"), mapping=aes(group=interaction(ROI, diagnostic), fill=ROI), stackdir="down", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.75) +
#         geom_dotplot(data=eye.diagnosticity.subj.analysis %>% filter(ROI!="Augen"), mapping=aes(group=interaction(ROI, diagnostic), fill=ROI), stackdir="up", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.75) +
#         geom_point(size=6, position=dodge) + geom_errorbar(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96), size=2, position=dodge) +
#         #scale_x_discrete(labels=c("Eyes", "Mouth/Nose")) +
#         scale_color_discrete(labels=c("Augen", "Mund/Nase")) +
#         #scale_color_discrete(name="Diagnosticity of ROI", labels=c("Diagnostic", "Non-Diagnostic")) +
#         #scale_color_viridis_d(labels=c("Augen", "Mund/Nase")) +
#         ylab("Relative Verweildauer (%)") + xlab("Diagnostische Region") +
#         ylim(c(0, 100)) + myGgTheme)

eye.diagnosticity.subj %>% filter(Diagnosticity=="Non-Diagnostic") %>% mutate(relDwell.z = scale(relDwell)[,1]) %>% arrange(desc(abs(relDwell.z))) #exclude subject 18 & 42
#eye.diagnosticity.subj %>% filter(Diagnosticity=="Non-Diagnostic", ROI=="Mouth/Nose") %>% mutate(relDwell.z = scale(relDwell)[,1]) %>% arrange(desc(abs(relDwell.z))) #exclude subject 42?

#2-way interaction (n.s.)
#eye.diagnosticity %>% group_by(diagnostic, ROI) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% summarise(relDwell.anyROI = sum(relDwell))

#Eye Bias descriptive
eye.diagnosticity %>% filter(subject %in% exclusions.eye.dwell == F) %>% #manual exclusion because of extreme dwell
  group_by(ROI, subject) %>% summarise(relDwell = mean(relDwell, na.rm=T)) %>% summarise(relDwell.se = se(relDwell, na.rm=T), relDwell = mean(relDwell)) %>% select(relDwell, everything())

#Diagnosticity Bias descriptive
eye.diagnosticity %>% filter(subject %in% exclusions.eye.dwell == F) %>% #manual exclusion because of extreme dwell
  group_by(Diagnosticity, subject) %>% summarise(relDwell = mean(relDwell, na.rm=T)) %>% summarise(relDwell.se = se(relDwell, na.rm=T), relDwell = mean(relDwell)) %>% select(relDwell, everything())

#SPAI x Diagnosticity
eye.diagnosticity.spai = eye.diagnosticity %>% group_by(subject, Diagnosticity, SPAI) %>% 
  summarise(relDwell.se = se(relDwell, na.rm=T), relDwell = mean(relDwell, na.rm=T))
eye.diagnosticity.spai %>% group_by(Diagnosticity) %>% summarise(cor.test(relDwell, SPAI, alternative="two.sided") %>% apa::cor_apa(r_ci=T, print=F))
eye.diagnosticity.spai %>% ggplot(aes(x=SPAI, y=relDwell, color=Diagnosticity, fill=Diagnosticity, group=Diagnosticity)) +
  facet_wrap(vars(Diagnosticity)) +
  #facet_wrap(vars(ROI), labeller=labeller(ROI=c("Eyes" = "Augen", "Mouth/Nose" = "Mund/Nase"))) +
  scale_color_viridis_d() + scale_fill_viridis_d() +
  geom_errorbar(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96)) +
  stat_smooth(method="lm", color = "black") +
  geom_point(size=4) + 
  ylab("Relative Dwell Time (%)") + myGgTheme

#SPAI x Diagnosticity in German
eye.diagnosticity.spai = eye.diagnosticity %>% group_by(subject, Diagnosticity, SPAI) %>% 
  summarise(relDwell.se = se(relDwell, na.rm=T), relDwell = mean(relDwell, na.rm=T))
eye.diagnosticity.spai %>% group_by(Diagnosticity) %>% summarise(cor.test(relDwell, SPAI, alternative="two.sided") %>% apa::cor_apa(r_ci=T, print=F))
eye.diagnosticity.spai <- eye.diagnosticity.spai %>%
  mutate(Diagnosticity =
                    case_when(
                      Diagnosticity == "Diagnostic" ~ "Diagnostisch",
                      Diagnosticity == "Non-Diagnostic" ~ "Nicht-Diagnostisch"
                    ))%>%
  rename(Diagnostizität = Diagnosticity)

eye.diagnosticity.spai %>% ggplot(aes(x=SPAI, y=relDwell, color=Diagnostizität, fill=Diagnostizität, group=Diagnostizität)) +
  facet_wrap(vars(Diagnostizität)) +
  #facet_wrap(vars(ROI), labeller=labeller(ROI=c("Eyes" = "Augen", "Mouth/Nose" = "Mund/Nase"))) +
  scale_color_viridis_d() + scale_fill_viridis_d() +
  geom_errorbar(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96)) +
  stat_smooth(method="lm", color = "black") +
  geom_point(size=4) + 
  ylab("Relative Verweildauer (%)") + myGgTheme

#STAI x Diagnosticity
eye.diagnosticity.stai = eye.diagnosticity %>% group_by(subject, Diagnosticity, STAI) %>% 
  summarise(relDwell.se = se(relDwell, na.rm=T), relDwell = mean(relDwell, na.rm=T))
eye.diagnosticity.stai %>% group_by(Diagnosticity) %>% summarise(cor.test(relDwell, STAI, alternative="two.sided") %>% apa::cor_apa(r_ci=T, print=F))
eye.diagnosticity.stai %>% ggplot(aes(x=STAI, y=relDwell, color=Diagnosticity, fill=Diagnosticity, group=Diagnosticity)) +
  facet_wrap(vars(Diagnosticity)) +
  #facet_wrap(vars(ROI), labeller=labeller(ROI=c("Eyes" = "Augen", "Mouth/Nose" = "Mund/Nase"))) +
  scale_color_viridis_d() + scale_fill_viridis_d() +
  geom_errorbar(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96)) +
  stat_smooth(method="lm", color = "black") +
  geom_point(size=4) + 
  ylab("Relative Dwell Time (%)") + myGgTheme

# Hypotheses ROI switches ---------------------------------------------------
#exclusions.eye.switch = c(18, 62) #TODO don't exclude but prune subjects?
exclusions.eye.switch = c() #results without exclusions

# small ANOVA
eye.gen %>% 
  filter(subject %in% exclusions.eye.switch == F) %>% #manual exclusion because of extreme switching??
  ez::ezANOVA(dv=.(roiSwitch.m), wid=.(subject),
              #within=.(threat, Diagnosticity),
              between=.(SPAI), observed=SPAI,
              #between=.(STAI), observed=STAI,
              detailed=T, type=2) %>% apa::anova_apa(force_sph_corr=T)
eye.gen %>% with(cor.test(roiSwitch.m, SPAI, alternative="greater")) %>% correlation_out()

# bigger ANOVA
eye.diagnosticity %>% 
  filter(subject %in% exclusions.eye.switch == F) %>% #manual exclusion because of extreme switching??
  ez::ezANOVA(dv=.(roiSwitch), wid=.(subject),
              within =.(diagnostic, threat), #within_full= .(threat, trial),
              #within=.(threat, ROI, Diagnosticity),
              between=.(SPAI), observed=SPAI,
              #between=.(STAI), observed=STAI,
              detailed=T, type=2) %>% apa::anova_apa(force_sph_corr=T)

#SPAI main effect
eye.gen.spai = eye.gen %>% 
  filter(subject %in% exclusions.eye.switch == F) %>% #manual exclusion because of extreme latency
  group_by(SPAI, subject) %>% summarise(switches.se = se(roiSwitch.m, na.rm=T), roiSwitch.m = mean(roiSwitch.m, na.rm=T))
eye.gen.spai %>% with(cor.test(roiSwitch.m, SPAI, alternative="two.sided")) %>% correlation_out()
eye.gen.spai %>% 
  ggplot(aes(y=roiSwitch.m, x=SPAI, color=SPAI)) +
  geom_errorbar(aes(ymin=roiSwitch.m-switches.se*1.96, ymax=roiSwitch.m+switches.se*1.96), width=.05) +
  geom_smooth(method="lm", color="black") + geom_point(size=4) +
  ylab("Mean number of Switches") +
  scale_color_viridis_c() + myGgTheme + theme(legend.position = "none")

#Diagnostic Region main effect
eye.diagnosticity.diagnostic = eye.diagnosticity %>% 
  filter(subject %in% exclusions.eye.switch == F) %>% #manual exclusion because of extreme latency
  group_by(subject,SPAI,diagnostic) %>% summarise(switches.se = se(roiSwitch, na.rm=T), roiSwitch.m = mean(roiSwitch, na.rm=T))
eye.diagnosticity.diagnostic  %>% t.test(roiSwitch.m ~ diagnostic, ., paired=T) %>% apa::t_apa(es_ci=T)
eye.diagnosticity.diagnostic %>%
  ggplot(aes(y=roiSwitch.m, x=SPAI, color= as.factor(diagnostic), fill=as.factor(diagnostic), group=diagnostic)) +
  #facet_wrap(vars(diagnostic)) +
  #geom_errorbar(aes(ymin=roiSwitch.m-switches.se*1.96, ymax=roiSwitch.m+switches.se*1.96), width=.05) +
  scale_color_discrete(labels=c("Eyes", "Mouth/Nose")) +
  geom_smooth(method="lm", color="black") + 
  geom_point(size=4, alpha = 0.3) +
  ylab("Mean number of Switches") +
  #scale_color_viridis_c() + 
  myGgTheme +
  labs(color='Diagnostic Region', fill = 'Diagnostic Region' ) 

#threat main effect
eye.diagnosticity.threat.subject = eye.diagnosticity %>% 
  filter(subject %in% exclusions.eye.switch == F) %>% #manual exclusion because of extreme latency
  group_by(subject, threat) %>% summarise(switches.se = se(roiSwitch, na.rm=T), roiSwitch.m = mean(roiSwitch, na.rm=T))
for (i in (min(as.numeric(eye.diagnosticity.threat.subject$threat))+1):max(as.numeric(eye.diagnosticity.threat.subject$threat))) {
  levels = c(i-1, i)
  cat(paste0("\n\nComparing levels: ", paste(levels, collapse=" vs. "), "\n"))
  eye.diagnosticity.threat.subject %>% filter(threat %in% levels) %>% 
    t.test(roiSwitch.m ~ threat, ., paired=T) %>% apa::t_apa(es_ci=T)
}
eye.diagnosticity.threat <- eye.diagnosticity %>% 
  filter(subject %in% exclusions.eye.switch == F) %>% #manual exclusion because of extreme latency
  group_by(threat) %>% summarise(switches.se = se(roiSwitch, na.rm=T), roiSwitch.m = mean(roiSwitch, na.rm=T))
print(eye.diagnosticity.threat.plot <- eye.diagnosticity.threat %>% ggplot(aes(x=threat, y=roiSwitch.m, color=threat, group=NA)) + 
        geom_dotplot(data=eye.diagnosticity.threat.subject, mapping=aes(group=threat, fill=threat), binaxis="y", alpha=.25, color="black", stackratio=1, stackdir="centerwhole", dotsize=.5) +
        #geom_path(data=heart.ga.gen %>% filter(threat %in% c(1, 6)), aes(group=NA), color = "black", size=1.5) + #generalization line
        geom_errorbar(aes(ymin=roiSwitch.m-switches.se*1.96, ymax=roiSwitch.m+switches.se*1.96), size=1.5) +
        geom_line(size=1) + 
        geom_point(size=4.5) +
        scale_color_manual(values=colors)+#, guide=guide_legend(reverse=T)) + scale_y_reverse() +
        #scale_fill_manual(values=colors, guide=guide_legend(reverse=T)) +
        scale_fill_manual(values=rep("grey", 6), guide=guide_legend(reverse=T)) +
        scale_x_discrete(labels=c("CS-", paste0("GS", 1:4), "CS+")) +
        ylab("Mean number of Switches") + xlab("Threat") + labs(color="Threat") +
        myGgTheme + theme(
          #aspect.ratio = 1,
          legend.position = "none",
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, "pt"))))
eye.diagnosticity.threat %>% select(threat, roiSwitch.m, switches.se) #descriptive values

# Hypotheses Latency ------------------------------------------------------
#exclusions.eye.ms = c(18, 62) #TODO don't exclude but prune subjects?
exclusions.eye.ms = c() #results without exclusions

eye.diagnosticity.ms %>% filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
  ez::ezANOVA(dv=.(ms), wid=.(subject),
              within=.(ROI, Diagnosticity), within_full=.(threat, trial),
              #within=.(threat, ROI, Diagnosticity),
              between=.(SPAI), observed=SPAI,
              #between=.(STAI), observed=STAI,
              detailed=T, type=2) %>% apa::anova_apa(force_sph_corr=T)

eye.diagnosticity.ms.subj = eye.diagnosticity.ms %>% mutate(ms = ms / 1000) %>% #readability of plot -> convert to sec
  group_by(diagnostic, Diagnosticity, ROI, subject) %>% summarise(ms=mean(ms, na.rm=T))
eye.diagnosticity.ms.subj.analysis = eye.diagnosticity.ms.subj %>% filter(subject %in% exclusions.eye.ms == F) #manual exclusion because of extreme latency
print(eye.main.ms <- eye.diagnosticity.ms.subj.analysis %>% summarise(ms.se=se(ms, na.rm=T), ms=mean(ms, na.rm=T)) %>% 
        ggplot(aes(x=diagnostic, y=ms, color=ROI, shape=Diagnosticity)) + 
        #ggplot(aes(x=diagnostic, y=ms, color=as.numeric(ROI)==as.numeric(diagnostic))) + 
        geom_dotplot(data=eye.diagnosticity.ms.subj.analysis %>% filter(ROI=="Eyes"), mapping=aes(group=interaction(ROI, diagnostic), fill=ROI), stackdir="down", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.75) +
        geom_dotplot(data=eye.diagnosticity.ms.subj.analysis %>% filter(ROI!="Eyes"), mapping=aes(group=interaction(ROI, diagnostic), fill=ROI), stackdir="up", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.75) +
        geom_point(size=6, position=dodge) + geom_errorbar(aes(ymin=ms-ms.se*1.96, ymax=ms+ms.se*1.96), size=2, position=dodge) +
        scale_color_discrete(labels=c("Eyes", "Mouth/Nose")) +
        #scale_color_discrete(name="Diagnosticity of ROI", labels=c("Diagnostic", "Non-Diagnostic")) +
        #scale_color_viridis_d(labels=c("Eyes", "Mouth/Nose")) +
        ylab("Latency to ROI (sec)") + xlab("Diagnostic Region") + myGgTheme)

# # Figure in German 
#eye.diagnosticity.ms.subj.analysis <- eye.diagnosticity.ms.subj.analysis %>%
#   mutate(ROI =
#          case_when(
#            ROI == "Eyes"~ "Augen",
#            ROI == "Mouth/Nose" ~ "Mund/Nase"
#          )) %>%
#   mutate(Diagnosticity =
#            case_when(
#              Diagnosticity == "Diagnostic" ~ "Diagnostisch",
#              Diagnosticity == "Non-Diagnostic" ~ "Nicht-Diagnostisch"
#            ))%>%
#   mutate(diagnostic =
#            case_when(
#              diagnostic == "Eyes" ~ "Augen",
#              diagnostic == "Mouth/Nose" ~ "Mund/Nase"
#            ))%>%
#   rename(Diagnostizität = Diagnosticity)
# print(eye.main.ms <- eye.diagnosticity.ms.subj.analysis %>% summarise(ms.se=se(ms, na.rm=T), ms=mean(ms, na.rm=T)) %>%
#         ggplot(aes(x=diagnostic, y=ms, color=ROI, shape=Diagnostizität)) +
#         #ggplot(aes(x=diagnostic, y=ms, color=as.numeric(ROI)==as.numeric(diagnostic))) +
#         geom_dotplot(data=eye.diagnosticity.ms.subj.analysis %>% filter(ROI=="Augen"), mapping=aes(group=interaction(ROI, diagnostic), fill=ROI), stackdir="down", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.75) +
#         geom_dotplot(data=eye.diagnosticity.ms.subj.analysis %>% filter(ROI!="Augen"), mapping=aes(group=interaction(ROI, diagnostic), fill=ROI), stackdir="up", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.75) +
#         geom_point(size=6, position=dodge) + geom_errorbar(aes(ymin=ms-ms.se*1.96, ymax=ms+ms.se*1.96), size=2, position=dodge) +
#         scale_color_discrete(labels=c("Augen", "Mund/Nase")) +
#         #scale_color_discrete(name="Diagnosticity of ROI", labels=c("Diagnostic", "Non-Diagnostic")) +
#         #scale_color_viridis_d(labels=c("Eyes", "Mouth/Nose")) +
#         ylab("Latenz zu ROI (s)") + xlab("Diagnostische Region") + myGgTheme)


eye.diagnosticity.ms.subj %>% mutate(ms.z = scale(ms)[,1]) %>% arrange(desc(abs(ms.z))) #exclude 18 & 62
# eye.diagnosticity.ms.subj %>% filter(Diagnosticity=="Diagnostic") %>% arrange(desc(ms)) %>% mutate(ms.z = scale(ms)[,1]) #exclude subject 18?
# eye.diagnosticity.ms.subj %>% filter(Diagnosticity=="Non-Diagnostic") %>% arrange(ms) %>% mutate(ms.z = scale(ms)[,1])

#ROI x Diagnosticity
eye.diagnosticity.ms %>% group_by(Diagnosticity, ROI) %>% summarise(ms=mean(ms, na.rm=T))
#it takes particularly long to look into non-diagnostic mouth/nose but not for non-diagnostic eyes

#SPAI main effect
eye.diagnosticity.ms.spai = eye.diagnosticity.ms %>% 
  filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
  mutate(ms = ms / 1000) %>% #readability of y-axis
  group_by(SPAI, subject) %>% summarise(ms.se = se(ms, na.rm=T), ms = mean(ms, na.rm=T))
eye.diagnosticity.ms.spai %>% with(cor.test(ms, SPAI, alternative="two.sided")) %>% correlation_out()
eye.diagnosticity.ms.spai %>% 
  ggplot(aes(y=ms, x=SPAI, color=SPAI)) +
  geom_errorbar(aes(ymin=ms-ms.se*1.96, ymax=ms+ms.se*1.96), width=.05) +
  geom_smooth(method="lm", color="black") + geom_point(size=4) +
  ylab("Average Time to ROIs (sec)") +
  scale_color_viridis_c() + myGgTheme + theme(legend.position = "none")

#SPAI x Diagnosticity
eye.diagnosticity.ms.spaiXdia = eye.diagnosticity.ms %>% 
  filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
  mutate(ms = ms / 1000) %>% #readability of y-axis
  group_by(subject, SPAI, Diagnosticity) %>% summarise(ms.se = se(ms, na.rm=T), ms = mean(ms, na.rm=T))
eye.diagnosticity.ms.spaiXdia %>% group_by(Diagnosticity) %>% 
  summarise(cor.test(ms, SPAI, alternative="two.sided") %>% apa::cor_apa(r_ci=T, print=F))
eye.diagnosticity.ms.spaiXdia %>% 
  ggplot(aes(y=ms, x=SPAI, color=SPAI)) +
  facet_wrap(vars(Diagnosticity)) +
  geom_errorbar(aes(ymin=ms-ms.se*1.96, ymax=ms+ms.se*1.96), width=.05) +
  geom_smooth(method="lm", color="black") + geom_point(size=4) +
  ylab("Average Time to ROIs (sec)") +
  scale_color_viridis_c() + myGgTheme + theme(legend.position = "none")

# #SPAI x Diagnosticity in German
# eye.diagnosticity.ms.spaiXdia = eye.diagnosticity.ms %>% 
#   filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
#   mutate(ms = ms / 1000) %>% #readability of y-axis
#   group_by(subject, SPAI, Diagnosticity) %>% summarise(ms.se = se(ms, na.rm=T), ms = mean(ms, na.rm=T))
# eye.diagnosticity.ms.spaiXdia %>% group_by(Diagnosticity) %>% 
#   summarise(cor.test(ms, SPAI, alternative="two.sided") %>% apa::cor_apa(r_ci=T, print=F))
# eye.diagnosticity.ms.spaiXdia <- eye.diagnosticity.ms.spaiXdia %>%
#   mutate(Diagnosticity =
#            case_when(
#              Diagnosticity == "Diagnostic" ~ "Diagnostisch",
#              Diagnosticity == "Non-Diagnostic" ~ "Nicht-Diagnostisch"
#            ))%>%
#   rename(Diagnostizität = Diagnosticity)
# eye.diagnosticity.ms.spaiXdia %>% 
#   ggplot(aes(y=ms, x=SPAI, color=SPAI)) +
#   facet_wrap(vars(Diagnostizität)) +
#   geom_errorbar(aes(ymin=ms-ms.se*1.96, ymax=ms+ms.se*1.96), width=.05) +
#   geom_smooth(method="lm", color="black") + geom_point(size=4) +
#   ylab("Durchschnittliche Latenz zu ROIs (s)") +
#   scale_color_viridis_c() + myGgTheme + theme(legend.position = "none")

#Sqrt of Latency to ROIs -----------------------------------------------------------------------------------------------
#exclusions.eye.ms = c(18, 62) #TODO don't exclude but prune subjects?
exclusions.eye.ms = c() #results without exclusions

eye.diagnosticity.ms.sqrt <- eye.diagnosticity.ms%>%
  mutate(ms = sqrt(ms))

eye.diagnosticity.ms.sqrt %>% filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
  ez::ezANOVA(dv=.(ms), wid=.(subject),
              within=.(ROI, Diagnosticity), within_full=.(threat, trial),
              #within=.(threat, ROI, Diagnosticity),
              between=.(SPAI), observed=SPAI,
              #between=.(STAI), observed=STAI,
              detailed=T, type=2) %>% apa::anova_apa(force_sph_corr=T)

eye.diagnosticity.ms.sqrt.subj = eye.diagnosticity.ms.sqrt %>% mutate(ms = ms / 1000) %>% #readability of plot -> convert to sec
  group_by(diagnostic, Diagnosticity, ROI, subject) %>% summarise(ms=mean(ms, na.rm=T))
eye.diagnosticity.ms.sqrt.subj.analysis = eye.diagnosticity.ms.sqrt.subj %>% filter(subject %in% exclusions.eye.ms == F) #manual exclusion because of extreme latency

print(eye.main.ms <- eye.diagnosticity.ms.sqrt.subj.analysis %>% summarise(ms.se=se(ms, na.rm=T), ms=mean(ms, na.rm=T)) %>% 
        ggplot(aes(x=diagnostic, y=ms, color=ROI, shape=Diagnosticity)) + 
        #ggplot(aes(x=diagnostic, y=ms, color=as.numeric(ROI)==as.numeric(diagnostic))) + 
        geom_dotplot(data=eye.diagnosticity.ms.sqrt.subj.analysis %>% filter(ROI=="Eyes"), mapping=aes(group=interaction(ROI, diagnostic), fill=ROI), stackdir="down", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.75) +
        geom_dotplot(data=eye.diagnosticity.ms.sqrt.subj.analysis %>% filter(ROI!="Eyes"), mapping=aes(group=interaction(ROI, diagnostic), fill=ROI), stackdir="up", binaxis="y", alpha=.25, color="black", stackratio=1, dotsize=.75) +
        geom_point(size=6, position=dodge) + geom_errorbar(aes(ymin=ms-ms.se*1.96, ymax=ms+ms.se*1.96), size=2, position=dodge) +
        scale_color_discrete(labels=c("Eyes", "Mouth/Nose")) +
        #scale_color_discrete(name="Diagnosticity of ROI", labels=c("Diagnostic", "Non-Diagnostic")) +
        #scale_color_viridis_d(labels=c("Eyes", "Mouth/Nose")) +
        ylab(expression("Latency to ROI (" * sqrt(ms) * ")")) + xlab("Diagnostic Region") + myGgTheme)

#SPAI main effect
eye.diagnosticity.ms.sqrt.spai = eye.diagnosticity.ms.sqrt %>% 
  filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
  mutate(ms = ms / 1000) %>% #readability of y-axis
  group_by(SPAI, subject) %>% summarise(ms.se = se(ms, na.rm=T), ms = mean(ms, na.rm=T))
eye.diagnosticity.ms.sqrt.spai %>% with(cor.test(ms, SPAI, alternative="two.sided")) %>% correlation_out()
eye.diagnosticity.ms.sqrt.spai %>% 
  ggplot(aes(y=ms, x=SPAI, color=SPAI)) +
  geom_errorbar(aes(ymin=ms-ms.se*1.96, ymax=ms+ms.se*1.96), width=.05) +
  geom_smooth(method="lm", color="black") + geom_point(size=4) +
  ylab(expression("Latency to ROI (" * sqrt(ms) * ")")) +
  scale_color_viridis_c() + myGgTheme + theme(legend.position = "none")

#SPAI x Diagnosticity
eye.diagnosticity.ms.sqrt.spaiXdia = eye.diagnosticity.ms.sqrt %>% 
  filter(subject %in% exclusions.eye.ms == F) %>% #manual exclusion because of extreme latency
  mutate(ms = ms / 1000) %>% #readability of y-axis
  group_by(subject, SPAI, Diagnosticity) %>% summarise(ms.se = se(ms, na.rm=T), ms = mean(ms, na.rm=T))
eye.diagnosticity.ms.sqrt.spaiXdia %>% group_by(Diagnosticity) %>% 
  summarise(cor.test(ms, SPAI, alternative="two.sided") %>% apa::cor_apa(r_ci=T, print=F))
eye.diagnosticity.ms.sqrt.spaiXdia %>% 
  ggplot(aes(y=ms, x=SPAI, color=SPAI)) +
  facet_wrap(vars(Diagnosticity)) +
  geom_errorbar(aes(ymin=ms-ms.se*1.96, ymax=ms+ms.se*1.96), width=.05) +
  geom_smooth(method="lm", color="black") + geom_point(size=4) +
  ylab(expression("Latency to ROI (" * sqrt(ms) * ")")) +
  scale_color_viridis_c() + myGgTheme + theme(legend.position = "none")


# Exploration -------------------------------------------------------------
eye.anova.big = eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% mutate(bin = as.factor(bin)) %>% 
  ez::ezANOVA(dv=.(relDwell), wid=.(subject), 
              within=.(Diagnosticity, ROI, threat, bin), within_full=.(trial),
              detailed=T, type=2)
eye.anova.big %>% apa::anova_apa(force_sph_corr=T)

# Diagnosticity x ROI x bins (analysis time + average "bin") [n.s. but comparison to first study + summary of 2-way interactions with bin]
eye.diagnosticity.bins.average.ga = eye.diagnosticity.bins.average %>% 
  group_by(diagnostic, Diagnosticity, ROI, bin, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
  summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
  ungroup() %>% mutate(relDwell=relDwell*100, relDwell.se=relDwell.se*100, #improves readability of axes
                       bin = case_when(bin > ratingStart/1000 & Diagnosticity=="Diagnostic" & ROI!="Eyes" ~ bin - .125,
                                       bin > ratingStart/1000 & Diagnosticity!="Diagnostic" & ROI=="Eyes" ~ bin + .125,
                                       TRUE ~ bin)) #create a specific dodge manually
print(eye.time <- eye.diagnosticity.bins.average.ga %>% filter(bin <= ratingStart/1000) %>% 
        #ggplot(aes(x=bin, y=relDwell, color=ROI, shape=diagnostic)) +
        ggplot(aes(x=bin, y=relDwell, color=ROI, shape=Diagnosticity)) + 
        #geom_errorbar(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96), size=2) +
        geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96, fill=ROI), color=NA, alpha=.15) +
        geom_line(size=2) + geom_point(size=6) + 
        geom_vline(xintercept=4.375) + #hacked :/
        geom_errorbar(data=eye.diagnosticity.bins.average.ga %>% filter(bin > ratingStart/1000), 
                      aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96, color=ROI), size=2, width=.4) +
        geom_point(data=eye.diagnosticity.bins.average.ga %>% filter(bin > ratingStart/1000), size=6) +
        scale_color_discrete(labels=c("Eyes", "Mouth/Nose")) + scale_fill_discrete(labels=c("Eyes", "Mouth/Nose")) + 
        #scale_shape_discrete(name="Diagnosticity of ROI", labels=c("Non-Diagnostic", "Diagnostic")) + guides(shape = guide_legend(reverse = TRUE)) +
        ylab("Relative Dwell Time (%)") + xlab("Trial Time (sec)") + #labs(shape="Diagnostic") +
        scale_x_continuous(breaks=c(1:4, 4.75), labels=c(1:4, "Avg")) + myGgTheme) #hacked :/

#plotly::ggplotly(eye.time)

# Diagnosticity x ROI x bins (full trial time)
eye.diagnosticity.bins %>% group_by(diagnostic, Diagnosticity, ROI, bin, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>%
  ungroup() %>% mutate(relDwell = relDwell*100) %>% #improves readability of axes
  group_by(diagnostic, Diagnosticity, ROI, bin) %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>%
  ggplot(aes(x=bin, y=relDwell, color=ROI, shape=Diagnosticity)) +
  geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96, fill=ROI), color=NA, alpha=.15) +
  geom_line(size=2) + geom_point(size=6) + 
  geom_vline(xintercept=(ratingStart + binResolution/2) / 1000) +
  scale_color_discrete(labels=c("Eyes", "Mouth/Nose")) + scale_fill_discrete(labels=c("Eyes", "Mouth/Nose")) +
  ylab("Relative Dwell Time (%)") + xlab("Trial Time (sec)") + myGgTheme #+ labs(shape="Diagnostic"))

# Diagnosticity x ROI x bin x threat (n.s.)
# eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% 
#   group_by(diagnostic, Diagnosticity, ROI, bin, threat, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
#   ungroup() %>% mutate(relDwell = relDwell*100) %>% #improves readability of axes
#   group_by(diagnostic, Diagnosticity, ROI, bin, threat) %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
#   mutate(threat = recode(threat, `1` = "CS-", `2` = "GS1", `3` = "GS2", `4` = "GS3", `5` = "GS4", `6` = "CS+")) %>% 
#   ggplot(aes(x=bin, y=relDwell, color=ROI, shape=Diagnosticity)) + 
#   facet_wrap(vars(threat)) +
#   geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96, fill=ROI), color=NA, alpha=.15) +
#   geom_point(size=6) + geom_line(size=2) +
#   scale_color_discrete(labels=c("Eyes", "Mouth/Nose")) + scale_fill_discrete(labels=c("Eyes", "Mouth/Nose")) + 
#   ylab("Relative Dwell Time (%)") + xlab("Trial Time (sec)") + myGgTheme #+ labs(shape="Diagnostic"))

# Diagnosticity x bin
eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% 
  group_by(Diagnosticity, bin, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
  ungroup() %>% mutate(relDwell = relDwell*100) %>% #improves readability of axes
  group_by(Diagnosticity, bin) %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
  ggplot(aes(x=bin, y=relDwell, shape=Diagnosticity)) + 
  geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96), color=NA, fill="black", alpha=.15) +
  geom_point(size=6) + geom_line(size=2) +
  ylab("Relative Dwell Time (%)") + xlab("Trial Time (sec)") + myGgTheme #+ labs(shape="Diagnostic"))

# ROI x bin
eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% 
  group_by(ROI, bin, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
  ungroup() %>% mutate(relDwell = relDwell*100) %>% #improves readability of axes
  group_by(ROI, bin) %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
  ggplot(aes(x=bin, y=relDwell, color=ROI)) + 
  geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96, fill=ROI), color=NA, alpha=.15) +
  geom_point(size=6) + geom_line(size=2) +
  scale_color_discrete(labels=c("Eyes", "Mouth/Nose")) + scale_fill_discrete(labels=c("Eyes", "Mouth/Nose")) + 
  ylab("Relative Dwell Time (%)") + xlab("Trial Time (sec)") + myGgTheme #+ labs(shape="Diagnostic"))

# bin main effect
eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% 
  group_by(bin, subject) %>% summarise(relDwell.se = se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
  ungroup() %>% mutate(relDwell = relDwell*100*2, relDwell.se = relDwell.se*100*2) %>% #improves readability of axes
  group_by(bin) %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
  ggplot(aes(x=bin, y=relDwell)) + 
  geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96), color=NA, fill="black", alpha=.15) +
  geom_point(size=6) + geom_line(size=2) +
  ylab("Relative Dwell Time on ROIs (%)") + xlab("Trial Time (sec)") + myGgTheme #+ labs(shape="Diagnostic"))

# ROI x threat x bin (n.s.)
# eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% 
#   group_by(ROI, bin, threat, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
#   ungroup() %>% mutate(relDwell = relDwell*100) %>% #improves readability of axes
#   group_by(ROI, bin, threat) %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
#   mutate(threat = recode(threat, `1` = "CS-", `2` = "GS1", `3` = "GS2", `4` = "GS3", `5` = "GS4", `6` = "CS+")) %>% 
#   ggplot(aes(x=bin, y=relDwell, color=ROI)) + 
#   facet_wrap(vars(threat)) +
#   geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96, fill=ROI), color=NA, alpha=.15) +
#   geom_point(size=6) + geom_line(size=2) +
#   scale_color_discrete(labels=c("Eyes", "Mouth/Nose")) + scale_fill_discrete(labels=c("Eyes", "Mouth/Nose")) + 
#   ylab("Relative Dwell Time (%)") + xlab("Trial Time (sec)") + myGgTheme #+ labs(shape="Diagnostic"))


#interactions with threat
# Diagnosticity x threat
eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% 
  group_by(Diagnosticity, threat, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
  ungroup() %>% mutate(relDwell = relDwell*100) %>% #improves readability of axes
  group_by(Diagnosticity, threat) %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
  mutate(threat = recode(threat, `1` = "CS-", `2` = "GS1", `3` = "GS2", `4` = "GS3", `5` = "GS4", `6` = "CS+")) %>% 
  ggplot(aes(x=threat, y=relDwell, shape=Diagnosticity, group=Diagnosticity)) + 
  geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96), color=NA, fill="black", alpha=.15) +
  geom_point(size=6) + geom_line(size=2) +
  ylab("Relative Dwell Time (%)") + xlab("Threat Level") + myGgTheme #+ labs(shape="Diagnostic"))

# ROI x threat
eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% 
  group_by(ROI, threat, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
  ungroup() %>% mutate(relDwell = relDwell*100) %>% #improves readability of axes
  group_by(ROI, threat) %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
  mutate(threat = recode(threat, `1` = "CS-", `2` = "GS1", `3` = "GS2", `4` = "GS3", `5` = "GS4", `6` = "CS+")) %>% 
  ggplot(aes(x=threat, y=relDwell, color=ROI, fill=ROI, group=ROI)) + 
  geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96), color=NA, alpha=.15) +
  geom_point(size=6) + geom_line(size=2) +
  ylab("Relative Dwell Time (%)") + xlab("Threat Level") + myGgTheme #+ labs(shape="Diagnostic"))

#interaction with threat AND bin
# threat x bin (color)
eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% 
  group_by(bin, threat, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
  ungroup() %>% mutate(relDwell = relDwell*100*2) %>% #improves readability of axes
  group_by(bin, threat) %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
  mutate(threat = recode(threat, `1` = "CS-", `2` = "GS1", `3` = "GS2", `4` = "GS3", `5` = "GS4", `6` = "CS+")) %>% 
  ggplot(aes(x=bin, y=relDwell, color=threat, fill=threat, shape=threat, group=threat)) + 
  geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96), color=NA, alpha=.15) +
  geom_line(size=2) + geom_point(size=6) +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors) + 
  labs(color="Threat", fill="Threat", shape="Threat") +
  ylab("Relative Dwell Time on ROIs (%)") + xlab("Trial Time (sec)") + myGgTheme #+ labs(shape="Diagnostic"))

# threat x bin
eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% 
  group_by(bin, threat, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
  ungroup() %>% mutate(relDwell = relDwell*100*2) %>% #improves readability of axes
  group_by(bin, threat) %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
  mutate(threat = recode(threat, `1` = "CS-", `2` = "GS1", `3` = "GS2", `4` = "GS3", `5` = "GS4", `6` = "CS+")) %>% 
  ggplot(aes(x=bin, y=relDwell)) + 
  facet_wrap(vars(threat)) +
  geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96), color=NA, fill="black", alpha=.15) +
  geom_point(size=6) + geom_line(size=2) +
  ylab("Relative Dwell Time on ROIs (%)") + xlab("Trial Time (sec)") + myGgTheme #+ labs(shape="Diagnostic"))

# Diagnosticity x threat x bin
eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% 
  group_by(Diagnosticity, bin, threat, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
  ungroup() %>% mutate(relDwell = relDwell*100) %>% #improves readability of axes
  group_by(Diagnosticity, bin, threat) %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
  mutate(threat = recode(threat, `1` = "CS-", `2` = "GS1", `3` = "GS2", `4` = "GS3", `5` = "GS4", `6` = "CS+")) %>% 
  ggplot(aes(x=bin, y=relDwell, shape=Diagnosticity)) + 
  facet_wrap(vars(threat)) +
  geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96), color=NA, fill="black", alpha=.15) +
  geom_point(size=6) + geom_line(size=2) +
  ylab("Relative Dwell Time (%)") + xlab("Trial Time (sec)") + myGgTheme #+ labs(shape="Diagnostic"))

# ROI x threat x bin
eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% 
  group_by(ROI, bin, threat, subject) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
  ungroup() %>% mutate(relDwell = relDwell*100) %>% #improves readability of axes
  mutate(threat = recode(threat, `1` = "CS-", `2` = "GS1", `3` = "GS2", `4` = "GS3", `5` = "GS4", `6` = "CS+")) %>% 
  group_by(ROI, threat, bin) %>% summarise(relDwell.se=se(relDwell, na.rm=T), relDwell=mean(relDwell, na.rm=T)) %>% 
  ggplot(aes(x=bin, y=relDwell, color=ROI)) + 
  facet_wrap(vars(threat)) +
  geom_ribbon(aes(ymin=relDwell-relDwell.se*1.96, ymax=relDwell+relDwell.se*1.96, fill=ROI), color=NA, alpha=.15) +
  geom_point(size=6) + geom_line(size=2) +
  scale_color_discrete(labels=c("Eyes", "Mouth/Nose")) + scale_fill_discrete(labels=c("Eyes", "Mouth/Nose")) + 
  ylab("Relative Dwell Time (%)") + xlab("Trial Time (sec)") + myGgTheme #+ labs(shape="Diagnostic"))

# diagnosticity x ROI x bin (cp. 2-way interaction in confirmatory ANOVA) (n.s.)
# eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% group_by(diagnostic, bin, ROI) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
#   summarise(relDwell.anyROI = sum(relDwell)) %>% spread(bin, relDwell.anyROI)
# eye.diagnosticity.bins %>% filter(bin <= ratingStart/1000) %>% 
#   group_by(diagnostic, bin, subject, ROI) %>% summarise(relDwell=mean(relDwell, na.rm=T)) %>% 
#   summarise(relDwell.anyROI = sum(relDwell)) %>% #collapse over ROIs (-> sum instead of mean)
#   spread(diagnostic, relDwell.anyROI) %>% 
#   mutate(EyeBias=(Eyes-`Mouth/Nose`)*100) %>% 
#   summarise(EyeBias.se = se(EyeBias, na.rm=T), EyeBias = mean(EyeBias, na.rm=T)) %>% 
#   #select(-EyeBias.se) %>% spread(bin, EyeBias) #table of values
#   ggplot(aes(x=bin, y=EyeBias)) + geom_hline(yintercept=0) + geom_path(size=2) + 
#   geom_ribbon(aes(ymin=EyeBias-EyeBias.se*1.96, ymax=EyeBias+EyeBias.se*1.96), color=NA, fill="black", alpha=.15) +
#   #geom_point(size=6, aes(color=EyeBias>0)) + 
#   geom_point(aes(color=EyeBias), size=6) +
#   ylab("Eye Bias (%)") + xlab("Trial Time (sec)") +
#   scale_color_viridis_c()




# Wide format for correlations --------------------------------------------
eyes.wide = eye.gen.diagnostic %>% select(-contains(".se")) %>% 
  gather(measure, value, -c("subject", "diagnostic")) %>% unite(temp, measure, diagnostic) %>% spread(temp, value)
measures = c("dwell", "dwell.non", "dwell.rois", "ms.diag", "ms.nondiag", "diagFirst", "roiSwitch")
names(eyes.wide) = c("subject", paste0("Gen_", do.call(paste0, expand.grid(c("eyes_", "mn_"), measures %>% sort()))))
# names(eyes.wide)[-1] = names(eyes.wide)[-1] %>% {ifelse(grepl("Eyes", ., fixed=T), "eyes", "mn")} %>% 
#   paste("Gen", ., names(eyes.wide)[-1] %>% gsub("_.*", "", .) %>% gsub(".m", "", ., fixed=T), sep="_")

eyes.wide = eye.gen %>% select(-contains(".se")) %>% merge(eyes.wide, by="subject", all=T) %>% tibble()
names(eyes.wide)[1 + 1:length(measures)] = paste0("Gen_all_", measures)
#write_rds(eyes.wide, "eyes.wide.rds" %>% paste0(path.rds, .))
