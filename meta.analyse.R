# Meta-Analysis
# get HR, SCL, and saccade length in response to threat (6) and neutral-slightly aversive (1-5) stimluli
library(tidyverse)

# Determine Median for Median Split---------------------------------------------

data.wide %>% summarise(
  median = median(SPAI)
)

# HR ---------------------------------------------------------------------------

heart = read_rds("heart.rds" %>% paste0(path.rds, .))

heart.meta <- heart %>%
  mutate(hr = HRchange + hrbl,
         time = as.numeric(time),
         SAD = ifelse(SPAI > 2.61, 1,0),
         condition = ifelse(threat == 6, 1,0)) %>%
  filter(time <= 4, 
         trial >= 20)%>%
  group_by(subject, SPAI, SAD, condition)%>%
  summarise(
    hr = mean(hr, na.rm = T)
  )

means.hr <- heart.meta %>%
  group_by(SAD, condition)%>%
  summarise(
    n = n(),
    hr_mean = mean(hr),
    sd = sd(hr)
  )

# SCL --------------------------------------------------------------------------

eda.df$subject <- as.numeric(eda.df$subject)
eda.meta <- inner_join(eda.df, heart.meta%>%select(subject, SPAI), by = "subject", keep = NULL)

eda.meta <- eda.meta %>%
  mutate(SAD = ifelse(SPAI > 2.61, 1,0),
         condition = ifelse(threat == 6, 1,0)) %>%
  filter(trial >= 20)%>%
  group_by(subject, SAD, condition)%>%
  summarise(
    cr = mean(ln_cr, na.rm = T)
  )

means.eda <- eda.meta %>%
  group_by(SAD, condition)%>%
  summarise(
    n = n(),
    cr_mean = mean(cr),
    sd = sd(cr)
  )

# saccade length ---------------------------------------------------------------

eye.meta <- eye %>%
  mutate(SAD = ifelse(SPAI > 2.61, 1,0),
         condition = ifelse(threat == 6, 1,0)) %>%
  filter(trial >= 20)%>%
  group_by(subject, SAD, condition)%>%
  summarise(
    scan = mean(scanPath, na.rm = T)
  )

means.eye <- eye.meta %>%
  group_by(SAD, condition)%>%
  summarise(
    n = n(),
    scan_mean = mean(scan),
    sd = sd(scan)
  )

