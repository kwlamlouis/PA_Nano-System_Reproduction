R.version
citation()
citation("lme4")
citation("tidyverse")
citation("ggridges")
citation("ggplot2")
citation("extrafont")
citation("DHARMa")
citation("emmeans")
citation("glmmTMB")

#Change language to enlgish 
Sys.setlocale("LC_ALL", "English")

library(tidyverse)
library(ggridges)
library(ggplot2)
library(extrafont)
library(lme4)
library(DHARMa)  #testDispersion, simulateResiduals
library(emmeans)
library(glmmTMB) #Beta-Binomial, NB, "ar1" function


####<1. RT: Rank>####
#"Rank" represents one lunar cycle, total of 28 days
Rank_RT <- read_csv("C:/Users/USER/Downloads/0. Louis_Laptop/(2) R/1. R test/1. Nano Project/1. Rank_RT.csv")
Rank_RT$M <-as.factor(Rank_RT$M)
Rank_RT$Rank <- as.numeric(Rank_RT$Rank)
Rank_RT$Lunar_Month <-as.factor(Rank_RT$Lunar_Month)
Rank_RT$lunar_day <-as.numeric(Rank_RT$lunar_day)
Rank_RT$colony_id <-as.factor(Rank_RT$colony_id)
Rank_RT$Tank <- as.factor(Rank_RT$Tank)
Rank_RT$RO <- as.integer(Rank_RT$RO)
head(Rank_RT)

#Find weighted.MLD
#For RT, remove "Temp" first as "Temp" has more NAs than "RO" (in "RO, NA = no collection)
#If not,  "na.omit" code will remove rows that have "RO" but NA in "Temp" -> affecting MLD
Rank_RT.noTemp <-
  Rank_RT%>%
  select(-"Temp")

#1. Find monthly sum RO
RT_sum <-
  Rank_RT.noTemp%>%
  group_by(colony_id, Source, M)%>%
  mutate(month.sum = sum(RO, na.rm=T))

#2. calculate RT weight based on RO/month.sum
RT_weight <- 
  RT_sum%>%
  mutate(RT_weight = (RO/month.sum))

#Remove NA (NAs because no collection)
RT_weight_no.NA <- na.omit(RT_weight)

#3. MLD
#(1) MLD of each colony --> for "1.2 Analysis" LMER analysis
RT.MLD <- 
  RT_weight_no.NA%>%
  group_by(colony_id, Source, Tank, M)%>%
  summarise(W.MLD = weighted.mean(Rank,RT_weight))
print(RT.MLD, n=79)

#(2) MLD of each Source per month --> to add "MLD" in plot manually
RT.MLD2 <- 
  RT_weight_no.NA%>%
  group_by(M, Source)%>%
  summarise(W.MLD = weighted.mean(Rank,RT_weight))


####*1.1 Plot ####
#x = Rank -> then transform "x value" manually back to corresponding lunar day. 
#OL
OL_RT_subset<-
  RT_weight_no.NA%>%
  filter(Source =="OL")

#set reproductive cycles in order
OL_RT_subset$Lunar_Month <- factor(OL_RT_subset$M, 
                                   levels=c("M1", "M2", "M3", "M4", "M5", "M6"))
OL.plot <-
  ggplot(data=OL_RT_subset, aes(y=M, x=Rank)) +
  geom_density_ridges(aes(height=after_stat(density),
                          weight=RT_weight),    
                      scale= 0.95,
                      stat="density", fill="palevioletred1", bw=0.5)+    #specify bandwidth
  scale_x_continuous(limits=c(1, 30), breaks=c(1,5,10,15,20,25))+
  geom_vline(xintercept=5, linetype="dashed", 
             color = "black", linewidth=0.5)+
  theme_classic()+
  xlab("Lunar day")+
  ylab("Reproductive cycle")+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
OL.plot

#IL
IL_RT_subset<-
  RT_weight_no.NA%>%
  filter(Source =="IL")

#set reproductive cycles in order
IL_RT_subset$Lunar_Month <- factor(IL_RT_subset$M, 
                                   levels=c("M1", "M2", "M3", "M4", "M5", "M6"))

IL.plot <-
  ggplot(data=IL_RT_subset, aes(y=M, x=Rank)) +
  geom_density_ridges(aes(height=after_stat(density),
                          weight=RT_weight),    
                      scale= 0.95,
                      stat="density", fill="dodgerblue4", bw=0.5)+
  scale_x_continuous(limits=c(1, 30), breaks=c(1,5,10,15,20,25))+
  geom_vline(xintercept=5, linetype="dashed", 
             color = "black", linewidth=0.5)+
  theme_classic()+
  xlab("Lunar day")+
  ylab("Reproductive cycle")+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
IL.plot


####*1.2 Analysis####
#Linear mixed effect model
WT.MLD.lmer <- lmer(W.MLD~Source*M + (1|Tank) + (1|colony_id), data=RT.MLD)
summary(WT.MLD.lmer)

#check assumptions
plot(fitted(WT.MLD.lmer),residuals(WT.MLD.lmer))
hist(residuals(WT.MLD.lmer))
qqnorm(residuals(WT.MLD.lmer))
simulationOutput.WT.MLD.lmer <- simulateResiduals(fittedModel = WT.MLD.lmer, plot = T)
testDispersion(simulationOutput.WT.MLD.lmer)

#Post hoc
#[1] Within same month, IL vs OL
W.MLD.posthoc <- emmeans(WT.MLD.lmer, pairwise~Source|M)
#[2] Within same source, month vs month
W.MLD.posthoc2 <- emmeans(WT.MLD.lmer, pairwise~M|Source)


####*1.3 Weighted MLD and temp####
#Monthly MLD vs Monthly Temp.
#Getting MLD
#1. Remove NA in "RO"
Rank_RT.noNA <- 
  Rank_RT%>%
  filter(!is.na(RO))

#2. Find monthly sum RO
RT_sum.2 <-
  Rank_RT.noNA%>%
  group_by(colony_id, Source, M)%>%
  mutate(month.sum = sum(RO, na.rm=T))

#3. calculate RT weight based on RO/month.sum
RT_weight.2 <- 
  RT_sum.2%>%
  mutate(RT_weight = (RO/month.sum))

#4. MLD and monthly temp
RT.MLD.2 <- 
  RT_weight.2%>%
  group_by(colony_id, Source, Tank, M)%>%
  summarise(MM = mean(Temp), W.MLD = weighted.mean(Rank,RT_weight))

#Subset: Outlet each colony
OL_MLD_subset <-
  RT.MLD.2 %>%
  filter(Source =="OL")

#Inlet: each colony
IL_MLD_subset <-
  RT.MLD.2 %>%
  filter(Source =="IL")

####*1.3.1 Plot####
#Correlation: to see if there's a linear relationship

#OL
OL_MLD.Temp.plot <- 
  ggplot(OL_MLD_subset, aes(MM, W.MLD)) + 
  geom_point()+
  geom_smooth()+
  scale_y_continuous(limits=c(0, 25), breaks=c(5,10,15,20, 25))+
  labs(x = "Temperature (¢J) ", y="Mean lunar day (MLD)")+
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
OL_MLD.Temp.plot
#Spearman 
cor.test(OL_MLD_subset$MM, OL_MLD_subset$W.MLD,  method = "spearman",exact=FALSE)

#IL
IL_MLD.Temp.plot <- 
  ggplot(IL_MLD_subset, aes(MM,W.MLD)) + 
  geom_point()+
  geom_smooth()+
  labs(x = "Temperature (¢J) ", y="Mean lunar day (MLD)")+
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
IL_MLD.Temp.plot
#Spearman
cor.test(IL_MLD_subset$MM, IL_MLD_subset$W.MLD,  method = "spearman",exact=FALSE)



####<1. Autocorrelation>####
library(purrr)
library(pracma)
library(dtplyr)
library(minpack.lm)

#Set working directory
setwd("C:\\Users\\USER\\Downloads\\0. Louis_Laptop\\(2) R\\1. R test\\1. Nano Project")
acf_data <- read_csv("1.1 RT_Autocorrelation.csv")

acf_data$M <- as.factor(acf_data$M )
acf_data$colony_id <-as.factor(acf_data$colony_id)
acf_data$RO <- as.integer(acf_data$RO)
head(acf_data)
str(acf_data)

#Compute ACF per colony
#[The autocorrelation function (ACF) gives you a sense of repeating patterns by 
#showing how strongly today's larval count is correlated with the count from X days ago.]

acf_results <- 
  acf_data %>%
  group_by(colony_id) %>%
  group_split()%>%
  map(~ {
    x <- .x$RO  # larval count time series
    acf_res <- acf(x, lag.max = 100, plot = FALSE)
    data.frame(
      Colony = unique(.x$colony_id),
      Lag = acf_res$lag[, 1, 1],
      ACF = acf_res$acf[, 1, 1]
    )%>%
      #exclude lag 0
      filter(Lag != 0)
  }) %>%
  bind_rows()


#____________________________________________
#Example: print one group before fitting for checking
acf_results %>%
  group_by(Colony) %>%
  group_split() %>%
  .[[1]]
#x = lag, y = acf for next session
#___________________________________________

#Define fitting function for Damped Sine Model
#df: referring to just one colony's lag-ACF data at a time

fit_damped_sine <- function(df) {
  df <- df %>%
    rename(x = Lag, y = ACF) %>%
    mutate(x = as.numeric(x), y = as.numeric(y))
  
  fit <- tryCatch({
    nlsLM(
      y ~ a * sin(b * x + d) * exp(-c * x),
      data = df,
      start = list(a = 1, b = 2*pi/28, c = 0.01, d = 0),
      control = nls.lm.control(maxiter = 500)
    )
  }, error = function(e) NULL)
  
  if (!is.null(fit)) {                     #From here, this code is to create data frame for plotting 
    df$fit <- predict(fit, newdata = df)
    df$Colony <- unique(df$Colony)  # keep this for plotting
    return(df)
  } else {
    return(NULL)  # still safe for map_df
  }
}


#Fit the model to each colony: taking the ACFs for each colony, and sending each colony's lag-ACF subset (df) into fit_damped_sine
fitted_results <- acf_results %>%
  group_by(Colony) %>%
  group_split() %>%
  map_df(fit_damped_sine)


#Plot results
ggplot(fitted_results, aes(x = x)) +
  geom_point(aes(y = y)) +
  geom_line(aes(y = fit), color = "blue") +
  facet_wrap(~ Colony, scales = "free_y") +
  labs(x = "Lag (days)", y = "ACF", title = "Damped sine fit to ACF per colony") +
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))

#____________________________________________________
#I6 is missing
#Why? Their  ACF values are too flat, erratic, or don¡¦t resemble a sine-like pattern
#____________________________________________________
#If want to check which one was dropped
all_colonies <- unique(acf_results$Colony)
fitted_colonies <- unique(fitted_results$Colony)
missing_colonies <- setdiff(all_colonies, fitted_colonies)
missing_colonies
#[1] "I7" "O6" "O7"
#____________________________________________________

##Now, to extract a,b,c,d need to change the code a bit
fit_damped_sine2 <- function(df) {
  df <- df %>%
    rename(x = Lag, y = ACF) %>%
    mutate(x = as.numeric(x), y = as.numeric(y))
  
  fit <- tryCatch({
    nlsLM(
      y ~ a * sin(b * x + d) * exp(-c * x),
      data = df,
      start = list(a = 1, b = 2*pi/28, c = 0.01, d = 0),
      control = nls.lm.control(maxiter = 500)
    )
  }, error = function(e) NULL)
  
  list(Colony = unique(df$Colony), fit = fit)   #Changed to this
}

#
fit_results2 <- acf_results %>%
  group_by(Colony) %>%
  group_split() %>%
  map(fit_damped_sine2)


#
extract_params <- function(result) {
  Colony <- result$Colony
  fit <- result$fit
  
  if (is.null(fit)) {
    tibble(Colony = Colony, a = NA, b = NA, c = NA, d = NA, period = NA)
  } else {
    coef_vals <- coef(fit)
    tibble(
      Colony = Colony,
      a = coef_vals["a"],
      b = coef_vals["b"],
      c = coef_vals["c"],
      d = coef_vals["d"],
      period = 2 * pi / coef_vals["b"]
    )
  }
}

param_df <- map_dfr(fit_results2, extract_params)

#Note: trust the plot first, and use the model to summarize or group patterns, not to replace the visual check.


####Extract lag day for each colony to know at which lag (in days) each colony shows repeated peak reproduction
#1. split by colony
colony_list <- split(acf_data, acf_data$colony_id)

#2. Get peak
get_all_peak_lags <- function(df, lag.max = 100, threshold = 0.1) {
  RO <- df %>% arrange(acf_rank) %>% pull(RO)
  if (length(RO) < 3 || all(is.na(RO)) || sum(RO, na.rm = TRUE) == 0) return(NA)
  
  acf_res <- acf(RO, lag.max = lag.max, plot = FALSE, na.action = na.pass)
  acf_vals <- as.vector(acf_res$acf)[-1]  # skip lag 0
  lags <- as.vector(acf_res$lag)[-1]
  
  # Find local maxima above threshold
  peak_indices <- which(
    acf_vals > threshold &
      acf_vals > dplyr::lag(acf_vals) &
      acf_vals > dplyr::lead(acf_vals)
  )
  
  if (length(peak_indices) == 0) return(NA)
  
  peak_lags <- lags[peak_indices]
  return(paste(round(peak_lags), collapse = ","))
}

# Run across all colonies
acf_peaks_df <- map_df(names(colony_list), function(colony_id) {
  lag_str <- get_all_peak_lags(colony_list[[colony_id]])
  tibble(Colony = colony_id, Peak_Lags = lag_str)
})

#Save as csv
write.csv(acf_peaks_df, "acf_peak_lags_by_colony.csv", row.names = FALSE)


#Plot ACF for each reproductive cycle
#-------------------------#
acf_data %>%
  group_by(colony_id) %>%
  group_split() %>%
  walk(~{
    acf(.x$RO, lag.max = 100, main = paste("ACF -", unique(.x$colony_id)))
  })
#-------------------------#
acf_data_I6<-
  acf_data%>%
  filter(colony_id=="I6")

acf(acf_data_I6) 



####<2. RO>####
#Monthly RO per colony
RO <- read_csv("C:/Users/USER/Downloads/0. Louis_Laptop/(2) R/1. R test/1. Nano Project/2. RO.csv")
summary(RO)

RO$Lunar_Month_Year <- factor(RO$Lunar_Month_Year, levels=c("Oct_2023", "Nov_2023", "Dec_2023", "Jan_2024", "Mar_2024"))
RO$Measured_month <-as.factor(RO$Measured_month)
RO$Source <- as.factor(RO$Source)
RO$Tank_ID <- as.factor(RO$Tank_ID)
RO$ID <-as.factor(RO$ID)
RO$RO<-as.integer(RO$RO)
head(RO)


####*2.1 Plot####
system_boxplot <- 
  ggplot(data=RO, aes(x=Measured_month, y=RO, fill=Source)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), 
             size=1.5) +
  scale_fill_manual(values = c("skyblue","pink"))+
  labs(y = "Number of planulae (ind.)", x="Reproductive cycle")+
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
system_boxplot


####*2.2 Analysis####
RO$Lunar_Month_Year <- as.factor(RO$Lunar_Month_Year)
head(RO)

#take a look at reproductive output
plot(RO$RO)
hist(RO$RO)

#m1 (Poisson mixed-effects model)
m1 <- glmmTMB(RO ~ Source*Measured_month+ 
                (1|Tank_ID) + 
                ar1(Measured_month + 0|ID),
              data = RO, 
              family = poisson)

#check assumptions
plot(fitted(m1),residuals(m1))
hist(residuals(m1))
qqnorm(residuals(m1))
simulationOutput.m1 <- simulateResiduals(fittedModel = m1, plot = T)
testDispersion (simulationOutput.m1)

#Posthoc
#[1] Between populations within each cycle
RO_posthoc1 <- emmeans(m1, pairwise~Source|Measured_month)
#[2] Within each population across cycles
RO_posthoc2 <- emmeans(m1, pairwise~Measured_month|Source)

####<2. RO & Temperature>####
#Daily RO per source vs Daily Temp.
#Data
RO_Temp <- read_csv("C:/Users/USER/Downloads/0. Louis_Laptop/(2) R/1. R test/1. Nano Project/2. RO_Temp.csv")

RO_Temp$Source<- as.factor(RO_Temp$Source)
RO_Temp$Lunar<- as.factor(RO_Temp$Lunar)
RO_Temp$RO<- as.integer(RO_Temp$RO)
head(RO_Temp)


####*2.1 Plot####
#Nonlinear plot
ROvsTemp.plot <- 
  ggplot(RO_Temp, aes(Temp,RO)) + 
  geom_point(aes(color=Source))+
  scale_color_manual(values = c("IL"="skyblue", "OL"="pink"))+
  geom_smooth()+
  labs(x = "Temperature (¢J) ", y="Number of planulae (ind.)")+
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right", axis.title.x = element_blank())+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
ROvsTemp.plot


####*2.2 Analysis####
#Spearman
cor.test(RO_Temp$RO, RO_Temp$Temp,  method = "spearman",exact=FALSE)


####*2.3 <Which temperature tx had the highest RO>####
#Grouping "Temperature" into 6 groups from "2. RO_Temp" data
#Group 1: 22~23.9
#Group 2: 24~24.9
#Group 3: 25~25.9
#Group 4: 26~26.9
#Group 5: 27~28.5
Temp2 <-  read_csv("C:/Users/USER/Downloads/0. Louis_Laptop/(2) R/1. R test/1. Nano Project/2. GroupedTemp_RO.csv")

Temp2$Group <- as.factor(Temp2$Group)
Temp2$Source <- as.factor(Temp2$Source)
Temp2$Group2<-as.factor(Temp2$Group2)


####*2.3.1 Plot####
Temp_RO2 <- 
  ggplot(data=Temp2, aes(x=Group, y=RO, group=Group)) + 
  geom_boxplot(fill="grey")+
  geom_point(aes(color=Source))+
  scale_color_manual(values = c("IL"="skyblue", "OL"="pink"))+
  labs(y = "Number of planulae (ind.)", x="")+
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
Temp_RO2

####*2.3.2 Analysis####
#Kruskal-Wallis Test
kruskal.test(RO ~ Group2, data = Temp2)



####<3. Size>####
Size  <- read_csv("C:/Users/USER/Downloads/0. Louis_Laptop/(2) R/1. R test/1. Nano Project/3. Size.csv")

Size$Lunar <- as.Date(Size$Lunar, format ='%Y/%m/%d')
Size$Month <- as.factor(Size$Month)
Size$Site <- as.factor(Size$Site)
Size$ID <- as.factor(Size$ID)
Size$Tank <- as.factor(Size$Tank)
Size$RO <-as.integer(Size$RO)
head(Size)

Summary.size <-
  Size %>%
  group_by(Lunar, Month, Site)%>%
  summarise(mean=mean(TLE), sd=sd(TLE), n=n(), se =sd/sqrt(n) )
Summary.size


####*3.1 Plot####
Size_lineplot <-
  ggplot(Summary.size, aes(x=Lunar, y=mean, color=Site)) + 
  geom_line(size=1)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5, size=0.7)+
  ylim(0,300)+
  scale_color_manual(values = c("skyblue","pink"))+
  labs(x = "Lunar month", y="Total linear extension (TLE, mm)")+
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right")+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
Size_lineplot


####*3.2 Analysis####
#Linear mixed effect  model
Adult.TLE.lmer <- lmer(TLE~Site*Month + (1|Tank) + (1|ID), data=Size)

#Check assumptions
plot(fitted(Adult.fvfm.lmer),residuals(Adult.fvfm.lmer))
hist(residuals(Adult.fvfm.lmer))
qqnorm(residuals(Adult.fvfm.lmer))
testDispersion(Adult.TLE.lmer)
simulationOutput.Frag.TLE <- simulateResiduals(fittedModel = Adult.TLE.lmer, plot = T)

#Posthoc
Adult.TLE.lmer.contrast1 <- emmeans(Adult.TLE.lmer, pairwise~Site|Month)



####<3 Size vs RO>####
####*3.1 Plot#### 
TLE_RO.plot <- 
  ggplot(Size, aes(TLE,RO)) + 
  geom_point(aes(color=Site)) + 
  geom_smooth()+
  scale_color_manual(values = c("OL" = "pink", "Inlet" = "skyblue"))+
  labs(x = "Total linear extension (TLE, mm)", y="Number of planulae (ind.)")+
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right")+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
TLE_RO.plot

####*3.2 Analysis####
#Spearman
cor.test(Size$TLE, Size$RO,  method = "spearman",exact=FALSE)



####<4. FvFm>####
#adult
FVFM <- read_csv("C:/Users/USER/Downloads/0. Louis_Laptop/(2) R/1. R test/1. Nano Project/4. FVFM.csv")

FVFM$Lunar <- as.Date(FVFM$Lunar, format ='%Y/%m/%d')
FVFM$Month <- as.factor(FVFM$Month)
FVFM$Site <- as.factor(FVFM$Site)
FVFM$ID <- as.factor(FVFM$ID)
FVFM$Tank <- as.factor(FVFM$Tank)
FVFM$Stage <- as.factor(FVFM$Stage)
head(FVFM)


####*4.1 Plot####
#Calculating Average, SD, SE
Summary <-
  FVFM %>%
  group_by(Lunar, Site)%>%
  summarise(mean=mean(fvfm), sd=sd(fvfm), n=n(), se =sd/sqrt(n) )
Summary

#Plot
fvfm_lineplot <-
  ggplot(Summary, aes(x=Lunar, y=mean, color=Site)) + 
  geom_line(size=0.8)+
  geom_point(size=2.2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  scale_color_manual(values = c("skyblue","pink"))+
  ylim(0.5,1)+
  labs(y = "Fv/Fm", x="Lunar month")+
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
fvfm_lineplot


####*4.2 Analysis####
#Linear mixed effect model
Adult.fvfm.lmer <- lmer(fvfm~Site*Month + (1|Tank) + (1|ID), data=FVFM)

#Check assumptions
plot(fitted(Adult.fvfm.lmer),residuals(Adult.fvfm.lmer))
hist(residuals(Adult.fvfm.lmer))
qqnorm(residuals(Adult.fvfm.lmer))
testDispersion(Adult.fvfm.lmer)
simulationOutput.Adult.fvfm <- simulateResiduals(fittedModel = Adult.fvfm.lmer, plot = T)

#Posthoc
Contrast.fvfm <- emmeans(Adult.fvfm.lmer, pairwise~Site|Month)




####<5. Temp>####
Temp <- read_csv("0. Louis_Laptop/(2) R/1. R test/1. Nano Project/5. Temp.csv")
Temp$Date <- as.Date(Temp$Date, format='%Y/%m/%d')
Temp$Tank <- as.factor(Temp$Tank)
head(Temp)

ggplot(Temp, aes(x=Date, y=Daily_mean))+
  geom_line()+
  theme_bw()+
  facet_grid(Tank~.)+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right", axis.title.x = element_blank())+
  labs(y="Temperature (¢XC)", x="")

####<6. Water chem>####
Water_chem <- read_csv("0. Louis_Laptop/(2) R/1. R test/1. Nano Project/6. Water chem.csv")
Water_chem$Date <- as.factor(Water_chem$Date)
Water_chem$Tank <- as.factor(Water_chem$Tank)
head(Water_chem)

#alk
ggplot(Water_chem, aes(x=Date, y=alk, group=Tank))+
  geom_line()+
  theme_bw()+
  facet_grid(Tank~.)+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right", axis.title.x = element_blank())+
  labs(y="", x="")+
  ylim(0,7)

#pH
ggplot(Water_chem, aes(x=Date, y=pH, group=Tank))+
  geom_line()+
  theme_bw()+
  facet_grid(Tank~.)+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right", axis.title.x = element_blank())+
  labs(y="", x="")+
  ylim(0,9)

#Ca
ggplot(Water_chem, aes(x=Date, y=Ca, group=Tank))+
  geom_line()+
  theme_bw()+
  facet_grid(Tank~.)+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right", axis.title.x = element_blank())+
  labs(y="", x="")+
  ylim(100,500)

#Mg
ggplot(Water_chem, aes(x=Date, y=Mg, group=Tank))+
  geom_line()+
  theme_bw()+
  facet_grid(Tank~.)+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right", axis.title.x = element_blank())+
  labs(y="", x="")+
  ylim(1000,1500)

#PHOS
ggplot(Water_chem, aes(x=Date, y=PHOS, group=Tank))+
  geom_line()+
  theme_bw()+
  facet_grid(Tank~.)+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right", axis.title.x = element_blank())+
  labs(y="", x="")+
  ylim(0,3)

#Amonia
ggplot(Water_chem, aes(x=Date, y=AMMO, group=Tank))+
  geom_line()+
  theme_bw()+
  facet_grid(Tank~.)+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right", axis.title.x = element_blank())+
  labs(y="", x="")+
  ylim(0,3)

#Nitrite
ggplot(Water_chem, aes(x=Date, y=NITRITE, group=Tank))+
  geom_line()+
  theme_bw()+
  facet_grid(Tank~.)+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right", axis.title.x = element_blank())+
  labs(y="", x="")+
  ylim(0,3)

#Nitrate
ggplot(Water_chem, aes(x=Date, y=NITRATE, group=Tank))+
  geom_line()+
  theme_bw()+
  facet_grid(Tank~.)+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right", axis.title.x = element_blank())+
  labs(y="", x="")+
  ylim(0,3)

####________####
####Save as HD####
setwd("C:/Users/USER/Downloads")
png("acf.png", units="mm", width=174, height=190, res=1200)
ggplot(fitted_results, aes(x = x)) +
  geom_point(aes(y = y)) +
  geom_line(aes(y = fit), color = "blue") +
  facet_wrap(~ Colony, scales = "free_y") +
  labs(x = "Lag (days)", y = "ACF") +
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"))
dev.off()
