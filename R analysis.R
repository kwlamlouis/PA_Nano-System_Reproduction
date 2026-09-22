R.version
citation()
citation("tidyverse")
citation("ggridges")
citation("ggplot2")
citation("extrafont")
citation("lme4")
citation("lmerTest")
citation("DHARMa")
citation("emmeans")
citation("glmmTMB")
citation("car")
citation("mgcv")
#Change language to enlgish 
Sys.setlocale("LC_ALL", "English")

library(tidyverse)
library(ggridges)
library(ggplot2)
library(extrafont)
windowsFonts(Arial = windowsFont("Arial"))
library(lme4)
library(lmerTest) #to get p values from "summary(lm_model)", " ANOVA III with Satterthwaite’s approximation for degrees of freedom""
library(DHARMa)  #testDispersion, simulateResiduals
library(emmeans)
library(glmmTMB) #Beta-Binomial, NB, "ar1" function
library(car) #ANOVA III with Wald chi-square tests
library(mgcv) #gam

####<1. RT: Rank>####
#"Rank" represents one lunar cycle, total of 28 days
Rank_RT  <- read_csv("C:/Users/lamkw/Downloads/3. R test/Nano Project (Louis Laptop)/1. Rank_RT.csv")

Rank_RT$M <-as.factor(Rank_RT$M)
Rank_RT$Rank <- as.numeric(Rank_RT$Rank)
Rank_RT$Lunar_Month <-as.factor(Rank_RT$Lunar_Month)
Rank_RT$lunar_day <-as.numeric(Rank_RT$lunar_day)
Rank_RT$Source <- as.factor(Rank_RT$Source)
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


####*1.1.1 Plot ####
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
                      stat="density", fill="#D55E00", bw=0.5)+    #specify bandwidth
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
                      stat="density", fill="#0072B2", bw=0.5)+
  scale_x_continuous(limits=c(1, 30), breaks=c(1,5,10,15,20,25))+
  scale_y_discrete(limits = rev)+
    geom_vline(xintercept=5, linetype="dashed", 
             color = "black", linewidth=0.5)+
  theme_classic()+
  xlab("Lunar day")+
  ylab("Reproductive cycle")+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
IL.plot


####*1.1.2 Analysis####
#Linear mixed effect model
WT.MLD.lmer <- lmer(W.MLD~Source*M + (1|Tank) + (1|colony_id), data=RT.MLD)

#check assumptions
plot(fitted(WT.MLD.lmer),residuals(WT.MLD.lmer))
hist(residuals(WT.MLD.lmer))
qqnorm(residuals(WT.MLD.lmer))
simulationOutput.WT.MLD.lmer <- simulateResiduals(fittedModel = WT.MLD.lmer, plot = T)
testDispersion(simulationOutput.WT.MLD.lmer)

summary(WT.MLD.lmer)
anova(WT.MLD.lmer)

#Post hoc
#[1] Within same month, IL vs OL
W.MLD.posthoc <- emmeans(WT.MLD.lmer, pairwise~Source|M)
#[2] Within same source, month vs month
W.MLD.posthoc2 <- emmeans(WT.MLD.lmer, pairwise~M|Source)


####<1.2 W.MLD against Temp>####
#Individual MLD vs corresponding mean monthly temp.
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
  summarise(
    MM = mean(Temp, na.rm = TRUE), 
    W.MLD = weighted.mean(Rank, RT_weight, na.rm = TRUE)
  )



####*1.2.1 Plot####
# 1. GENERATE THE PREDICTED VALUES FROM THE MODEL
# create a smooth grid of temperatures from 23.5°C to 29.5°C for both populations
temp_grid <- seq(23.5, 29.5, length.out = 100)

# Create a data frame for predictions
predict_df <- data.frame(
  MM = rep(temp_grid, 2),
  Source = rep(c("IL", "OL"), each = 100)
)

# Calculate the predicted W.MLD using the fixed effects estimates
# Formula: W.MLD = Intercept + (Beta1 * MM) + (Beta2 * MM^2) + (Beta3 if Source == OL)
intercept <- -413.0212
beta1     <- 31.3424
beta2     <- -0.5783
beta_ol   <- -1.2729


# Calculate predictions on the model's true internal sequential day scale (1-28)
predict_df$pred_seq_day <- with(predict_df, 
                                intercept + (beta1 * MM) + (beta2 * (MM^2)) + ifelse(Source == "OL", beta_ol, 0)
)


# 2. CREATE A TRANSLATION FUNCTION FOR THE PLOT LABELS
# This function dynamically changes sequential numbers (1-28) into the lunar days (26-30, 1-23)
lunar_label_converter <- function(seq_days) {
  onset_lunar_day <- 26
  actual_lunar <- (seq_days - 1 + onset_lunar_day)
  # Wrap around if it exceeds lunar day 30
  actual_lunar <- ifelse(actual_lunar > 30, actual_lunar - 30, actual_lunar)
  return(as.character(actual_lunar))
}


# 3. GENERATE THE VISUALLY SEAMLESS PLOT
ggplot() +
  geom_point(data = RT.MLD.2, aes(x = MM, y = W.MLD, color = Source, shape = Source), 
             alpha = 0.5, size = 2.5, position = position_jitter(width = 0.05, height = 0)) +
  
  geom_line(data = predict_df, aes(x = MM, y = pred_seq_day, color = Source), 
            linewidth = 1.2) +
  theme_classic() +
  scale_color_manual(values = c("IL" = "#0072B2", "OL" = "#D55E00")) +
  scale_shape_manual(values = c("IL" = 16, "OL" = 17)) +
  scale_y_continuous(
    limits = c(1, 25), 
    breaks = seq(1, 25, by = 5),
    labels = lunar_label_converter
  ) +
  labs(
    x = "Mean monthly temperature (°C)",
    y = "Lunar Day",
    color = "Coral source",
    shape = "Coral source"
  ) +
  theme(
    text = element_text(size = 12, family = "Arial"),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )


#### *1.2.2 Analysis####
head(RT.MLD.2)
RT.MLD.2$Source<-as.factor(RT.MLD.2$Source)

MLD_temp_lmm <- lmer(
  W.MLD ~ poly(MM, 2, raw = TRUE) * Source + (1 | Tank) + (1 | colony_id),
  data = RT.MLD.2,
  REML = TRUE
)

MLD_temp_lmm_reduced <- lmer(
  W.MLD ~ poly(MM, 2, raw = TRUE) + Source + (1 | Tank) + (1 | colony_id),
  data = RT.MLD.2, 
  REML = TRUE
)

anova(MLD_temp_lmm, MLD_temp_lmm_reduced)   
#                     npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
#MLD_temp_lmm_reduced    7 361.94 378.53 -173.97    347.94                        -> use this one
#MLD_temp_lmm            9 361.47 382.79 -171.73    343.47 4.4741  2     0.1068
# => identical AIC, but smaller BIC for "MLD_temp_lmm_reduced"

sim_res <- simulateResiduals(MLD_temp_lmm_reduced, plot = TRUE)
anova(MLD_temp_lmm_reduced, type = 3)   
#Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#poly(MM, 2, raw = TRUE) 177.668  88.834     2 60.881 18.9851 3.911e-07 ***
#Source                   26.386  26.386     1 13.297  5.6391   0.03325 * 

summary(MLD_temp_lmm_reduced)
#Fixed effects:
#                          Estimate Std. Error        df t value Pr(>|t|)   
#poly(MM, 2, raw = TRUE)1   31.3424     9.8164   61.0424   3.193  0.00223 **
#poly(MM, 2, raw = TRUE)2   -0.5783     0.1878   61.0489  -3.079  0.00311 **
#SourceOL                   -1.2729     0.5360   13.2971  -2.375  0.03325 * 


####<1.3 Autocorrelation>####
library(purrr)
library(pracma)
library(dtplyr)
library(minpack.lm)
acf_data  <- read_csv("C:/Users/lamkw/Downloads/3. R test/Nano Project (Louis Laptop)/1.1 RT_Autocorrelation.csv")

acf_data$M <- as.factor(acf_data$M )
acf_data$colony_id <-as.factor(acf_data$colony_id)
acf_data$RO <- as.integer(acf_data$RO)

#--------------------------------------------------#
# 1. Calculate ACF for each colony
#--------------------------------------------------#

acf_results <- acf_data %>%
  group_by(colony_id) %>%
  group_modify(~ {
    
    x <- .x %>%
      arrange(acf_rank) %>%
      pull(RO)
    
    acf_res <- acf(
      x,
      lag.max = 100,
      plot = FALSE,
      na.action = na.pass
    )
    
    tibble(
      Colony = unique(.y$colony_id),
      Lag = as.vector(acf_res$lag)[-1],
      ACF = as.vector(acf_res$acf)[-1]
    )
    
  }) %>%
  ungroup()


#--------------------------------------------------#
# 2. Identify empirical ACF peaks
#--------------------------------------------------#

get_acf_peaks <- function(df, threshold = 0.1) {
  
  peak_idx <- which(
    df$ACF > threshold &
      df$ACF > lag(df$ACF) &
      df$ACF > lead(df$ACF)
  )
  
  if (length(peak_idx) == 0) {
    return(tibble(
      Peak_Lag = NA_real_,
      Peak_ACF = NA_real_
    ))
  }
  
  tibble(
    Peak_Lag = df$Lag[peak_idx],
    Peak_ACF = df$ACF[peak_idx]
  )
}

acf_peak_table <- acf_results %>%
  group_by(Colony) %>%
  group_modify(~ get_acf_peaks(.x)) %>%
  ungroup()


#--------------------------------------------------#
# 3. Create one summary row per colony
#--------------------------------------------------#

acf_peak_summary <- acf_peak_table %>%
  group_by(Colony) %>%
  summarise(
    Peak_Lags = if (all(is.na(Peak_Lag))) {
      NA_character_
    } else {
      paste(round(Peak_Lag[!is.na(Peak_Lag)]), collapse = ", ")
    },
    
    Peak_ACF = if (all(is.na(Peak_ACF))) {
      NA_character_
    } else {
      paste(
        round(Peak_ACF[!is.na(Peak_ACF)], 3),
        collapse = ", "
      )
    },
    
    .groups = "drop"
  )

# SAVE
write_csv(acf_peak_summary,"acf_peak_summary.csv")


#--------------------------------------------------#
# 4. Fit damped sine curve for Figure 3
#--------------------------------------------------#

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

fitted_results <- acf_results %>%
  group_by(Colony) %>%
  group_split() %>%
  map_df(fit_damped_sine)


#Plot results
ggplot(fitted_results, aes(x = Lag)) +
  geom_point(aes(y = ACF)) +
  geom_line(aes(y = fit), color = "blue") +
  facet_wrap(~ Colony, scales = "fixed") +
  labs(x = "Lag (days)", y = "ACF", title = "Damped sine fit to ACF per colony") +
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))

write.csv(model_summary_df, "acf_model_summary.csv", row.names = FALSE)


#### <1.4 RO vs Temp>####
# Monthly sum reproductive output (RO)/colony and MEAN temperature (Temp) per tank
ROvsTemp_monthly <- Rank_RT %>%
  group_by(M, Tank, Source, colony_id) %>%
  summarise(
    Monthly_Sum_RO = sum(RO, na.rm = TRUE),
    MM_Temp        = mean(Temp, na.rm = TRUE), # Skips daily NAs to calculate monthly mean
    .groups = "drop"
  )

#### 1.4.1 Plot####
## 1. GENERATE PREDICTIONS FROM THE FITTED GAMM

# Create a smooth temperature grid (100 points) for both IL and OL
temp_seq <- seq(
  min(ROvsTemp_monthly$MM_Temp, na.rm = TRUE),
  max(ROvsTemp_monthly$MM_Temp, na.rm = TRUE),
  length.out = 100
)

pred_grid <- expand.grid(
  MM_Temp = temp_seq,
  Source  = factor(c("IL", "OL")),
  # Exclude random effects by assigning dummy levels (exclude terms in predict.gam)
  Tank      = ROvsTemp_monthly$Tank[1],
  colony_id = ROvsTemp_monthly$colony_id[1]
)

# Predict on the link (log) scale excluding random effects, then transform to response scale
pred_raw <- predict(
  ROvsTemp_gamm,
  newdata = pred_grid,
  type = "link",
  se.fit = TRUE,
  exclude = c("s(Tank)", "s(colony_id)")
)

# Convert log predictions and 95% CIs to the original scale (number of planulae)
pred_grid$fit <- ROvsTemp_gamm$family$linkinv(pred_raw$fit)
pred_grid$lwr <- ROvsTemp_gamm$family$linkinv(pred_raw$fit - 1.96 * pred_raw$se.fit)
pred_grid$upr <- ROvsTemp_gamm$family$linkinv(pred_raw$fit + 1.96 * pred_raw$se.fit)


## 2. PLOT FIGURE 5a USING MODEL PREDICTIONS

fig5a <- ggplot() +
  # Raw observations
  geom_point(
    data = ROvsTemp_monthly, 
    aes(x = MM_Temp, y = Monthly_Sum_RO, color = Source),
    size = 3.5, alpha = 0.7
  ) + 
  # 95% Confidence Interval Ribbons from GAMM model
  geom_ribbon(
    data = pred_grid, 
    aes(x = MM_Temp, ymin = lwr, ymax = upr, fill = Source),
    alpha = 0.2
  ) +
  # Fitted GAMM Curves
  geom_line(
    data = pred_grid, 
    aes(x = MM_Temp, y = fit, color = Source),
    linewidth = 1.2
  ) +
  facet_wrap(~Source) + 
  scale_color_manual(values = c("IL" = "#0072B2", "OL" = "#D55E00")) +
  scale_fill_manual(values = c("IL" = "#0072B2", "OL" = "#D55E00")) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Mean monthly tank temperature (°C)",
       y = "Number of planulae (ind.)")

print(fig5a)


  

#### 1.4.2 Analysis####
# Fit GAMM with smooth temperature curves by Source and colony random effect
ROvsTemp_gamm <- gam(
  Monthly_Sum_RO ~ Source + s(MM_Temp, by = Source, k = 6) + 
    s(Tank, bs = "re") + s(colony_id, bs = "re"),
  data = ROvsTemp_monthly,
  family = nb()
)

#Check assumption
gam.check(ROvsTemp_gamm)

summary(ROvsTemp_gamm)

#calculating random effect contribution (%)
gam.vcomp(ROvsTemp_gamm)


####<2. RO>####
#Monthly RO per colony
RO <- read_csv("C:/Users/lamkw/Downloads/3. R test/Nano Project (Louis Laptop)/2. RO.csv")

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
  scale_fill_manual(values = c("#0072B2","#D55E00"))+
  labs(y = "Number of planulae (ind.)", x="Reproductive cycle", fill = "Coral source")+
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

summary(m1)
Anova(m1, type = "III") 

#Posthoc
#[1] Between populations within each cycle
RO_posthoc1 <- emmeans(m1, pairwise~Source|Measured_month)
#[2] Within each population across cycles
RO_posthoc2 <- emmeans(m1, pairwise~Measured_month|Source)


####<3. TLE>####
Size  <- read_csv("C:/Users/lamkw/Downloads/3. R test/Nano Project (Louis Laptop)/3. Size.csv")

Size$Lunar <- as.Date(Size$Lunar, format ='%m/%d/%Y')
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


####*3.1.1 Plot####
Size_lineplot <-
  ggplot(Summary.size, aes(x=Lunar, y=mean, color=Site)) + 
  geom_line(size=1)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5, size=0.7)+
  scale_x_date(limits = c(as.Date("2023-09-01"), as.Date("2024-02-15"))) + 
  ylim(0,300)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  labs(x = "Lunar month", y="Total linear extension (TLE, mm)", color = "Coral source")+
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right")+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
Size_lineplot


####*3.1.2 Analysis####
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



####<3.2 RO vs TLE>####
####* 3.2.1 Plot####
## 1. GENERATE PREDICTIONS FROM THE SIZE/TLE GAMM MODEL ##

# Create a smooth sequence of TLE values (100 points) across the observed range
tle_seq <- seq(
  min(Size$TLE, na.rm = TRUE),
  max(Size$TLE, na.rm = TRUE),
  length.out = 100
)

# Build a prediction grid for both IL and OL sites
pred_grid_size <- expand.grid(
  TLE   = tle_seq,
  Site  = factor(c("IL", "OL")),
  # Supply dummy levels for random effects (will be excluded during prediction)
  Tank = Size$Tank[1],
  ID   = Size$ID[1]
)

# Predict on the log-link scale excluding random intercepts for Tank and ID
pred_raw_size <- predict(
  ROvsGrowth_gamm,
  newdata = pred_grid_size,
  type = "link",
  se.fit = TRUE,
  exclude = c("s(Tank)", "s(ID)")
)

# Back-transform log-scale predictions and 95% CIs to the count scale (RO)
pred_grid_size$fit <- ROvsGrowth_gamm$family$linkinv(pred_raw_size$fit)
pred_grid_size$lwr <- ROvsGrowth_gamm$family$linkinv(pred_raw_size$fit - 1.96 * pred_raw_size$se.fit)
pred_grid_size$upr <- ROvsGrowth_gamm$family$linkinv(pred_raw_size$fit + 1.96 * pred_raw_size$se.fit)


## 2. PLOT FIGURE 5b USING MODEL PREDICTIONS

fig5b <- ggplot() +
  # Raw observations
  geom_point(
    data = Size, 
    aes(x = TLE, y = RO, color = Site),
    size = 3, alpha = 0.7
  ) +
  # 95% Confidence Interval Ribbons from GAMM model
  geom_ribbon(
    data = pred_grid_size, 
    aes(x = TLE, ymin = lwr, ymax = upr, fill = Site),
    alpha = 0.15
  ) +
  # Fitted GAMM Population-Level Curves
  geom_line(
    data = pred_grid_size, 
    aes(x = TLE, y = fit, color = Site),
    linewidth = 1.2
  ) +
  facet_wrap(~Site) +
  scale_color_manual(values = c("IL" = "#0072B2", "OL" = "#D55E00")) +
  scale_fill_manual(values = c("IL" = "#0072B2", "OL" = "#D55E00")) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Total linear extension (TLE, mm)",
    y = "Number of planulae (ind.)"
  )

print(fig5b)



####* 3.2.2 Analysis####
ROvsGrowth_gamm <- gam(
  RO ~ Site + s(TLE, by = Site, k = 6) + 
    s(Tank, bs = "re") + s(ID, bs = "re"),
  data = Size,
  family = nb())

#Check assumption
gam.check(ROvsGrowth_gamm)

# Check the results to see if the curves are linear or non-linear
summary(ROvsGrowth_gamm)


####<4. FvFm>####
#adult
FVFM <- read_csv("C:/Users/lamkw/Downloads/3. R test/Nano Project (Louis Laptop)/4. FVFM.csv")

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
head(Summary)

#Plot
fvfm_lineplot <-
  ggplot(Summary, aes(x=Lunar, y=mean, color=Site)) + 
  geom_line(size=0.8)+
  geom_point(size=2.2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  scale_x_date(limits = c(as.Date("2023-09-01"), as.Date("2024-02-15"))) + 
  ylim(0.5,1)+
  labs(y = "Fv/Fm", x="Lunar month", color = "Coral source")+
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
  labs(y="Temperature (?XC)", x="")

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
setwd("C:\\Users\\lamkw\\Downloads\\3. R test\\Nano Project (Louis Laptop)\\Manuscript Figures")
png("Fig 6a.png", units="mm", width=250, height=150, res=1200)
ggplot(Summary.size, aes(x=Lunar, y=mean, color=Site)) + 
  geom_line(size=1)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5, size=0.7)+
  scale_x_date(limits = c(as.Date("2023-09-01"), as.Date("2024-02-15"))) + 
  ylim(0,300)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  labs(x = "Lunar month", y="Total linear extension (TLE, mm)", color = "Coral source")+
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right")+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
dev.off()

