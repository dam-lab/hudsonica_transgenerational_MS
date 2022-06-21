library(survminer)
library(survival)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(popbio)
library(car)
library(emmeans)
library(sjPlot)
library(lme4)




#################################################################################################################################################################################
#################################################################################################################################################################################
####################################################################### Survival Data analysis ##################################################################################
setwd("C:/Users/james/Documents/Grad_school/OA_hudsonica/Cost_experiments/")


# Read in data table
SurvData <- fread("C:/Users/james/Documents/Grad_school/OA_hudsonica/Cost_experiments/SurvDataTransplant.txt")



# Create uniqe sorting column to group each technical and biological replicate

SurvData <- unite(SurvData,#the data frame
                  unique, #the name of the new column
                  c(Treatment, Rep, Beak), #the existing columns you want to combine
                  remove = FALSE) 

# change number of individuals to numeric format
SurvData$nx <- as.numeric(SurvData$nx)

SurvData$Generation <- 11

SurvData$lx <- as.numeric(SurvData$lx)

SurvData1 <- filter(SurvData, SurvData$lx > 0)

#SurvData1 <- filter(SurvData1, SurvData1$Food == 600)

SurvData.agg <- aggregate(SurvData1$lx,
                          list(SurvData1$Treatment, SurvData1$Rep, SurvData1$Beak, SurvData1$Generation),
                          min)



colnames(SurvData.agg) <- c("Treatment", "Rep", "Beak", "Generation", "lx")

# Modify dataframe to repeat number of rows by the number of remaining individuals
# Other vectors now weighted to reflect the number of remaining individuals for that specific time point
SurvData.kaplan <- SurvData[rep(seq(nrow(SurvData)), 
                                       SurvData$nx),#the column which holds the value you wish to have everything repeated by
                                   ]
# Designate "Survival" for the individuals that survived to adulthood
# help from: https://stackoverflow.com/questions/51642416/changing-number-of-rows-based-on-changing-vector-in-r

SurvData.kaplan <- SurvData.kaplan %>%
  group_by(Treatment, Rep, Beak) %>%
  mutate(
    survival = case_when(
      Food == 0 ~ 1,# no starved treatments survived
      nx == min(nx) ~ 0, #need to use the "==" and VALUE OF 1 EQUALS DEATH HAPPENING!!!!!!
      nx > min(nx) ~ 1 )) %>% 
  group_by(unique) %>%
  mutate(keep = max(nx) - nx) %>% 
  group_by(unique, keep) %>% #need to group by unique and keep to get all the animals that died in between the start and the lowest amount of dead animals
  filter(survival == 1 & nx == min(nx) |
           survival == 0) %>% 
  filter(row_number() %in% 1:unique(keep) |
           survival == 0) %>% 
  select(-keep) %>% 
  ungroup()


# Remove starting point since they all started alive
SurvData.kaplan <- SurvData.kaplan[!SurvData.kaplan$time == 0,] 

# Create a category for unique treatments and food
SurvData.kaplan <- unite(SurvData.kaplan,#the data frame
                         Treat.food, #the name of the new column
                         c(Treatment, Food), #the existing columns you want to combine
                         remove = FALSE) 

## Plot only food limited and food replete treatments

# Remove starved individuals to plot 
SurvData.kaplan.no.starve <- SurvData.kaplan[!SurvData.kaplan$Food == 0,]

# New grouping vectors for line and developmental treatment
SurvData.kaplan.no.starve$line <- if_else(SurvData.kaplan.no.starve$Treatment == 1 | SurvData.kaplan.no.starve$Treatment == 6, 
                                          "AA", "HH")

SurvData.kaplan.no.starve$Treatment2 <- if_else(SurvData.kaplan.no.starve$Treatment == 1 | SurvData.kaplan.no.starve$Treatment == 5, 
                                                "AA", "HH")

SurvData.kaplan.no.starve$Food <- as.factor(SurvData.kaplan.no.starve$Food)

# Create a survival object
surv_object <- Surv(time = SurvData.kaplan.no.starve$time, event = SurvData.kaplan.no.starve$survival)

# Create a model to fit the object
surv_object_fit <- survfit(surv_object~Treatment+Food, data = SurvData.kaplan.no.starve)

surv_object

surv_object_fit


# log-rank test results for resulting model
surv_pvalue(surv_object_fit, SurvData.kaplan.no.starve)

# Two-way ANOVA comparing day-specific survivorship by Food abundance and Line
surv.anova <- aov(lx~Food*line, data = SurvData.kaplan.no.starve)
summary(surv.anova)



# Plot the fits
surv.plot <- ggsurvplot_facet(surv_object_fit,
                        data = SurvData.kaplan.no.starve,
                        facet.by = "Food",
                        #colors are the "treatments"
                        palette = c("cornflowerblue", #AAAA
                                    "brown2", #HHHH
                                    "cornflowerblue", #AAHH
                                    "brown2" #HHAA
                        ),
                        #linetypes are the "lines"
                        linetype =c("solid", #AAAA 250
                                    "solid", #AAAA 800
                                    "dashed", #HHHH 250
                                    "dashed", #HHHH 800
                                    "dashed", #AAHH 250
                                    "dashed", #AAHH 800
                                    "solid", #HHAA 250
                                    "solid" #HHAA 800
                        ),
                        conf.int = TRUE,
                        break.time.by = 10
                        
)


surv.plot
surv.plot.large <- ggpar(surv.plot, font.x = font.size, font.y = font.size, font.tickslab = font.size)

surv.plot.large

aspect.ratio <- 1.15
ggsave(filename = "survplot_square.pdf", surv.plot.large, height = 8.5, width = 8.5, units = "in", device = "pdf")



#################################################################################################################################################################################
#################################################################################################################################################################################
####################################################################### Egg Production and Hatching #############################################################################


cost.dir <- "C:/Users/james/Documents/Grad_school/OA_hudsonica/Cost_experiments/"
Cost.epr <- fread(paste(cost.dir, "Feeding_cost_epr.txt", sep = ""))

#Cost.epr1 <- filter(Cost.epr, Cost.epr$Food == 600)
#Cost.epr1 <- filter(Cost.epr1, Cost.epr1$Treatment < 5)


# Pairwise comparison of EPR results
Cost.epr$Food.f <- as.factor(as.numeric(Cost.epr$Food))
Cost.epr$Treatment.f <- as.factor(as.numeric(Cost.epr$Treatment))
pairwise.t.test(Cost.epr$EPR, Cost.epr$Food.f:Cost.epr$Treatment.f, p.adj = "none")


l <- lm(EPR~as.factor(Treatment)*as.factor(Food), Cost.epr)
Anova(l)



summary(l)
tab_model(l)


epr.emm <- emmeans(l, pairwise ~ Treatment | Food)
epr.tukey.groups <- CLD(epr.emm)
epr.pw.groups <- CLD(epr.emm, p.adjust.methods = "none")

fwrite(epr.pw.groups, file = paste(cost.dir,"Cost_epr_pw_groups.txt", sep = ""), sep = "\t")
fwrite(epr.tukey.groups, file = paste(cost.dir,"Cost_epr_tukey_groups.txt", sep = ""), sep = "\t")


Cost.epr.no.starve <- filter(Cost.epr, Food > 1)


Cost.epr.no.starve <- unite(Cost.epr.no.starve,
                            Treat.Bin,
                            c(Treatment, Bin),
                            remove = F)

l2 <- lmer(EPR~as.factor(Treatment)*as.factor(Food) + (1|Treat.Bin), Cost.epr.no.starve)
Anova(l2)



summary(l2)
tab_model(l2)


epr.emm.no.starve <- emmeans(l2, pairwise ~ Treatment | Food)
pairs(epr.emm.no.starve)


epr.emm.no.starve.T <- pairs(epr.emm.no.starve)
epr.emm.no.starve.T <- epr.emm.no.starve.T$emmeans
epr.emm.no.starve.T <- tidy(epr.emm.no.starve.T)

epr.emm.no.starve <- pairs(epr.emm.no.starve, adjust = "none")
epr.emm.no.starve <- epr.emm.no.starve$emmeans
epr.emm.no.starve <- tidy(epr.emm.no.starve)


epr.emm.no.starve2 <- emmeans(l2, pairwise ~ factor(Food) | factor(Treatment))

epr.emm.no.starve2.Tukey <- pairs(epr.emm.no.starve2)
epr.emm.no.starve2.Tukey <- epr.emm.no.starve2.Tukey$emmeans
epr.emm.no.starve2.Tukey <- tidy(epr.emm.no.starve2.Tukey)

epr.emm.no.starve2.pw <- pairs(epr.emm.no.starve2, adjust = "none")
epr.emm.no.starve2.pw <- epr.emm.no.starve2.pw$emmeans
epr.emm.no.starve2.pw <- tidy(epr.emm.no.starve2.pw)



fwrite(epr.emm.no.starve, paste0(cost.dir, "EPR_cost_pw.txt"), sep = "\t")
fwrite(epr.emm.no.starve.T, paste0(cost.dir, "EPR_cost_Tukey.txt"), sep = "\t")
fwrite(epr.emm.no.starve2.pw, paste0(cost.dir, "EPR_cost_pw2.txt"), sep = "\t")
fwrite(epr.emm.no.starve2.Tukey, paste0(cost.dir, "EPR_cost_Tukey2.txt"), sep = "\t")


epr.tukey.groups.no.starve <- CLD(epr.emm.no.starve)
epr.pw.groups.no.starve <- CLD(epr.emm.no.starve, p.adjust.methods = "none")

fwrite(epr.pw.groups.no.starve, file = paste(cost.dir,"Cost_epr_pw_groups_no_starve.txt", sep = ""), sep = "\t")
fwrite(epr.tukey.groups.no.starve, file = paste(cost.dir,"Cost_epr_tukey_groups_no_starve.txt", sep = ""), sep = "\t")


## Line x environment interactions

Cost.epr.no.starve$Line <- case_when(Cost.epr.no.starve$Treatment == 1 ~ "AM",
                                     Cost.epr.no.starve$Treatment == 4 ~ "GH",
                                     Cost.epr.no.starve$Treatment == 5 ~ "GH",
                                     Cost.epr.no.starve$Treatment == 6 ~ "AM")

Cost.epr.no.starve$Environment <- case_when(Cost.epr.no.starve$Treatment == 1 ~ "AM",
                                           Cost.epr.no.starve$Treatment == 4 ~ "GH",
                                           Cost.epr.no.starve$Treatment == 5 ~ "AM",
                                           Cost.epr.no.starve$Treatment == 6 ~ "GH")


l3 <- lmer(EPR ~ factor(Line)*factor(Environment) + (1|Treat.Bin), data = Cost.epr.no.starve)

summary(l3)
Anova(l3)
tab_model(l3)


# 600 food only
l4 <- lmer(EPR ~ factor(Line)*factor(Environment) + (1|Treat.Bin), data = subset(Cost.epr.no.starve, Food > 250))

summary(l4)
Anova(l4)
tab_model(l4)


# 250 food only
l5 <- lmer(EPR ~ factor(Line)*factor(Environment) + (1|Treat.Bin), data = subset(Cost.epr.no.starve, Food <600))

summary(l5)
Anova(l5)
tab_model(l5)
## Plot EPR results


Cost.epr.mean <- Cost.epr %>% 
  group_by(Treatment, Food) %>% 
  summarise(mean = mean(EPR, na.rm = TRUE),
            sd = sd(EPR, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

fwrite(Cost.epr.mean, file = paste0(cost.dir,"Cost_epr_mean.txt"), sep = "\t")

Cost.epr.mean$Line <- case_when(Cost.epr.mean$Treatment == 1 ~ "AM",
                                Cost.epr.mean$Treatment == 4 ~ "GH",
                                Cost.epr.mean$Treatment == 5 ~ "GH",
                                Cost.epr.mean$Treatment == 6 ~ "AM")

Cost.epr.mean$Treatment2 <- case_when(Cost.epr.mean$Treatment == 1 ~ "AM",
                                      Cost.epr.mean$Treatment == 4 ~ "GH",
                                      Cost.epr.mean$Treatment == 5 ~ "AM",
                                      Cost.epr.mean$Treatment == 6 ~ "GH")


Cost.epr.mean$Food <-as.numeric(Cost.epr.mean$Food)

# Plot only food limited and food replete conditions
Cost.epr.mean.no.starve <- filter(Cost.epr.mean, Food > 1)


Cost.epr.mean.no.starve$Food <- as.factor(Cost.epr.mean.no.starve$Food)


# Create variables for grouping when plotting
Cost.epr.mean.no.starve$Order <- case_when(Cost.epr.mean.no.starve$Treatment == 1 ~ "A",
                                           Cost.epr.mean.no.starve$Treatment == 4 ~ "C",
                                           Cost.epr.mean.no.starve$Treatment == 5 ~ "D",
                                           Cost.epr.mean.no.starve$Treatment == 6 ~ "B")



# Create variables for grouping when plotting
Cost.epr.mean.no.starve$Order <- case_when(Cost.epr.mean.no.starve$Treatment == 1 ~ "A",
                                           Cost.epr.mean.no.starve$Treatment == 4 ~ "C",
                                           Cost.epr.mean.no.starve$Treatment == 5 ~ "D",
                                           Cost.epr.mean.no.starve$Treatment == 6 ~ "B")



Cost.epr.mean.no.starve <- unite(Cost.epr.mean.no.starve,
                                 Treat.food,
                                 c(Treatment2, Food),
                                 remove = F)

Cost.epr.mean.no.starve <- unite(Cost.epr.mean.no.starve,
                                 Line.food,
                                 c(Line, Food),
                                 remove = F)


Cost.epr.mean.no.starve$Treat.food <- factor(Cost.epr.mean.no.starve$Treat.food, levels = c("AM_600", "GH_600", "AM_250", "GH_250"))

font.size <- 32

pd=0.5
EPR.point <- ggplot(data = Cost.epr.mean.no.starve) +
  
  
  geom_errorbar(aes(x=Treat.food, 
                    ymin=lower.ci, 
                    ymax=upper.ci,
                    group = Line.food), 
                colour = "black", 
                size=2, 
                width = 0.4, 
                position = position_dodge(width = pd)
  )+
  
  geom_line(aes(x=Treat.food, 
                y=mean,
                group = Line.food),
            size = 2,
            position = position_dodge(width = pd)
  )+
  
  
  geom_point(aes(x=Treat.food, 
                 y=mean, 
                 colour = factor(Treatment2),
                 shape = factor(Line),
                 group = Line.food),
             size = 13.5,
             position = position_dodge(width = pd)
  )+
  
  
  geom_vline(xintercept = 2.5, linetype = "dotted") +
  
  
  
  scale_shape_manual(name = "Line",
                     values = c(16, 17))+
  
  scale_color_manual(name = "Environment",
                     values = c("cornflowerblue", "brown2"))+
  
  scale_x_discrete(labels = c("AM_600" = "AM",
                              "GH_600" = "GH",
                              "AM_250" = "AM",
                              "GH_250" = "GH"))+
  
  theme_classic()+
  theme(legend.text.align = 0,
        legend.text = element_text(size = font.size),
        legend.title = element_text(size = font.size),
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size),
        axis.text.x = element_text(size = font.size, colour = "black"),
        axis.text.y = element_text(size = font.size, colour = "black"),
        legend.position = "bottom")+
  xlab(expression ("Environment"))+
  ylab(expression (atop("Egg Production Rate",(day^-1~female^-1))))



EPR.point

ggsave(filename = "EPR_point_cost.pdf", plot = EPR.point, path = cost.dir, height = 9.75, width = 8.5, units = "in")


### Hatching success 
Cost.hf.mean <- Cost.epr %>% 
  group_by(Treatment, Food) %>% 
  summarise(mean = mean(HF, na.rm = TRUE),
            sd = sd(HF, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

fwrite(Cost.hf.mean, file = paste0(cost.dir,"Cost_hf_mean.txt"), sep = "\t")

Cost.hf.mean$Line <- case_when(Cost.hf.mean$Treatment == 1 ~ "AM",
                                Cost.hf.mean$Treatment == 4 ~ "GH",
                                Cost.hf.mean$Treatment == 5 ~ "GH",
                                Cost.hf.mean$Treatment == 6 ~ "AM")

Cost.hf.mean$Treatment2 <- case_when(Cost.hf.mean$Treatment == 1 ~ "AM",
                                      Cost.hf.mean$Treatment == 4 ~ "GH",
                                      Cost.hf.mean$Treatment == 5 ~ "AM",
                                      Cost.hf.mean$Treatment == 6 ~ "GH")


Cost.hf.mean$Food <-as.numeric(Cost.hf.mean$Food)

# Plot only food limited and food replete conditions
Cost.hf.mean.no.starve <- filter(Cost.hf.mean, Food > 1)


Cost.hf.mean.no.starve$Food <- as.factor(Cost.hf.mean.no.starve$Food)


# Create variables for grouping when plotting
Cost.hf.mean.no.starve$Order <- case_when(Cost.hf.mean.no.starve$Treatment == 1 ~ "A",
                                           Cost.hf.mean.no.starve$Treatment == 4 ~ "C",
                                           Cost.hf.mean.no.starve$Treatment == 5 ~ "D",
                                           Cost.hf.mean.no.starve$Treatment == 6 ~ "B")



# Create variables for grouping when plotting
Cost.hf.mean.no.starve$Order <- case_when(Cost.hf.mean.no.starve$Treatment == 1 ~ "A",
                                           Cost.hf.mean.no.starve$Treatment == 4 ~ "C",
                                           Cost.hf.mean.no.starve$Treatment == 5 ~ "D",
                                           Cost.hf.mean.no.starve$Treatment == 6 ~ "B")



Cost.hf.mean.no.starve <- unite(Cost.hf.mean.no.starve,
                                 Treat.food,
                                 c(Treatment2, Food),
                                 remove = F)

Cost.hf.mean.no.starve <- unite(Cost.hf.mean.no.starve,
                                 Line.food,
                                 c(Line, Food),
                                 remove = F)


Cost.hf.mean.no.starve$Treat.food <- factor(Cost.hf.mean.no.starve$Treat.food, levels = c("AM_600", "GH_600", "AM_250", "GH_250"))

HF.point <- ggplot(data = Cost.hf.mean.no.starve) +
  
  
  geom_errorbar(aes(x=Treat.food, 
                    ymin=lower.ci, 
                    ymax=upper.ci,
                    group = Line.food), 
                colour = "black", 
                size=2, 
                width = 0.4, 
                position = position_dodge(width = pd)
  )+
  
  geom_line(aes(x=Treat.food, 
                y=mean,
                group = Line.food),
            size = 2,
            position = position_dodge(width = pd)
  )+
  
  
  geom_point(aes(x=Treat.food, 
                 y=mean, 
                 colour = factor(Treatment2),
                 shape = factor(Line),
                 group = Line.food),
             size = 13.5,
             position = position_dodge(width = pd)
  )+
  
  
  geom_vline(xintercept = 2.5, linetype = "dotted") +
  
  
  
  scale_shape_manual(name = "Line",
                     values = c(16, 17))+
  
  scale_color_manual(name = "Environment",
                     values = c("cornflowerblue", "brown2"))+
  
  scale_x_discrete(labels = c("AM_600" = "AM",
                              "GH_600" = "GH",
                              "AM_250" = "AM",
                              "GH_250" = "GH"))+
  
  theme_classic()+
  theme(legend.text.align = 0,
        legend.text = element_text(size = font.size),
        legend.title = element_text(size = font.size),
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size),
        axis.text.x = element_text(size = font.size, colour = "black"),
        axis.text.y = element_text(size = font.size, colour = "black"),
        legend.position = "bottom")+
  xlab(expression ("Environment"))+
  ylab("Hatching Success")



HF.point

ggsave(filename = "HF_point_cost.pdf", plot = HF.point, path = cost.dir, height = 9.75, width = 8.5, units = "in")


#################################################################################################################################################################################
#################################################################################################################################################################################
####################################################################### Calculating Fitness ########################################################################



SurvData$lx <- as.numeric(SurvData$lx)

## Extract the sex ratio data for scaling fecundity
SurvData$M.Ratio <- as.numeric(SurvData$M.Ratio)
SurvData$F.Ratio <- as.numeric(SurvData$F.Ratio)

SurvDataSex <- SurvData %>%
  filter(F.Ratio >= 0 | M.Ratio >= 0) %>%
  group_by(Treatment, Rep, Beak) %>%
  summarise(F.Ratio = last(F.Ratio))

SurvDataSex.Mean <- SurvDataSex %>%
  group_by(Treatment) %>%
  summarise(F.Ratio = mean(F.Ratio, na.rm = TRUE)) %>%
  as.data.frame(SurvDataSex.Mean)

SurvDataSex.Mean




# Create a vector formatted for day-specific survivorship as it fits along the off-diagonal of a Leslie Matrix
# See Caswell, H. 2001. Matrix Population Models for further details


SurvData1 <- SurvData %>%
  group_by(Treatment) %>%
  mutate(lx_diag = if_else(time == 0, 1, as.numeric(lx/lag(lx, default = first(lx))))) %>% ## always start with 100%
  mutate(lx_diag = if_else(lx_diag <= 1.000000, lx_diag, lag(lx_diag, n=1, default =last(lx_diag)))) %>%
  mutate(days = if_else(time == 0, 1, as.numeric(time-lag(time)))) # create a new column for the number of days spent at the respective survivorships

# Survivorship is reflective of prior day. Therefore, there can be no 0 survivorship in this vector.
# If all animals die, then the vector ends and the matrix is truncated at the appropriate time
SurvData1 <- filter(SurvData1, lx_diag > 0) 


# Check if there is any super survivorship. There can be no survivorship > 1.
# No individuals can be lost and reappear


any(SurvData1$lx_diag > 1)



SurvData1 <- as.data.frame(SurvData1)


# elongate the data frame to make it reflect actual days spent over the experiment. This essentially changes the matrix to a daily matrix.

SurvData1 <- SurvData1[rep(seq(nrow(SurvData1)), SurvData1$days),] 

# remove starved samples
SurvData1 <- SurvData1[!SurvData1$Food == 0,]





#################################################################################################################################################################################
#################################################################################################################################################################################
####################################################################### Development Time ########################################################################################


Dev.time <- filter(SurvData, SurvData$Cdev > 0)
Dev.time <- Dev.time[rep(seq(nrow(Dev.time)), Dev.time$Cdev),]

Dev.time.600 <- Dev.time %>% 
  filter(Food == 600) %>% 
  filter(Treatment <5) %>% 
  select(,-6) %>% 
  mutate(Generation = 11)


Ndev.time <- filter(SurvData, SurvData$Ndev > 0)
Ndev.time <- Ndev.time[rep(seq(nrow(Ndev.time)), Ndev.time$Ndev),]

Ndev.time.600 <- Ndev.time %>% 
  filter(Food == 600) %>% 
  filter(Treatment <5) %>% 
  select(,-6) %>% 
  mutate(Generation = 11)

Dev.time.sum <- Dev.time %>%
  group_by(Treatment, Rep, Food) %>%
  summarise(mean = mean(time, na.rm = TRUE),
            sd = sd(time, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)


Dev.time.sum <- Dev.time.sum[!Dev.time.sum$Food == 0,]

Dev.time.sum <- unite(Dev.time.sum,#the data frame
                      unique, #the name of the new column
                      c(Treatment, Rep, Food), #the existing columns you want to combine
                      remove = FALSE)
Dev.time.sum.short <- Dev.time.sum[,-c(6:10)]


## for plotting with F0 - F4
Devtimesum <- Dev.time %>%
  group_by(Treatment, Food) %>%
  summarise(mean = mean(time, na.rm = TRUE),
            sd = sd(time, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

fwrite(Devtimesum, file = paste0(cost.dir,"Dev_time_sum_cost.txt"), sep = "\t")

## Add the dev.time to the survival table
SurvData1 <- inner_join(SurvData1,
                       Dev.time.sum.short,
                       by = c("Treatment", "Rep", "Food"))

SurvData1 <- SurvData1 %>% rename(dev.time = mean) # the new name of the column has to come first

SurvData1 <- unite(SurvData1,#the data frame
                   unique2, #the name of the new column
                   c(Treatment, Food), #the existing columns you want to combine
                   remove = FALSE)

Survival.list <- split(SurvData1, f = SurvData1$unique2) # don't create the list with a list of organizers. Then you create a complete matrix of samples which include some tables with no data



# lists within lists
Survival.list <- lapply(Survival.list, function (x) split(x, f = x$Rep))

View(Survival.list)

#################################################################################################################################################################################
#################################################################################################################################################################################
######################################################################### EPR data #########################################################################


Cost.epr$fecundity <- Cost.epr$EPR*Cost.epr$HF

# These sex ratio values are based on previously calculated values in survival data

Cost.epr$F.ratio <- case_when(Cost.epr$Treatment == 4 ~ 0.40,
                            Cost.epr$Treatment == 5 ~ 0.43,
                            Cost.epr$Treatment == 1 ~ 0.31,
                            Cost.epr$Treatment == 6 ~ 0.41)


# Caluclate per capita sex-specific fecundity

Cost.epr$sex.spec.fecundity <- Cost.epr$fecundity*Cost.epr$F.ratio

# remove starved animals
Cost.epr <- Cost.epr[!Cost.epr$Food == 0,]

Cost.epr <- unite(Cost.epr,#the data frame
                      unique, #the name of the new column
                      c(Treatment, Food), #the existing columns you want to combine
                      remove = FALSE)


### NOTE: lists of survivorship and EPR data MUST have same names and same length of items

EPR.list <- split(Cost.epr, f = Cost.epr$unique)

View(EPR.list)


names(Survival.list)
names(EPR.list)


names.surv <- names(Survival.list)
names.surv <- as.list(names.surv)
names.surv

# Create a dummy table to add data to

lambda.results.cost <- data.frame(Variables = "dummy", Rep = 0, surv = 0, epr = 0, hf = 0, dev.time = 0, lambda = 0)



for (i in 1:length(Survival.list)) { #list of names for each gen and treatment
  
  diag.ls <- Survival.list[[i]]
  
  epr.df <- EPR.list[[i]] # select the appropriate data frame to pull epr values from
  
  for (j in 1: length(diag.ls)) {
    
    diag.df <-  diag.ls[[j]] # temporary list used to index at each point
    
    df.rep <- mean(diag.df$Rep)
    
    diag.v <- diag.df$lx_diag # extract the vector of interest for the diagonal
    
    leslie.matrix <- diag(diag.v) # create the matrix without the top row added
    
    dev.time.value <- as.integer(mean(diag.df$dev.time)) # calculate the development time
    
    zero <- matrix(0, nrow = 1, ncol = (dev.time.value)) # create an empty, one-row matrix to add the epr values to. This will be added to the leslie matrix as the top row.
    
    epr.vector <- epr.df$sex.spec.fecundity # select the vector with sex spec fecundity
    
    epr.vector2 <- epr.df$EPR
    
    hf.vector <- epr.df$HF
    
    epr.count <- as.integer(dim(leslie.matrix)[1]-dev.time.value) # calculate the number of days to incorporate egg production based on development time
    
    if (epr.count < 1) {
      
      epr.count <- 1
      
    } 
    
    
    for(k in 1:length(epr.vector)) {
      
      epr.value <- epr.vector[k] # use the sex specific fecundity values in sequence with the same survival matrix
      
      epr.value2 <- epr.vector2[k]
      
      hf.value <- hf.vector[k]
      
      if (is.na(epr.value) == TRUE) {
        
        epr.value <- 0 ## use an arrow when assigning numbers, not a "=="
        
      } 
      
      
      
      fecundity.row <- t(c(zero, rep(epr.value, epr.count))) # combine the zero row and the epr.value for the matrix and transpose it to make it a row
      
      if (ncol(fecundity.row) > ncol(leslie.matrix)) { # if the dev.time is somehow greater than the survival matrix, then only make the last column reflective of epr
        
        delete <- ncol(fecundity.row)-ncol(leslie.matrix) # find the number of days that the dev time is greater than the survivorship
        
        fecundity.row <- fecundity.row[,-c(1:delete)] # delete those days from the fecundity row to make it the same size as the leslie matrix
        
      }
      
      leslie.matrix1 <- rbind(fecundity.row, leslie.matrix) # add the fecundity row to the matrix
      
      matrix.final <- leslie.matrix1[-nrow(leslie.matrix1),] # eliminate the last row to make it a square
      
      eigen.calcs <- eigen.analysis(matrix.final, zero = FALSE) # calculate the eigen values
      
      lambda.value <- eigen.calcs$lambda1 # extract the dominant eigen values
      
      lambda.row <- data.frame(Variables = names.surv[i], 
                               Rep = df.rep, 
                               surv = min(diag.df$lx), 
                               epr = epr.value2,
                               hf = hf.value,
                               dev.time = dev.time.value,
                               lambda = lambda.value) # create a 1x2 data frame to add to the end of the final data frame
      # data frame has to have the same colnames in order to rbind
      
      colnames(lambda.row) <- colnames(lambda.results.cost) # make sure the data frames have the same names
      
      lambda.results.cost <- bind_rows(lambda.results.cost, lambda.row) # append the data frame with new results
      
      
    }
    
    
    
  }
}

######################################################################################################################################################

lambda.results.cost <- fread("lambda_results_cost_devtime_with_rep_hudsonica.txt")


View(lambda.results.cost)
lambda.results.cost <- separate(lambda.results.cost, "Variables", into = c("Treatment","Food"))


# remove the first row as a last step
lambda.results.cost <- lambda.results.cost[-1,] 


View(lambda.results.cost)

fwrite(lambda.results.cost, file = "lambda_results_cost_f11_surv_epr_hf_devtime.txt", sep = "\t")
# calculate malthusian parameter
#lambda.results.cost$Malthusian <- log10(lambda.results.cost$lambda)


# create continuous generation vector for anova and plotting
lambda.results.cost$Food.c <- as.numeric(as.character(lambda.results.cost$Food)) 

lambda.results.cost$Food <- as.factor(as.numeric(lambda.results.cost$Food))
lambda.results.cost$Treatment <- as.factor(lambda.results.cost$Treatment)


# create the model to use for pairwise comparisons -- generation must be a factor, not a continuous variable
lambda.model.2 <- lm(lambda~Food*Treatment, data=lambda.results.cost)


emmip(lambda.model.2, Treatment ~ Food)

cost.emm_1 <- emmeans(lambda.model.2, pairwise ~ Treatment | Food)

# Tukey comparison
p <- pairs(cost.emm_1)
p


lambda.results.cost <- unite(lambda.results.cost,#the data frame
                                  Treat.food, #the name of the new column
                                  c(Treatment, Food), #the existing columns you want to combine
                                  remove = FALSE)


library(broom)
# non-adjusted pairwise t-test
p2 <- tidy(pairwise.t.test(lambda.results.cost$lambda, lambda.results.cost$Treat.food, p.adjust.method = "none"))


fwrite(lambda.results.cost, file = "lambda_results_cost_devtime_with_rep_hudsonica.txt", sep = "\t")
fwrite(p2, file = "Lambda_cost_pairwise_ttest_devtime_hudsonica.txt", sep = "\t")

# model for writing the tukey results
lambda.model.3 <- aov(lambda~Treat.food, data = lambda.results.cost)
p3 <- TukeyHSD(lambda.model.3, "Treat.food")
p3.tukey <- tidy(p3)
fwrite(p3.tukey, file = "lambda_cost_devtime_hudsonica_Tukey.txt", sep = "\t")



library(emmeans)
library(tidyr)


lambda.results.cost <- unite(lambda.results.cost,#the data frame
                             Treat.rep, #the name of the new column
                             c(Treatment, Rep), #the existing columns you want to combine
                             remove = FALSE)

l <- lmer(lambda~factor(Treatment)*factor(Food) + (1|Treat.rep), 
        # REML = FALSE,
        # family = gaussian,
        data = lambda.results.cost)

lambda.emm <- emmeans(l, pairwise ~ factor(Treatment) | factor(Food))

lambda.emm.pairs <- pairs(lambda.emm, adjust = "none")
lambda.emm.pairs <- tidy(lambda.emm.pairs$emmeans)

lambda.emm.Tukey <- pairs(lambda.emm)
lambda.emm.Tukey <- tidy(lambda.emm.Tukey$emmeans)



lambda.emm2 <- emmeans(l, pairwise ~ factor(Food) | factor(Treatment))

lambda.emm.pairs2 <- pairs(lambda.emm2, adjust = "none")
lambda.emm.pairs2 <- tidy(lambda.emm.pairs2$emmeans)

lambda.emm.Tukey2 <- pairs(lambda.emm2)
lambda.emm.Tukey2 <- tidy(lambda.emm.Tukey2$emmeans)


fwrite(lambda.emm.pairs, paste0(cost.dir, "Lambda_contrasts1.txt"), sep = "\t")
fwrite(lambda.emm.pairs2, paste0(cost.dir, "Lambda_contrasts2.txt"), sep = "\t")
fwrite(lambda.emm.Tukey, paste0(cost.dir, "Lambda_Tukey1.txt"), sep = "\t")
fwrite(lambda.emm.Tukey2, paste0(cost.dir, "Lambda_Tukey2.txt"), sep = "\t")





lambda.tukey.groups <- CLD(lambda.emm)

fwrite(lambda.tukey.groups, file = "Lambda_cost_devtime_hudsonica_tukey_groups.txt", sep = "\t")



l2 <- lm(lambda~unique, lambda.results.cost)

library(agricolae)

test <- LSD.test(l2,trt = "unique")
lambda.lsd.groups <- test$groups
lambda.lsd.groups$unique <- rownames(lambda.lsd.groups)

fwrite(lambda.lsd.groups, file = paste0(cost.dir,"lambda_pairwise_devtime_hudsonica_groups.txt"), sep = "\t")



# evaluate line x environment interactions
lambda.results.cost$Line <- case_when(lambda.results.cost$Treatment == 1 ~ "AM",
                                      lambda.results.cost$Treatment == 4 ~ "GH",
                                      lambda.results.cost$Treatment == 5 ~ "GH",
                                      lambda.results.cost$Treatment == 6 ~ "AM")


lambda.results.cost$Environment <- case_when(lambda.results.cost$Treatment == 1 ~ "AM",
                                             lambda.results.cost$Treatment == 4 ~ "GH",
                                             lambda.results.cost$Treatment == 5 ~ "AM",
                                             lambda.results.cost$Treatment == 6 ~ "GH")


library(glmmTMB)
fit_zigauss <- glmmTMB(lambda~factor(Line)*factor(Environment)+(1|Treat.rep),
                         #family = truncated_poisson,
                         data = lambda.results.cost,
                         
                         ziformula = ~.)


summary(fit_zigauss)
tab_model(fit_zigauss)
Anova(fit_zigauss)


fit_zigauss.2 <- glmmTMB(lambda~factor(Line)*factor(Environment)+(1|Treat.rep),
                       #family = truncated_poisson,
                       data = subset(lambda.results.cost, Food == "600"),
                       
                       ziformula = ~.)


summary(fit_zigauss.2)
tab_model(fit_zigauss.2)
Anova(fit_zigauss.2)



fit_zigauss.3 <- glmmTMB(lambda~factor(Line)*factor(Environment)+(1|Treat.rep),
                         #family = truncated_poisson,
                         data = subset(lambda.results.cost, Food == "250"),
                         
                         ziformula = ~.)


summary(fit_zigauss.3)
tab_model(fit_zigauss.3)
Anova(fit_zigauss.3)

###### Lambda Point Plot #####


lambda.results.cost$Order <- case_when(lambda.results.cost$Treatment == 1 ~ "A",
                                  lambda.results.cost$Treatment == 4 ~ "C",
                                  lambda.results.cost$Treatment == 5 ~ "D",
                                  lambda.results.cost$Treatment == 6 ~ "B")


lambda.results.cost$Line <- case_when(lambda.results.cost$Treatment == 1 ~ "AM",
                                   lambda.results.cost$Treatment == 4 ~ "GH",
                                   lambda.results.cost$Treatment == 5 ~ "GH",
                                   lambda.results.cost$Treatment == 6 ~ "AM")

lambda.results.cost$Treatment2 <- case_when(lambda.results.cost$Treatment == 1 ~ "AM",
                                         lambda.results.cost$Treatment == 4 ~ "GH",
                                         lambda.results.cost$Treatment == 5 ~ "AM",
                                         lambda.results.cost$Treatment == 6 ~ "GH")

lambda.results.cost$Food.f<- as.factor(lambda.results.cost$Food)


setwd(cost.dir)
fwrite(lambda.results.cost, "lambda_results_cost_devtime_with_rep.txt", sep = "\t")



lambda.list.mean <- lambda.results.cost %>%
  group_by(Treatment, Food) %>%
  summarise(mean = mean(lambda, na.rm = TRUE),
            sd = sd(lambda, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)

fwrite(lambda.list.mean, file = "lambda_mean_devtime_hudsonica.txt", sep = "\t")

lambda.list.mean$Line <- case_when(lambda.list.mean$Treatment == 1 ~ "AM",
                                           lambda.list.mean$Treatment == 4 ~ "GH",
                                           lambda.list.mean$Treatment == 5 ~ "GH",
                                           lambda.list.mean$Treatment == 6 ~ "AM")

lambda.list.mean$Treatment2 <- case_when(lambda.list.mean$Treatment == 1 ~ "AM",
                                                lambda.list.mean$Treatment == 4 ~ "GH",
                                                lambda.list.mean$Treatment == 5 ~ "AM",
                                                lambda.list.mean$Treatment == 6 ~ "GH")


lambda.list.mean$Order <- case_when(lambda.list.mean$Treatment == 1 ~ "A",
                                            lambda.list.mean$Treatment == 4 ~ "C",
                                            lambda.list.mean$Treatment == 5 ~ "D",
                                            lambda.list.mean$Treatment == 6 ~ "B")

font.size <- 32


lambda.list.mean <- unite(lambda.list.mean,
                          Treat.food,
                          c(Treatment2, Food),
                          remove = F)

lambda.list.mean <- unite(lambda.list.mean,
                          Line.food,
                          c(Line, Food),
                          remove = F)


lambda.list.mean$Treat.food <- factor(lambda.list.mean$Treat.food, levels = c("AM_600", "GH_600", "AM_250", "GH_250"))

pd=0.5
Lambda.point <- ggplot(data = lambda.list.mean) +
  
  
  geom_errorbar(aes(x=Treat.food, 
                    ymin=lower.ci, 
                    ymax=upper.ci,
                    group = Line.food), 
                colour = "black", 
                size=2, 
                width = 0.4, 
                position = position_dodge(width = pd)
  )+
  
  geom_line(aes(x=Treat.food, 
                y=mean,
                group = Line.food),
            size = 2,
            position = position_dodge(width = pd))+
  
  
  geom_point(aes(x=Treat.food, 
                 y=mean, 
                 colour = factor(Treatment2),
                 shape = factor(Line),
                 group = Line.food),
             size = 13.5,
             position = position_dodge(width = pd))+
  
  
  geom_vline(xintercept = 2.5, linetype = "dotted") +
  
  
  
  scale_shape_manual(name = "Line",
                     values = c(16, 17))+
  
  scale_color_manual(name = "Environment",
                     values = c("cornflowerblue", "brown2"))+
  
  scale_x_discrete(labels = c("AM_600" = "AM",
                              "GH_600" = "GH",
                              "AM_250" = "AM",
                              "GH_250" = "GH"))+
  
  theme_classic()+
  theme(legend.text.align = 0,
        legend.text = element_text(size = font.size),
        legend.title = element_text(size = font.size),
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size),
        axis.text.x = element_text(size = font.size, colour = "black"),
        axis.text.y = element_text(size = font.size, colour = "black"),
        legend.position = "bottom")+
  xlab(expression ("Environment"))+
  ylab("Lambda")



Lambda.point

ggsave(filename = "Cost.lambda.point.plot.Ahud.pdf", path = cost.dir, Lambda.point, width = 8.5, height = 9.75,  units = "in", device = "pdf")


cost.plot.combined <- ggarrange(EPR.point, surv.plot.large, Lambda.point, ncol = 3, nrow = 1, common.legend = T)

ggsave(filename = "Cost_plot_combined_Ahud.pdf", path = cost.dir, cost.plot.combined, width = 28, height = 9,  units = "in", device = "pdf")



lambdaBoxplot <- ggplot(data = lambda.results.cost, aes(x = Food))+
  geom_boxplot(lwd = 1.1, aes(y = lambda, 
                              fill = factor(Treatment2), 
                              color = factor(Line), 
                              group = paste(Food, Order)),
               alpha = 0.3,
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  geom_point(data = lambda.list.mean, aes(y=mean, 
                                          color=factor(Treatment2),
                                          shape = factor(Line),
                                          group = paste(Food, Order)),
             size = 6, 
             #shape = factor(Line),
             position = position_dodge(width = 0.75),
  ) +
  
  #geom_errorbar(data = SurvTot.mean, aes(ymin=lower.ci, ymax=upper.ci), size=0.3, width = 0,
  #             position = position_dodge(width = 2))+
  
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  
  scale_fill_manual(values = c("cornflowerblue", "brown2"))+
  
  scale_color_manual(name = "Treatment",
                     values = c("cornflowerblue", "brown2"))+
  
  scale_shape_manual(values = c(18, 17))+
  theme_classic()+
  labs(y="lambda", x=expression ("Food Concentration ("*mu*"g C/L)"))+
  scale_x_discrete(breaks = c(600, 250))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  #ylim(1.402, 1.405)+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


lambdaBoxplot
ggsave(filename = "Malthusian.box.plot.wide_new_analysis_devtime.pdf", lambdaBoxplot, height = 8.5, width = 15, units = "in", device = "pdf")
