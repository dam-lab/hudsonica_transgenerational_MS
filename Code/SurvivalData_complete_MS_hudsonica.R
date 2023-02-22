
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(broom)
library(survminer)
library(survival)
library(survMisc)
library(sjPlot)
library(mgcv)
library(lme4)
library(itsadug)
library(car)
library(emmeans)



SurvTot <- fread("SurvDataFiles/Survival_data_total.txt")

##### Creating survival curves for each generation and treatment #####

### Original: Pop_growth_epr_and_survival.R
Surv.files <- filter(SurvTot, nx > 0)

Surv.files$nx <- as.numeric(Surv.files$nx)

Surv.files <- Surv.files[rep(seq(nrow(Surv.files)), 
                             Surv.files$nx),#the column which holds the value you wish to have everything repeated by
                         ]


Surv.files <- Surv.files %>%
  group_by(Generation, Treatment, Rep, Beak) %>%
  mutate(
    survival = case_when(
      # Food == 0 ~ 1,# all starved treatments did not survive
      nx == min(nx) ~ 0, #need to use the "==" and VALUE OF 1 EQUALS DEATH HAPPENING!!!!!!
      nx > min(nx) ~ 1 )) %>% 
  group_by(Generation, #have to add Generation when you are separating by generation!!!!
           Treatment, 
           Rep, 
           Beak) %>%
  mutate(keep = max(nx) - nx) %>% 
  group_by(Generation, 
           Treatment, 
           Rep, 
           Beak, 
           keep) %>% #need to group by unique and keep to get all the animals that died in between the start and the lowest amount of dead animals
  filter(survival == 1 & nx == min(nx) |
           survival == 0) %>% 
  filter(row_number() %in% 1:unique(keep) |
           survival == 0) %>% 
  select(-keep) %>% 
  ungroup()

#get rid of starting point since they all started alive

Surv.files <- Surv.files %>% 
  group_by(Generation, Treatment, Rep, Beak) %>% 
  filter(nx != max(nx))


#Surv.files <- Surv.files[!Surv.files$nx == 25,]
Surv.files <- as.data.frame(Surv.files)

# put the table in order
Surv.files$Generation <- as.numeric(Surv.files$Generation)
Surv.files <- Surv.files[order(Surv.files$Generation),]



Surv.files <- unite(Surv.files,#the data frame
                    unique, #the name of the new column
                    c(Generation, Treatment, Rep, Beak), #the existing columns you want to combine
                    remove = FALSE)



# must be a data frame for the grouping variable to be sorted appropriately
Surv.files <- as.data.frame(Surv.files)


#Surv.files.split <- split(Surv.files, f = Surv.files$Gen.Treat)




# create survival lists for all the data grouped by Gen, Treat, Rep, and Beak
Surv.object <- Surv.files %>%
  surv_group_by(grouping.vars = "unique") %>%
  separate(data = ., unique, into = c("Generation", "Treatment", "Rep", "Beak"))

Surv.object


# create the survival fit curves for statistical analysis 

Surv.fits <- surv_fit(Surv(time, #the time indicator
                           survival # the event indicator
)~Temp+pH, # the variables to test survival by
data = Surv.files, # the DATA FRAME with data you want to fit
group.by = "Generation") 

## survMisc package



# stats for each generation by temp and pH
Surv.stats <- surv_pvalue(Surv.fits, method = "log-rank", test.for.trend = TRUE)
Surv.stats


Surv.plots <-ggsurvplot(Surv.fits, 
                        color = "strata", 
                        conf.int = TRUE, 
                        palette = c("Green", # OA
                                    "Blue", # AM
                                    "Red", # OWA
                                    "Orange" # OW
                                    ), 
                        test.for.trend = TRUE,
                        pval = TRUE)

Surv.plots


##### Plots for mean survival across generations #####
SurvTot1 <- filter(SurvTot, SurvTot$lx > 0)
#SurvTot1 <- filter(SurvTot1, SurvTot1$Generation < 13)



SurvTot.agg <- aggregate(SurvTot1$lx, 
                         list(SurvTot1$Treatment, SurvTot1$Rep, SurvTot1$Beak, SurvTot1$Generation), #use a list so you don't have to use "interaction"
                         min)#use the minimum value for each of the individual beakers which should be the final survival

colnames(SurvTot.agg) <- c("Treatment", "Rep", "Beak", "Generation", "lx")

is.even <- function(x) x %% 2 == 0

SurvTot.agg$Generation.c <- as.numeric(SurvTot.agg$Generation)

SurvTot.mean <- SurvTot.agg %>%
  group_by(Treatment, Generation.c) %>%
  summarise(mean = mean(lx, na.rm = TRUE),
            sd = sd(lx, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)




pd=0.7
survPlotTotal <- ggplot(data = SurvTot.mean, aes(Generation.c, mean, color=factor(Treatment)))+
  geom_line(size=1,
            position = position_dodge(width = pd))+
  geom_point(size = 2, position = position_dodge(width = pd))+
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), size=0.3, width = 0,
                position = position_dodge(width = pd))+
  theme(legend.title = element_text(colour = "black", size=12))+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  )
  )+ 
  theme_classic()+
  labs(y="Survival", x="Generation")+
  scale_x_continuous(breaks = c(0,2,4))+
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00))+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  geom_text(x=15.5, y=0.51, size = 5, label = "*", color = "black")+
  geom_text(x=16, y=0.6, size = 5, label = "*", color = "black")+
  geom_text(x=26, y=0.3, size = 5, label = "*", color = "black")

survPlotTotal

ggsave(filename = paste0(surv.directory,"Survival_plot_F4.pdf"), plot = survPlotTotal, height = 101, width = 180, units = "mm")

##### Add F11 survival data #####

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
SurvData.agg$Generation.c <- as.numeric(SurvData.agg$Generation)

SurvTot.agg.F11 <- rbind(SurvTot.agg, SurvData.agg)
#SurvTot.agg <- filter(SurvTot.agg, SurvTot.agg$Treatment < 5)

SurvTot.agg.F11$Temp <- if_else(SurvTot.agg.F11$Treatment < 3, 13, 15)

SurvTot.agg.F11$pH <- if_else(is.even(SurvTot.agg.F11$Treatment), 7.8, 8.2)

Surv.lists <- list(SurvTot.agg, SurvTot.agg.F11)

for (i in length(1:Surv.lists)) {
  
  surv <- Surv.lists[[1]]
  
  surv$Generation.c <- as.numeric(surv$Generation)
  
  surv$Generation <- as.factor(as.numeric(surv$Generation))
  
  surv$Treatment <- as.factor(as.numeric(surv$Treatment))
  
  surv <- unite(surv,
                Treat.Rep,
                c(Treatment, Rep),
                remove = F)
  
}

SurvTot.agg.F11$Generation.c <- as.numeric(SurvTot.agg.F11$Generation)

SurvTot.agg.F11$Generation <- as.factor(as.numeric(SurvTot.agg.F11$Generation))

#SurvTot.agg.F11$Treatment <- as.factor(as.numeric(SurvTot.agg.F11$Treatment))


SurvTot.agg.F11 <- unite(SurvTot.agg.F11,
                     Treat.Rep,
                     c(Treatment, Rep),
                     remove = FALSE)

#SurvTot.agg$line <- if_else(SurvTot.agg$Treatment == 1 | SurvTot.agg$Treatment == 6, 
#                                          "AA", "HH")


SurvTot.mean.F11 <- SurvTot.agg.F11 %>%
  group_by(Treatment, Generation.c) %>%
  summarise(mean = mean(lx, na.rm = TRUE),
            sd = sd(lx, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)


fwrite(SurvTot.mean, file = paste(surv.directory,"Ahud_surv_mean.txt", sep = ""), sep = "\t")

SurvTot.mean[is.na(SurvTot.mean)] <- 0
SurvTot.mean.F11 <- SurvTot.mean.F11 %>% 
  filter(Treatment <5)

pd=0.7
text.size = 12
survPlotTotal.dashed <- ggplot()+
  geom_line(data = SurvTot.mean, aes(Generation.c, mean, color=factor(Treatment)),
            size=1,
            position = position_dodge(width = pd))+
  
  
  geom_point(data = SurvTot.mean.F11, aes(Generation.c, mean, color=factor(Treatment)),
                                     size = 2,
                                     position = position_dodge(width = pd))+
  
               
  geom_errorbar(data = SurvTot.mean.F11, aes(x = Generation.c, 
                                             ymin=lower.ci, 
                                             ymax=upper.ci, 
                                             color = factor(Treatment)), 
                size=0.3, 
                width = 0,
                position = position_dodge(width = pd))+
  
  geom_segment(aes(x=3.8, y=0.8311111, xend = 10.8, yend = 0.5523132), #AM dashed line
               color = "#0c7bdc", 
               linetype = "dashed",
               size = 1.3,
               alpha = 0.5)+
  
  geom_segment(aes(x=4.2, y=0.6458824, xend = 11.15, yend = 0.5555556), # OWA dashed line
               color = "#d55e00", 
               linetype = "dashed",
               size = 1.3,
               alpha = 0.5)+
  
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  )
  )+
  theme_classic()+
  labs(y="Survival", x="Generation")+
  scale_x_continuous(breaks = c(0,2,4,11))+
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00))+
  theme(legend.position = "none",
        axis.title.x = element_text(size = text.size), 
        axis.text.x = element_text(size = text.size),
        axis.title.y = element_text(size = text.size),
        axis.text.y = element_text(size = text.size))

survPlotTotal.dashed


ggsave(filename = "survplot_dashed_transparent.pdf", plot = survPlotTotal.dashed, path = surv.directory, width = 180, height = 101, units = "mm")



SurvTot.mean$Generation <- as.factor(as.numeric(SurvTot.mean$Generation.c))

SurvTot.agg.F11 <- SurvTot.agg.F11 %>% 
  filter(Treatment < 5)
survBoxplot <- ggplot(data = SurvTot.agg.F11, aes(x = factor(Generation)))+
  geom_boxplot(lwd = 1.1, aes(y = lx, fill = factor(Treatment)),
               alpha = 0.3,
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  geom_point(data = SurvTot.mean.F11, aes(factor(Generation.c), mean, color=factor(Treatment)),
             size = 2, 
             shape = 18,
             position = position_dodge(width = 0.75)
             ) +
  
  #geom_errorbar(data = SurvTot.mean, aes(ymin=lower.ci, ymax=upper.ci), size=0.3, width = 0,
   #             position = position_dodge(width = 2))+
  
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#0c7bdc", #AM
                               "#009e73", #OA
                               "#ffa500", #OW
                               "#d55e00" #OWA
  ))+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  ))+
  theme_classic()+
  labs(y="Survival", x="Generation")+
  scale_x_discrete(breaks = c(0,2,4, 11))+
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = text.size), 
        axis.text.x = element_text(size = text.size),
        axis.title.y = element_text(size = text.size),
        axis.text.y = element_text(size = text.size))


survBoxplot

ggsave(filename = "surv_boxplot_F11.pdf", plot = survBoxplot, path = surv.directory, height = 101, width = 180, units = "mm")




##### Statistics for Survival #####

## model for pairwise tukey comparisons
SurvTot.agg <- SurvTot.agg.F11 %>% 
  filter(Generation.c <11)
SurvTot.lm <- lmer(lx~Generation.c*Treatment+(1|Treat.Rep), 
                    REML = FALSE,
                    #family = gaussian,
                    data = SurvTot.agg)
summary(SurvTot.lm)

plot_model(SurvTot.lm, type = 'int')

tab_model(SurvTot.lm)

SurvTot.lm2 <- lmer(lx~factor(Generation)*factor(Treatment)+(1|Treat.Rep), 
                   REML = FALSE,
                   #family = gaussian,
                   data = SurvTot.agg)

surv.emm <- emmeans(SurvTot.lm2, pairwise ~ Generation | Treatment)

pairs(surv.emm, which=1, adjust = "none")

surv.emm2 <- emmeans(SurvTot.lm2, pairwise ~ Treatment | Generation)

pairs(surv.emm2, which = 1) # Tukey results
pairs(surv.emm2, which = 1, adjust = "none") # normal t-test comparisons
pairs(surv.emm2, which = 1, adjust = "bonf") # bonferonni corrections


Surv.pairwise <- tidy(pairwise.t.test(SurvTot.agg$lx, SurvTot.agg$Generation:SurvTot.agg$Treatment, p.adjust.method = "none"))
fwrite(Surv.pairwise, file = paste(surv.directory, "Statistics/Survival_pw.txt", sep = ""), sep = "\t")


## continuous model for a linear mixed model anova
SurvTot.lm.c <- lmer(lx~Generation.c*Treatment+(1|Treat.Rep), 
                   REML = FALSE,
                   #family = gaussian,
                   data = SurvTot.agg)

tab_model(SurvTot.lm.c)


## gam model
SurvTot.agg <- SurvTot.agg %>% 
  filter(Treatment < 5)


# s() indicates a smooth function
gam1 <- gam(lx ~ s(Generation.c, by = factor(Treatment), k = 3), data = SurvTot.agg)

plot(gam1, pages = 1)
plot_smooth(gam1, view="Generation.c", plot_all="Treatment", rug=FALSE)
#plot_smooth(gam1, view="Generation.c", plot_all="Beak", rug=FALSE)
gam.stats <- summary(gam1)




## three way anova
Surv.lm2 <- lmer(lx ~ Generation.c*Temp*pH +(1|Treat.Rep), data = SurvTot.agg)

Surv.3way.anova <- Anova(Surv.lm2)

Surv.3way.anova$factors <- rownames(Surv.3way.anova)

fwrite(Surv.3way.anova, file = paste(surv.directory, "Statistics/Survival_3way_anova.txt", sep = ""), sep = "\t")




##### stats with F11 #####

treats.1.4 <- c("1", "4")

SurvTot.agg.F11.1.4 <- SurvTot.agg.F11 %>% 
  filter(Treatment %in% treats.1.4)

SurvTot.lm.F11 <- lmer(lx~Generation.c*Treatment+(1|Treat.Rep), 
                   REML = FALSE,
                   #family = gaussian,
                   data = SurvTot.agg.F11.1.4)
summary(SurvTot.lm.F11)

plot_model(SurvTot.lm.F11, type = 'int')

tab_model(SurvTot.lm.F11)


SurvTot.lm.F11.2 <- lmer(lx~factor(Generation)*factor(Treatment)+(1|Treat.Rep), 
                    REML = FALSE,
                    #family = gaussian,
                    data = SurvTot.agg.F11.1.4)

surv.emm.f11 <- emmeans(SurvTot.lm.F11.2, pairwise ~ Generation | Treatment)

pairs(surv.emm.f11)

surv.emm2 <- emmeans(SurvTot.lm.F11.2, pairwise ~ Treatment | Generation)

pairs(surv.emm2, adjust = "none")


gam2 <- gam(lx ~ s(Generation.c, by = factor(Treatment), k = 3), data = SurvTot.agg.F11.1.4)

plot(gam2, pages = 1)
plot_smooth(gam2, view="Generation.c", plot_all="Treatment", rug=FALSE)
#plot_smooth(gam1, view="Generation.c", plot_all="Beak", rug=FALSE)
gam.stats2 <- summary(gam2)
gam.stats2




################################################################################################################################
########################################## Development Time ####################################################################


SurvTot.Cdev <- filter(SurvTot, SurvTot$Cdev > 0)


SurvTot.Cdev <- SurvTot.Cdev[rep(seq(nrow(SurvTot.Cdev)), SurvTot.Cdev$Cdev),]



SurvTot.Cdev <- SurvTot.Cdev %>% 
  select(,-c(Date,x,y)) %>% 
  relocate(Generation, .after = last_col())

Dev.time <- filter(SurvData, SurvData$Cdev > 0)
Dev.time <- Dev.time[rep(seq(nrow(Dev.time)), Dev.time$Cdev),]
Dev.time.600 <- Dev.time %>% 
  filter(Food == 600) %>% 
  filter(Treatment <5) %>% 
  select(,-c(unique,Food,Generation)) %>% 
  mutate(Generation =11)
  


Cdev.sum <- SurvTot.Cdev %>%
  group_by(Treatment, Generation) %>%
  summarise(mean = mean(time, na.rm = TRUE),
            sd = sd(time, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)



fwrite(Cdev.sum, file = paste(dev.directory, "Ahud_cdev_sum.txt", sep = ""), sep = "\t")




SurvTot.Cdev.f11 <- rbind(SurvTot.Cdev, Dev.time.600)
SurvTot.Cdev.f11$Generation.c <- as.numeric(SurvTot.Cdev.f11$Generation)
SurvTot.Cdev.f11$Generation <- as.factor(as.numeric(SurvTot.Cdev.f11$Generation.c))
fwrite(SurvTot.Cdev.f11, file = paste0(dev.directory, "Ahud_cdev_w_f11.txt"), sep = "\t")



Cdev.sum.f11 <- SurvTot.Cdev.f11 %>%
  group_by(Treatment, Generation.c) %>%
  summarise(mean = mean(time, na.rm = TRUE),
            sd = sd(time, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

Cdev.sum.f11$Generation <- as.factor(Cdev.sum.f11$Generation.c)

fwrite(Cdev.sum.f11, file = paste0(dev.directory, "Ahud_cdev_sum.txt"), sep = "\t")


CdevBoxplot <- ggplot(data = SurvTot.Cdev.f11, aes(x = factor(Generation)))+
  
  geom_boxplot(lwd = 1.1, aes(y = time, fill = factor(Treatment)),
               alpha = 0.3,
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  geom_point(data = Cdev.sum.f11, aes(Generation, mean, color=factor(Treatment)),
             size = 6, 
             shape = 18,
             position = position_dodge(width = 0.75)
  ) +
  
  #geom_errorbar(data = SurvTot.mean, aes(ymin=lower.ci, ymax=upper.ci), size=0.3, width = 0,
  #             position = position_dodge(width = 2))+
  
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#0c7bdc", #AM
                               "#009e73", #OA
                               "#ffa500", #OW
                               "#d55e00" #OWA
  ))+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  ))+
  theme_classic()+
  labs(y="Total Development time (days)", x="Generation")+
  scale_x_discrete(breaks = c(0,2,4,11))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  #ylim(1.402, 1.405)+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


CdevBoxplot


ggsave(filename = "Dev_time_boxplot_w_f11.pdf", plot = CdevBoxplot, path = dev.directory, width = 180, height = 101, units = "mm")



Cdev.sum$Generation.c <- as.numeric(Cdev.sum$Generation)

CdevPlotTotal <- ggplot(data = Cdev.sum, aes(Generation, mean, color=factor(Treatment)))+
  geom_line(size=1,
            position = position_dodge(width = pd))+
  geom_point(size = 2, position = position_dodge(width = pd))+
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), size=0.3, width=0, position = position_dodge(width = pd))+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  ))+
  theme_classic()+
  labs(y="Total Development Time (days)", x="Generation")+
  scale_x_continuous(breaks = c(0,2,4,11))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none",
        axis.title = element_text(size = text.size,
                                  color = "black"), 
        axis.text = element_text(size = text.size,
                                 color = "black"))

CdevPlotTotal

ggsave(filename = "Dev_time_plot_w_f11.pdf", plot = CdevPlotTotal, path = dev.directory, width = 180, height = 101, units = "mm")




##### F11 dashed plot #####


Cdev.sum.f11$Generation <- as.factor(as.numeric(Cdev.sum.f11$Generation.c))


CdevPlot.F11.dashed <- ggplot()+
  geom_line(data = subset(Cdev.sum.f11, Generation.c <5), aes(Generation.c, mean, color=factor(Treatment)),
            size=1,
            position = position_dodge(width = pd))+
  
  
  geom_point(data = Cdev.sum.f11, aes(Generation.c, mean, color=factor(Treatment)),
             size = 2,
             position = position_dodge(width = pd))+
  
  
  geom_errorbar(data = Cdev.sum.f11, aes(x = Generation.c, 
                                             ymin=lower.ci, 
                                             ymax=upper.ci, 
                                             color = factor(Treatment)), 
                size=0.3, 
                width = 0,
                position = position_dodge(width = pd))+
  
  geom_segment(aes(x=3.8, y=20.31551, xend = 10.8, yend = 20.75000), #AM dashed line
               color = "#0c7bdc", 
               linetype = "dashed",
               size = 1.3,
               alpha = 0.5)+
  
  geom_segment(aes(x=4.2, y=15.92357, xend = 11.1, yend = 17.11321), # OWA dashed line
               color = "#d55e00", 
               linetype = "dashed",
               size = 1.3,
               alpha = 0.5)+
  
  
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  ))+
  theme_classic()+
  labs(y="Total Development Time (days)", x="Generation")+
  scale_x_continuous(breaks = c(0,2,4,11))+
  #scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none",
        axis.title = element_text(size = text.size,
                                  color = "black"), 
        axis.text = element_text(size = text.size,
                                 color = "black"))

CdevPlot.F11.dashed

ggsave(filename = "Dev_time_plot_w_f11_dashed.pdf", plot = CdevPlot.F11.dashed, path = dev.directory, width = 180, height = 101, units = "mm")



#SurvTot.Cdev$Generation.c <- as.numeric(SurvTot.Cdev$Generation)

#SurvTot.Cdev$Generation <- as.factor(as.numeric(SurvTot.Cdev$Generation))



##### Statistics ######
SurvTot.Cdev$time <- as.numeric(SurvTot.Cdev$time)

SurvTot.Cdev$Treatment <- as.factor(SurvTot.Cdev$Treatment)

SurvTot.Cdev <- unite(SurvTot.Cdev,
                      Treat.Rep,
                      c(Treatment, Rep),
                      remove = FALSE)

SurvTot.Cdev$pH <- as.factor(as.numeric(SurvTot.Cdev$pH))
SurvTot.Cdev$Temp <- as.factor(as.numeric(SurvTot.Cdev$Temp))
SurvTot.Cdev$Generation.c <- as.numeric(SurvTot.Cdev$Generation)


## linear mixed models for pairwise comparisons
## have to use "generation" instead of continuous "generation.c" for pw comparisons

Cdev.lmm <- lmer(time ~ factor(Generation) * factor(Treatment) + (1|Treat.Rep),
                 REML = FALSE,
                 #family = gaussian,
                 data = SurvTot.Cdev)

tab_model(Cdev.lmm)


Cdev.emm <- emmeans(Cdev.lmm, pairwise ~ Generation | Treatment)

pairs(Cdev.emm, adjust = "none")

Cdev.emm2 <- emmeans(Cdev.lmm, pairwise ~ Treatment | Generation)

pairs(Cdev.emm2, adjust = "none")


Cdev.lmm2 <- lmer(time ~ Generation.c * factor(Treatment) + (1|Treat.Rep),
                 REML = FALSE,
                 #family = gaussian,
                 data = SurvTot.Cdev)



tab_model(Cdev.lmm2)

plot_model(Cdev.lmm2, 'int')


Cdev.pairwise <- tidy(pairwise.t.test(SurvTot.Cdev$time, SurvTot.Cdev$Generation:SurvTot.Cdev$Treatment, p.adjust.method = "none"))

Cdev.pairwise

#Cdev.pairwise <- filter(Cdev.pairwise, Cdev.pairwise$p.value < 0.05)
fwrite(Cdev.pairwise, file = paste(dev.directory,"Statistics/Dev_time_pairwise.txt", sep = ""), sep = "\t")

## model for three-way anova

Cdev.lm.2 <-  lmer(time~Generation.c*Temp*pH+(1|Treat.Rep), data = SurvTot.Cdev)

Cdev.anova.2 <- Anova(Cdev.lm.2)

Cdev.anova.2$factors <- rownames(Cdev.anova.2)

fwrite(Cdev.anova.2, file = paste(dev.directory,"Statistics/Dev_time_anova_ph_temp.txt", sep = ""), sep = "\t")

plot_model(Cdev.lm.2, 'int')

summary(Cdev.lm.2)
tab_model(Cdev.lm.2)


## gam model for analyzing development time across generations

Cdev.gam <- gam(time ~ s(Generation.c, by = factor(Treatment), k = 3), data = SurvTot.Cdev)

plot(Cdev.gam, pages = 1)
plot_smooth(Cdev.gam, view="Generation.c", plot_all="Treatment", rug=FALSE)
#plot_smooth(gam1, view="Generation.c", plot_all="Beak", rug=FALSE)
summary(Cdev.gam)



##### development with F11 ######
SurvTot.Cdev.f11 <- unite(SurvTot.Cdev.f11,
                      Treat.Rep,
                      c(Treatment, Rep),
                      remove = FALSE)

SurvTot.Cdev.f11.1 <- SurvTot.Cdev.f11 %>% 
  filter(Treatment == 1)
  
SurvTot.Cdev.f11.4 <- SurvTot.Cdev.f11 %>% 
  filter(Treatment == 4)


SurvTot.Cdev.f11.1.4 <- rbind(SurvTot.Cdev.f11.1, SurvTot.Cdev.f11.4)

SurvTot.Cdev.f11.1.4$pH <- as.factor(as.numeric(SurvTot.Cdev.f11.1.4$pH))
SurvTot.Cdev.f11.1.4$Temp <- as.factor(as.numeric(SurvTot.Cdev.f11.1.4$Temp))


## linear mixed models for pairwise comparisons
## have to use "generation" instead of continuous "generation.c" for pw comparisons

Cdev.lmm <- lmer(time ~ factor(Generation) * factor(Treatment) + (1|Treat.Rep),
                 REML = FALSE,
                 #family = gaussian,
                 data = SurvTot.Cdev.f11.1.4)

tab_model(Cdev.lmm)


Cdev.emm <- emmeans(Cdev.lmm, pairwise ~ Generation | Treatment)

pairs(Cdev.emm, adjust = "none")
pairs(Cdev.emm)

Cdev.emm2 <- emmeans(Cdev.lmm, pairwise ~ Treatment | Generation)

pairs(Cdev.emm2, adjust = "none")
pairs(Cdev.emm2)


Cdev.lmm2 <- lmer(time ~ Generation.c * factor(Treatment) + (1|Treat.Rep),
                  REML = FALSE,
                  #family = gaussian,
                  data = SurvTot.Cdev.f11.1.4)



tab_model(Cdev.lmm2)

plot_model(Cdev.lmm2, 'int')




Cdev.pairwise <- tidy(pairwise.t.test(SurvTot.Cdev.f11.1.4$time, SurvTot.Cdev.f11.1.4$Generation:SurvTot.Cdev.f11.1.4$Treatment, p.adjust.method = "none"))

Cdev.pairwise

#Cdev.pairwise <- filter(Cdev.pairwise, Cdev.pairwise$p.value < 0.05)
fwrite(Cdev.pairwise, file = paste(dev.directory,"Statistics/Dev_time_pairwise.txt", sep = ""), sep = "\t")

## model for three-way anova -- can't do this with just AM and OWA

#Cdev.lm.2 <-  lmer(time~Generation.c*Temp*pH+(1|Treat.Rep), data = SurvTot.Cdev.f11.1.4)

#Cdev.anova.2 <- Anova(Cdev.lm.2)

#Cdev.anova.2$factors <- rownames(Cdev.anova.2)

#fwrite(Cdev.anova.2, file = paste(dev.directory,"Statistics/Dev_time_anova_ph_temp.txt", sep = ""), sep = "\t")

#plot_model(Cdev.lm.2, 'int')

#summary(Cdev.lm.2)


## gam model for analyzing development time across generations

Cdev.gam <- gam(time ~ s(Generation.c, by = factor(Treatment), k = 3), data = SurvTot.Cdev.f11.1.4)

plot(Cdev.gam, pages = 1)
plot_smooth(Cdev.gam, view="Generation.c", plot_all="Treatment", rug=FALSE)
#plot_smooth(gam1, view="Generation.c", plot_all="Beak", rug=FALSE)
summary(Cdev.gam)




##########################################################################################################################################
######################################## Sex ratio #######################################################################################


## Extract the sex ratio data for scaling fecundity
SurvTot$M.Ratio <- as.numeric(SurvTot$M.Ratio)
SurvTot$F.Ratio <- as.numeric(SurvTot$F.Ratio)

SurvTotSex <- SurvTot %>%
  filter(F.Ratio >= 0 | M.Ratio >= 0) %>%
  group_by(Generation, Treatment, Rep, Beak) %>%
  summarise(F.Ratio = last(F.Ratio))

SurvTotSex$Generation.c <- as.numeric(SurvTotSex$Generation)

SurvTotSex$Generation <- as.factor(as.numeric(SurvTotSex$Generation))




gam.sexratio <- gam(F.Ratio ~ s(Generation.c, by = factor(Treatment), k = 3), data = SurvTotSex)
summary(gam.sexratio)


plot(gam.sexratio, pages = 1)
plot_smooth(gam.sexratio, view="Generation.c", plot_all="Treatment", rug=FALSE)
#plot_smooth(gam1, view="Generation.c", plot_all="Beak", rug=FALSE)
gam.stats.sexratio <- summary(gam.sexratio)

gam.stats.sexratio

# an HTML table of the results
tab_model(gam.sexratio)


SurvTotSex <- unite(SurvTotSex,#the data frame
                     Treat.Rep, #the name of the new column
                     c(Treatment, Rep), #the existing columns you want to combine
                     remove = FALSE)

SurvTotSex$Treatment <- as.factor(SurvTotSex$Treatment)


sexratio.lmm <- lmer(F.Ratio ~ Generation.c * Treatment + (1|Treat.Rep),
                     REML = FALSE,
                     #family = gaussian,
                     data = SurvTotSex)

summary(sexratio.lmm)
Anova(sexratio.lmm)

AIC(gam.sexratio, sexratio.lmm)

sexratio.emm <- emmeans(sexratio.lmm, pairwise ~ Generation | Treatment)



pairs(sexratio.emm)


sex.plot <- ggplot(data = SurvTotSex, aes(x = Generation.c,
                              y = F.Ratio/0.5,
                              color = factor(Treatment)))+
  
  geom_point() + 
  
  geom_smooth(method = "lm", se = F) +
  
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  ))+
  theme_classic()+
  labs(y="Female:Male", x="Generation")+
  scale_x_continuous(breaks = c(0,2,4))+

  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))
sex.plot
  
ggsave(filename = "Sex_ratio_plot.pdf", plot = sex.plot, path = sex.directory, height = 112, width = 180, units = "mm")

