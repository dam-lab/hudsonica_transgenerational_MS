library(ggplot2)
library(data.table)
library(tidyr)
library(broom)
library(dplyr)
library(lme4)
library(mgcv)
library(emmeans)
library(car)


EPRtot <- fread("EPR_HF_data_total_w_11.txt")


EPRtot <- EPRtot %>% 
  mutate(unite(., Gen.treat, c(Generation, Treatment), remove = F),
         unite(., Treat.Rep, c(Treatment, Rep), remove = F))


eprStatsAll <- EPRtot %>% 
  group_by(Treatment, Generation) %>% 
  summarise(Rate = mean(EPRtot, na.rm = TRUE),
            sd = sd(EPRtot, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = Rate - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = Rate + qt(1 - (0.05 / 2), n.count-1)*se)

new.df <- data.frame(Treatment = c(2,3),
                     Generation = c(11,11),
                     Rate = c(0,0),
                     sd = c(0,0),
                     n.count = c(0,0),
                     se = c(0,0),
                     lower.ci = c(0,0),
                     upper.ci = c(0,0))

eprStatsAll <- rbind(eprStatsAll, new.df)



hfStatsAll <- EPRtot %>% 
  group_by(Treatment, Generation) %>% 
  summarise(HF = mean(Hftot, na.rm = TRUE),
            sd = sd(Hftot, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = HF - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = HF + qt(1 - (0.05 / 2), n.count-1)*se)

new.df <- data.frame(Treatment = c(2,3),
                      Generation = c(11,11),
                      HF = c(0,0),
                      sd = c(0,0),
                      n.count = c(0,0),
                      se = c(0,0),
                      lower.ci = c(0,0),
                      upper.ci = c(0,0))


hfStatsAll <- rbind(hfStatsAll, new.df)

eprStatsAll$Trait <- "EPR"
hfStatsAll$Trait <- "HS"




fwrite(eprStatsAll, file = "Ahud_epr_mean.txt")
fwrite(hfStatsAll, file = "Ahud_hf_mean.txt")


#eprStatsAll <- eprStatsAll[-1,]
#hfStatsAll <- hfStatsAll[-1,]


#####combine all plots with HF on the same plot#####

colnames(hfStatsAll) <- c("Treatment2", "Generation2", "HF", "sd2", "n.count2", "se2", "lower.ci2", "upper.ci2", "Trait.2")

eprSumComplete <- cbind(eprStatsAll, hfStatsAll)
eprSumComplete$Treatment <- as.factor(eprSumComplete$Treatment)


## have to make the HF values on the same scale as EPR values
eprSumComplete$HF <- eprSumComplete$HF*max(eprSumComplete$Rate)
eprSumComplete$lower.ci2 <- eprSumComplete$lower.ci2*max(eprSumComplete$Rate)
eprSumComplete$upper.ci2 <- eprSumComplete$upper.ci2*max(eprSumComplete$Rate)

# create NA values for plotting purposes
eprSumComplete$Rate[eprSumComplete$Rate==0] <- NA
eprSumComplete$HF[eprSumComplete$HF==0] <- NA
eprSumComplete$lower.ci[eprSumComplete$lower.ci==0] <- NA
eprSumComplete$upper.ci[eprSumComplete$upper.ci==0] <- NA
eprSumComplete$lower.ci2[eprSumComplete$lower.ci2==0] <- NA
eprSumComplete$upper.ci2[eprSumComplete$upper.ci2==0] <- NA

#eprSumComplete <- eprSumComplete[complete.cases(eprSumComplete),]

eprSumAA <- subset(eprSumComplete, Treatment == "1")
eprSumAH <- subset(eprSumComplete, Treatment == "2")
eprSumHA <- subset(eprSumComplete, Treatment == "3")
eprSumHH <- subset(eprSumComplete, Treatment == "4")

##### box plot #####
EPRtot$Generation.c <- as.numeric(EPRtot$Generation)
EPRtot$Generation <- as.factor(as.numeric(EPRtot$Generation.c))

eprStatsAll$Generation.c <- as.numeric(eprStatsAll$Generation)
eprStatsAll$Generation <- as.factor(as.numeric(eprStatsAll$Generation.c))

eprBoxplot <- ggplot(data = EPRtot, aes(x = factor(Generation)))+
  geom_boxplot(lwd = 1.1, aes(y = EPRtot, fill = factor(Treatment)),
               alpha = 0.3,
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  geom_point(data = eprStatsAll, aes(factor(Generation), Rate, color=factor(Treatment)),
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
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), x="Generation")+
  scale_x_discrete(breaks = c(0,2,4,11))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  #ylim(1.402, 1.405)+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))


eprBoxplot

ggsave(filename = "EPR_boxplot.pdf", path = epr.directory, plot = eprboxplot, height = 101, width = 180, units = "mm")


hfStatsAll$Generation.c <- as.numeric(hfStatsAll$Generation2)
hfStatsAll$Generation2 <- as.factor(as.numeric(hfStatsAll$Generation.c))


hfboxplot <- ggplot(data = EPRtot, aes(x = factor(Generation)))+
  geom_boxplot(lwd = 1.1, aes(y = Hftot, fill = factor(Treatment)),
               alpha = 0.3,
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  geom_point(data = hfStatsAll, aes(Generation2, HF, color=factor(Treatment2)),
             size = 2, 
             shape = 18,
             position = position_dodge(width = 0.75)
  ) +
  
  #geom_errorbar(data = SurvTot.mean, aes(ymin=lower.ci, ymax=upper.ci), size=0.3, width = 0,
  #             position = position_dodge(width = 2))+
  
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
  labs(y="Hatching Success", x="Generation")+
  scale_x_discrete(breaks = c(0,2,4,11))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))


hfboxplot

ggsave(filename = "HS_boxplot.pdf", path = epr.directory, plot = hfboxplot, height = 101, width = 180, units = "mm")


##### separate bar graphs for each treatment ##########

font.size = 12

eprSumAA.1 <- eprSumAA[,1:8]
eprSumAA.2 <- eprSumAA[,9:16]
eprSumAA.3 <- eprSumAA[1:3,9:16]


## AA
epr.bar.graphs.AA.no.f11 <- ggplot()+
  geom_bar(data = subset(eprSumAA, Generation < 5), aes(Generation, Rate, fill = factor(Treatment)),
           stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1)+
  geom_errorbar(data = subset(eprSumAA, Generation < 5), 
                aes(x = Generation, 
                    ymin=lower.ci, 
                    ymax=upper.ci), 
                size=0.3, 
                width = 1,
                position = position_dodge())+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#0c7bdc"))+
  theme_classic()+
  labs(y=expression(atop(Egg~Production~Rate,(female^-1~day^-1))), 
       x="Generation")+
  scale_x_continuous(breaks = c(0,2,4))+
  
  geom_line(data = subset(eprSumAA, Generation < 5), aes(x = Generation2, y=HF), 
            size=0.3,
            colour="purple")+
 # geom_segment(aes(x=4, y=18.83580, xend = 11, yend = 22.29463), 
  #             color = "purple", 
   #            linetype = "dashed",
    #           size = 1.3)+
  
  geom_errorbar(data = subset(eprSumAA, Generation < 5), aes(x = Generation2, ymin=lower.ci2, ymax=upper.ci2), 
                size=0.3, 
                width = 1,
                position = position_dodge(0.3),
                colour = "purple")+
  scale_y_continuous(limits = c(0,30),
                     breaks = c(0,5,10,15,20,25,30),
                     sec.axis = sec_axis(~./max(eprStatsAll$Rate), 
                                         name = "Hatching Success",
                                         breaks = c(0,0.25,0.5,0.75,1, 1.25)))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = font.size), 
        axis.text.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size),
        axis.text.y = element_text(size = font.size))


epr.bar.graphs.AA.no.f11

ggsave(filename = "EPR_HS_AM.pdf", plot = epr.bar.graphs.AA.no.f11, path = epr.directory, height = 65, width = 100, units = "mm")


## AA with F11
epr.bar.graphs.AA <- ggplot()+
  geom_bar(data = eprSumAA, aes(Generation, Rate, fill = factor(Treatment)),
           stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8)+
  geom_errorbar(data = eprSumAA, 
                aes(x = Generation, 
                    ymin=lower.ci, 
                    ymax=upper.ci), 
                size=0.3, 
                position = position_dodge())+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#0c7bdc"))+
  theme_classic()+
  labs(y=expression(atop(Egg~Production~Rate,(female^-1~day^-1))), 
       x="Generation")+
  scale_x_continuous(breaks = c(0,2,4,11))+
  
  geom_line(data = eprSumAA.3, aes(x = Generation2, y=HF), 
            size=0.3,
            colour="purple")+
  geom_segment(aes(x=4, y=18.83580, xend = 11, yend = 22.29463), 
               color = "purple", 
               linetype = "dashed",
               size = 0.3)+
  
  geom_errorbar(data = eprSumAA.2, aes(x = Generation2, ymin=lower.ci2, ymax=upper.ci2), 
                size=0.3, 
                position = position_dodge(0.3),
                colour = "purple")+
  scale_y_continuous(limits = c(0,30),
                     breaks = c(0,5,10,15,20,25,30),
                     sec.axis = sec_axis(~./max(eprStatsAll$Rate), 
                                         name = "Hatching Success",
                                         breaks = c(0,0.25,0.5,0.75,1, 1.25)))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = font.size), 
        axis.text.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size),
        axis.text.y = element_text(size = font.size))


epr.bar.graphs.AA

ggsave(filename = "EPR_HS_AM_F11.pdf", plot = epr.bar.graphs.AA, path = epr.directory, height = 65, width = 100, units = "mm")

## OA

epr.bar.graphs.AH <- ggplot(data = eprSumAH, aes(Generation, Rate, fill = factor(Treatment)))+
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8)+
  geom_errorbar(aes(ymin=lower.ci, 
                    ymax=upper.ci), 
                size=0.3, 
                width = 1.8,
                position = position_dodge())+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#009e73"))+
  theme_classic()+
  labs(y=expression(atop(Egg~Production~Rate,(female^-1~day^-1))), 
       x="Generation")+
  scale_x_continuous(breaks = c(0,2,4,11))+
  geom_line(aes(y=HF), 
            size=0.3,
            colour="purple")+
  geom_errorbar(aes(ymin=lower.ci2, ymax=upper.ci2), 
                size=0.3, 
                width =1.8,
                position = position_dodge(0.3),
                colour = "purple")+
  scale_y_continuous(limits = c(0,30),
                     breaks = c(0,5,10,15,20,25,30),
                     minor_breaks = waiver(),
                     sec.axis = sec_axis(~./max(eprStatsAll$Rate), 
                                         name = "Hatching Success",
                                         breaks = c(0,0.25,0.5,0.75,1, 1.25)))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = font.size, color = "black"), 
        axis.text.x = element_text(size = font.size, color = "black"),
        axis.title.y = element_text(size = font.size, color = "black"),
        axis.text.y = element_text(size = font.size, color = "black")) 

epr.bar.graphs.AH ## saved as 8.8 x 11.75

ggsave(filename = "EPR_HS_OA_2.pdf", plot = epr.bar.graphs.AH, path = epr.directory, height = 65, width = 100, units = "mm")

## OW
epr.bar.graphs.HA <- ggplot(data = eprSumHA, aes(Generation, Rate, fill = factor(Treatment)))+
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8)+
  geom_errorbar(aes(ymin=lower.ci, 
                    ymax=upper.ci), 
                size=0.3,
                width = 1.8,
                position = position_dodge())+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#ffa500"))+
  theme_classic()+
  labs(y=expression(atop(Egg~Production~Rate,(female^-1~day^-1))), 
       x="Generation")+
  scale_x_continuous(breaks = c(0,2,4,11))+
  geom_line(aes(y=HF), 
            size=0.3,
            colour="purple")+
  geom_errorbar(aes(ymin=lower.ci2, ymax=upper.ci2), 
                size=0.3, 
                width = 1.8,
                position = position_dodge(0.3),
                colour = "purple")+
  scale_y_continuous(limits = c(0,30),
                     breaks = c(0,5,10,15,20,25,30),
                     sec.axis = sec_axis(~./max(eprStatsAll$Rate), 
                                         name = "Hatching Success",
                                         breaks = c(0,0.25,0.5,0.75,1, 1.25)))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = font.size, color = "black"), 
        axis.text.x = element_text(size = font.size, color = "black"),
        axis.title.y = element_text(size = font.size, color = "black"),
        axis.text.y = element_text(size = font.size, color = "black"))


epr.bar.graphs.HA

ggsave(filename = "EPR_HS_OW2.pdf", path = epr.directory, plot = epr.bar.graphs.HA, height = 65, width = 100, units = "mm")## HH

## OWA
eprSumHH.1 <- eprSumHH[,10:18]
eprSumHH.2 <- eprSumHH[1:3,10:18]

## without F11 data
epr.bar.graphs.HH.no.f11 <- ggplot()+
  geom_bar(data = subset(eprSumHH, Generation < 5), aes(Generation, Rate, fill = factor(Treatment)),
           stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1)+
  geom_errorbar(data = subset(eprSumHH, Generation < 5), 
                aes(x = Generation, 
                    ymin=lower.ci, 
                    ymax=upper.ci), 
                size=0.3, 
                width = 1,
                position = position_dodge())+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d55e00"))+
  theme_classic()+
  labs(y=expression(atop(Egg~Production~Rate,(female^-1~day^-1))), 
       x="Generation")+
  scale_x_continuous(breaks = c(0,2,4))+
  
  geom_line(data = subset(eprSumHH, Generation < 5), aes(x = Generation2, y=HF), 
            size=0.3,
            colour="purple")+
  # geom_segment(aes(x=4, y=18.83580, xend = 11, yend = 22.29463), 
  #             color = "purple", 
  #            linetype = "dashed",
  #           size = 1.3)+
  
  geom_errorbar(data = subset(eprSumHH, Generation < 5), aes(x = Generation2, ymin=lower.ci2, ymax=upper.ci2), 
                size=0.3, 
                width = 1,
                position = position_dodge(0.3),
                colour = "purple")+
  scale_y_continuous(limits = c(0,30),
                     breaks = c(0,5,10,15,20,25,30),
                     sec.axis = sec_axis(~./max(eprStatsAll$Rate), 
                                         name = "Hatching Success",
                                         breaks = c(0,0.25,0.5,0.75,1, 1.25)))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = font.size), 
        axis.text.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size),
        axis.text.y = element_text(size = font.size))

epr.bar.graphs.HH.no.f11
ggsave(filename = "EPR_HS_OWA.pdf", path = epr.directory, plot = epr.bar.graphs.HH.no.f11, width = 100, height = 65, units = "mm")

## with F11 data
epr.bar.graphs.HH <- ggplot()+
  geom_bar(data = eprSumHH, aes(Generation, Rate, fill = factor(Treatment)),
           stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8)+
  geom_errorbar(data = eprSumHH, 
                aes(x= Generation,
                    ymin=lower.ci, 
                    ymax=upper.ci), 
                size=0.3, 
                position = position_dodge())+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d55e00"))+
  theme_classic()+
  labs(y=expression(atop(Egg~Production~Rate,(female^-1~day^-1))), 
       x="Generation")+
  scale_x_continuous(breaks = c(0,2,4,11))+
  geom_line(data = eprSumHH.2, aes(x=Generation2, y=HF), 
            size=0.3,
            colour="purple")+
  
  geom_errorbar(data = eprSumHH.1, aes(x=Generation2,
                                       ymin=lower.ci2, ymax=upper.ci2), 
                size=0.3, 
                position = position_dodge(0.3),
                colour = "purple")+
  
  geom_segment(aes(x=4, y=23.47346, xend = 11, yend = 17.49660), 
               color = "purple", 
               linetype = "dashed",
               size = 0.3,
               alpha = 0.5)+
  
  scale_y_continuous(limits = c(0,30),
                     breaks = c(0,10,20,30),
                     sec.axis = sec_axis(~./max(eprStatsAll$Rate), 
                                         name = "Hatching Success",
                                         breaks = c(0,0.25,0.5,0.75,1.00, 1.25)))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))


epr.bar.graphs.HH
ggsave(filename = "EPR_HS_OWA_F11.pdf", plot = epr.bar.graphs.HH, path = epr.directory, height = 65, width = 100, units = "mm")



epr.plot.combined <- ggarrange(epr.bar.graphs.AA.no.f11+rremove("x.title"), epr.bar.graphs.AH+rremove("axis.title"),  
                               epr.bar.graphs.HA, epr.bar.graphs.HH.no.f11+rremove("y.title"),
                               ncol = 2, nrow = 2)


epr.plot.combined

ggsave(filename = "EPR_plot_combined_thin.pdf", plot = epr.plot.combined, path = epr.directory, width = 180, height = 130.25, units = "mm")



EPRtot$fecundity <- EPRtot$EPRtot*EPRtot$Hftot


#eprStatsAll$Rate[eprStatsAll$Rate==0] <- NA
#hfStatsAll$HF[hfStatsAll$HF==0] <- NA

eprPlotTotal <- ggplot(data = eprStatsAll, aes(Generation, EPRtot, color=factor(Treatment)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=EPRtot-ci, ymax=EPRtot+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), x="Generation")+
  scale_x_continuous(breaks = c(0,2,4))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

eprPlotTotal

hfPlotTotal <- ggplot(data = hfStatsAll, aes(Generation2, HF, color=factor(Treatment2)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=HF-ci2, ymax=HF+ci2), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y="Hatching Frequency", x="Generation")+
  scale_x_continuous(breaks = c(0,2,4))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


hfPlotTotal 



##### Fitness landscapes for epr #####



lambda.results <- fread("lambda_results_devtime_surv_epr_hf_sex_stand_rel.txt")
EPRtot.0 <- filter(EPRtot, Generation == 0)

last.gen <- max(lambda.results$Generation.c)

EPRtot.last <- filter(EPRtot, Generation == last.gen)

EPRtot.0.last <- rbind(EPRtot.0, EPRtot.last)


EPRtot.0.last$Treatment[EPRtot.0.last$Treatment=="AA"] <- 1
EPRtot.0.last$Treatment[EPRtot.0.last$Treatment=="AH"] <- 2
EPRtot.0.last$Treatment[EPRtot.0.last$Treatment=="HA"] <- 3
EPRtot.0.last$Treatment[EPRtot.0.last$Treatment=="HH"] <- 4


#EPRtot.0.last.HH <- filter(EPRtot.0.last, EPRtot.0.last$Treatment == 4)


lambda.results.0.full <- filter(lambda.results, Generation == 0)
lambda.results.last.full <- filter(lambda.results, Generation.c == last.gen)
lambda.results.0.last.full <- rbind(lambda.results.0.full, lambda.results.last.full)

lambda.results.0.last.full$Treatment2 <- case_when(lambda.results.0.last.full$Treatment == 1 ~ "AM",
                                      lambda.results.0.last.full$Treatment == 2 ~ "OA",
                                      lambda.results.0.last.full$Treatment == 3 ~ "OW",
                                      lambda.results.0.last.full$Treatment == 4 ~ "OWA")

# change the order of the factor to be in the order you want. Help from: https://stackoverflow.com/questions/5490638/how-to-change-the-order-of-facet-labels-in-ggplot-custom-facet-wrap-labels
lambda.results.0.last.full <- within(lambda.results.0.last.full, 
                        Treatment2 <- factor(Treatment2, 
                                             
                                             # put the specific levels in the order you want
                                             
                                             levels = c("AM", # 1
                                                        "OA", # 2
                                                        "OW", # 3
                                                        "OWA") # 4 
                                             
                        ))

scale <- 0.0166666666666666666666666666666666666667



## Panels are wrapped by Treatment
x.pos <- 62.5
y.pos <- 0.02


EPRtot.0.last.mean <- EPRtot.0.last %>% 
  group_by(Treatment, Generation) %>% 
  summarise(Rate = mean(EPRtot, na.rm = TRUE),
            sd = sd(EPRtot, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = Rate - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = Rate + qt(1 - (0.05 / 2), n.count-1)*se)

## calculate pearson correlation coefficient

lambda.results.0.last.full <- unite(lambda.results.0.last.full,
                Gen.treat,
                c(Generation, Treatment),
                remove = FALSE)


gp = group_by(lambda.results.0.last.full, Gen.treat)

summarise(gp, cor(epr, lambda))


ggplot()+
  geom_density(data = lambda.results.0.last.full, aes(x=epr, 
                                       y=after_stat(count),
                                       group = factor(Generation), 
                                       
                                       fill = factor(Generation)), alpha = 0.7, binwidth = 0.01)+
  
  #geom_freqpoly(data = EPRtot.0.25, aes(x=EPRtot,
   #                                     y=after_stat(count),
    #                                    group = factor(Generation),
     #                                   fill = factor(Generation)), alpha = 0.7, binwidth = 10)+
  
  ## add the points for where the mean epr value shows up
  
  #geom_point(data = EPRtot.0.25.mean, aes(x = Rate,
  #                                   y = y.pos,
   #                                  color = factor(Generation)),
    #         position = position_dodge(width = 15),
     #        size = 7) +
  
  
  #geom_point(data = lambda.mean.0.25, aes(x=x.pos, # x.pos needs to be 62.5
                                          
   #                                       y = mean/scale, # keep this as divided when just using the points and no landscape with scale = 25 for density or 0.6666665 for count
                                          
    #                                      color = factor(Generation)),
             
     #        position = position_dodge(width = 15),
      #       alpha = 0.7,
       #      size = 7)+
  
  #geom_errorbar(data = lambda.mean.0.25, aes(x=x.pos, # x.pos needs to be 62.5
                                             
   #                                         ymin=lower.ci/scale,
    #                                         ymax=upper.ci/scale, # keep this as divided when just using the points and no landscape with scale = 25 for density or 0.6666665 for count
                                            
     #                                       color = factor(Generation)),
      #          position = position_dodge(width = 15),
       #         size = 1.0,
        #        width = 8)+
  
  geom_smooth(data = lambda.results.0.last.full, aes(x = epr, y = lambda.rel*2, color = factor(Generation), group = factor(Generation)),
              method = "lm",# make sure we know that there is only 2 possible data values
              se = TRUE
              )+
  
  #geom_point(data = lambda.results.0.last.full, aes(x=epr, y = lambda.rel, color = factor(Generation), group = factor(Generation)))+
  
  
  
  theme_minimal()+
  scale_fill_manual(values = c("black", "red"),
                    
                    labels = c("F0", 
                               bquote(F*.(last.gen))) ## reference an object in the legend. Help from: https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
                    )+
  
  scale_color_manual(guide = "none",
                     values = c("black", "red"))+
  
  
  scale_x_continuous(name = expression(Egg~Production~Rate~(female^-1~day^-1)))+
  
  scale_y_continuous(labels = function (x) x*10,
                     sec.axis = dup_axis(trans = ~./25, # keep this as multiplied when just using the points and no landscape with scale = 25 for density or 0.6666665 for count
                                         
                                         name = expression(Relative~Fitness)))+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_blank(),
        legend.position = "bottom",
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.text = element_text(size = 20))+
  facet_wrap(~Treatment2)



## try it with Hatching frequency
HF.0.last.mean <- EPRtot.0.last %>% 
  group_by(Treatment, Generation) %>% 
  summarise(Freq = mean(Hftot, na.rm = TRUE),
            sd = sd(Hftot, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = Freq - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = Freq + qt(1 - (0.05 / 2), n.count-1)*se)



ggplot()+
  geom_density(data = lambda.results.0.last.full, aes(x=hf, 
                                                      y=after_stat(count)/25,
                                       group = factor(Generation), 
                                       
                                       fill = factor(Generation)), 
               alpha = 0.7, 
               binwidth = 0.01)+
  
  ## add the points for where the mean epr value shows up
  
  #geom_point(data = EPRtot.0.25.mean, aes(x = Rate,
  #                                   y = y.pos,
  #                                  color = factor(Generation)),
  #         position = position_dodge(width = 15),
  #        size = 7) +
  
  
  #geom_point(data = lambda.mean.0.25, aes(x=x.pos, # x.pos needs to be 62.5

#                                       y = mean/scale, # keep this as divided when just using the points and no landscape with scale = 25

#                                      color = factor(Generation)),

#        position = position_dodge(width = 15),
#       alpha = 0.7,
#      size = 7)+

#geom_errorbar(data = lambda.mean.0.25, aes(x=x.pos, # x.pos needs to be 62.5

#                                         ymin=lower.ci/scale,
#                                         ymax=upper.ci/scale, # keep this as divided when just using the points and no landscape with scale = 25

#                                       color = factor(Generation)),
#          position = position_dodge(width = 15),
#         size = 1.0,
#        width = 8)+

  geom_smooth(data = lambda.results.0.last.full, aes(x = hf, y = lambda.rel*20, color = factor(Generation), group = factor(Generation)),
           method = "lm"
           )+
  
  
  
  theme_minimal()+
  scale_fill_manual(values = c("black", "red"),
                    
                    labels = c("F0", 
                               bquote(F*.(last.gen))),
  )+
  
  scale_color_manual(guide = "none",
                     values = c("black", "red"))+
  
  
  scale_x_continuous(name = "Hatching Frequency")+
  
  

  
  scale_y_continuous(limits = c(0,40),
                     
                     name = "count",
                     
                     sec.axis = dup_axis(trans = ~./20, # keep this as multiplied when just using the points and no landscape with scale = 25
                                         
                                         name = "Relative Fitness"))+
  
  
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_blank(),
        legend.position = "bottom",
        legend.text.align = 0,
        legend.text = element_text(size = 20),
        legend.title = element_blank())+
  facet_wrap(~Treatment2)#+ylim(0,1000)







### SAVE THE HF FIGURE AS 6.6 x 16 INCHES #########################################################################################################

ggplot()+
  geom_smooth(data = lambda.results.0.last.full, aes(x = hf, y = lambda.rel, color = factor(Generation), group = factor(Generation)),
              method = "lm")+
  facet_wrap(~Treatment)

## Panels are wrapped by Generation
ggplot()+
  geom_density(data = EPRtot.0.last, aes(x=EPRtot, 
                                         #group = Generation, 
                                         fill = Treatment), alpha = 0.3, binwidth = 0.01)+
  geom_point(data = lambda.mean.0.last, aes(x=Treatment*scale, y = mean/scale, color = factor(Treatment)), 
             #alpha = 0.7,
             size = 7)+
  geom_errorbar(data = lambda.mean.0.last, aes(x=Treatment*scale,
                                               ymin=lower.ci/scale,
                                               ymax=upper.ci/scale,
                                               color = factor(Treatment)),
                size = 1.0,
                width = 8,
                alpha = 0.7)+
  #geom_line(data = lambda.mean.0.25, aes(x=Treatment.scale, y = mean.scale, color = factor(Treatment), group = Treatment))+
  theme_minimal()+
  scale_fill_manual(values = c("blue", "forestgreen", "orange", "red"),
                    labels = c("Control", 
                               expression("Effect of C"*O[2]), 
                               "Effect of Temperature", 
                               expression("Effect of Temperature and C"*O[2])))+
  scale_color_manual(guide = "none",
                     values = c("blue", "forestgreen", "orange", "red"))+
  
  
  scale_x_continuous(name = expression(Egg~Production~Rate~(female^-1~day^-1)))+
  
  scale_y_continuous(sec.axis = dup_axis(trans = ~.*scale, 
                                         name = expression(Mean~Absolute~Fitness~(generation^-1))))+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_blank(),
        legend.position = "bottom",
        legend.text.align = 0)+
  facet_wrap(~Generation)



###############################################################################################################################################################
##### Stats for EPR and HF #####

## Notes: gams can't be used as a three-way interaction. 
## To do a pairwise comparison, you need to create linear models to feed into emmeans.


library(dplyr)
library(broom)
library(car)
library(sjPlot)

EPRtot <- fread("EPR_HF_data_total_w_11.txt")

#create linear models that are tested against each other grouped by generation and analyzes EPR/HF by temp, pH, and temp*pH
EPRtot <- EPRtot[!is.na(EPRtot),]
EPRtot <- EPRtot[!is.na(HFtot),]
EPRtot$Generation[EPRtot$Generation==10] <-9


#EPRtot <- na.omit(EPRtot)
EPRtot$Generation.c <- as.numeric(EPRtot$Generation)
EPRtot$Generation <- as.factor(as.numeric(EPRtot$Generation))
#EPRtot$Treatment <- as.factor(EPRtot$Treatment)


EPRtot <- unite(EPRtot,
                Treat.Rep,
                c(Treatment, Rep),
                remove = FALSE)

EPRtot <- unite(EPRtot,
                Gen.treat,
                c(Generation, Treatment),
                remove = FALSE)


# remove for where there is no data
EPRtot <- subset(EPRtot, Gen.treat != "9_HA")

## create gam models to test if epr and hf change across generations

library(mgcv) # s() indicates a smooth function
gam1 <- gam(EPRtot ~ s(Generation.c, by = factor(Treatment), k = 3), data = subset(EPRtot, Generation.c < 5))

epr.gam.stats <- summary(gam1)
epr.gam.stats

tab_model(gam1)

AA_HH <- c(1,4)

EPRtot2 <- EPRtot %>% 
  filter(Treatment == AA_HH)
gam2 <- gam(EPRtot ~ s(Generation.c, by = factor(Treatment), k = 3), data = EPRtot2)

summary(gam2)
tab_model(gam2)




## generalized additive mixed models with replicates as mixed effects
mm <- gamm(EPRtot ~ s(Generation.c, by = Treatment, k = 3), random = list(Treat.Rep=~1), data = EPRtot)
summary(mm$lme)
summary(mm$gam)

#tab_model(mm$gam)

## gam for HF
gam3 <- gam(Hftot ~ s(Generation.c, by = factor(Treatment), k = 3), data = subset(EPRtot, Generation.c < 5))
hf.gam.stats <- summary(gam3)
hf.gam.stats

tab_model(gam3)


gam4 <- gam(Hftot ~ s(Generation.c, by = factor(Treatment), k = 3), data = EPRtot2)

summary(gam4)

# plot the gam models
library(itsadug)

gam.list <- list(gam1, gam2)

lapply(gam.list, function (x) plot_smooth(x, view = "Generation.c", plot_all = "Treatment", rug = FALSE))



#plot(gam1, pages = 1)
#plot_smooth(gam1, view="Generation.c", plot_all="Treatment", rug=FALSE)
#plot_smooth(gam1, view="Generation.c", plot_all="Beak", rug=FALSE)

#plot_smooth(mm$gam, view="Generation.c", plot_all="Treatment", rug=FALSE)


#plot(gam2, pages = 1)
#plot_smooth(gam2, view="Generation.c", plot_all="Treatment", rug=FALSE)

#anova.gam(gam2)






# create models for Tukey pairwise comparisons
library(emmeans)
library(lme4)
l <- lmer(EPRtot~Generation.c*factor(Treatment)+(1|Treat.Rep), 
          #REML = FALSE,
          #family = gaussian,
          data = subset(EPRtot, Generation.c <5))
l
summary(l)
Anova(l)
ranef(l)

plot_model(l, type = 'int')
tab_model(l)
Anova(l)




l2 <- lmer(EPRtot~factor(Generation)*factor(Treatment)+(1|Treat.Rep),
           data = subset(EPRtot, Generation.c <5))

tab_model(l2)

emmip(l2, Treatment ~ Generation)

epr.emm <- emmeans(l2, pairwise ~ Generation | Treatment)

pairs(epr.emm)


l3 <- lmer(Hftot~Generation.c*factor(Treatment)+(1|Treat.Rep), 
            #REML = FALSE,
            #family = gaussian,
            data = subset(EPRtot, Generation.c <5))

plot_model(l3, type = 'int')
tab_model(l3)
Anova(l3)


l4 <- lmer(Hftot~factor(Generation)*factor(Treatment)+(1|Treat.Rep), 
           #REML = FALSE,
           #family = gaussian,
           data = subset(EPRtot, Generation.c <5))


hf.emm <- emmeans(l4, pairwise ~ Generation | Treatment)

pairs(hf.emm)


### Models for AA and HH with F11
l5 <- lmer(Hftot~Generation.c*factor(Treatment)+(1|Treat.Rep), 
           #REML = FALSE,
           #family = gaussian,
           data = EPRtot2)


tab_model(l5)
Anova(l5)

l6 <- lmer(Hftot~factor(Generation)*factor(Treatment)+(1|Treat.Rep), 
           #REML = FALSE,
           #family = gaussian,
           data = EPRtot2)


hf.emm.2 <- emmeans(l6, pairwise ~ Treatment | Generation)
pairs(hf.emm.2)


tab_model(l6)
Anova(l6)


l7 <- lmer(EPRtot~Generation.c*factor(Treatment)+(1|Treat.Rep), 
           #REML = FALSE,
           #family = gaussian,
           data = EPRtot2)


tab_model(l7)
Anova(l7)

l8 <- lmer(EPRtot~factor(Generation)*factor(Treatment)+(1|Treat.Rep), 
           #REML = FALSE,
           #family = gaussian,
           data = EPRtot2)

epr.emm.2 <- emmeans(l8, pairwise ~ Treatment | Generation)
pairs(epr.emm.2)

tab_model(l8)
Anova(l8)



# create an object with a dataframe of the estimated marginal means and the groups for each factor
epr.tukey.groups <- CLD(epr.emm)
hf.tukey.groups <- CLD(hf.emm)


library(agricolae)

l9 <- lm(EPRtot ~ Gen.treat, data = subset(EPRtot, Generation.c < 5))
lsd <- LSD.test(l9, "Gen.treat")
epr.pw.groups <- lsd$groups
epr.pw.groups$variables <- rownames(epr.pw.groups)



l10 <- lm(Hftot ~ Gen.treat, data = subset(EPRtot, Generation.c < 5))
lsd2 <- LSD.test(l10, "Gen.treat")
hf.pw.groups <- lsd2$groups
hf.pw.groups$variables <- rownames(hf.pw.groups)


l11 <- lm(EPRtot ~ Gen.treat, data = EPRtot)
lsd3 <- LSD.test(l11, "Gen.treat")
epr.pw.groups.f11 <- lsd3$groups
epr.pw.groups.f11$variables <- rownames(epr.pw.groups.f11)


l12 <- lm(Hftot ~ Gen.treat, data = EPRtot)
lsd4 <- LSD.test(l12, "Gen.treat")
hf.pw.groups.f11 <- lsd4$groups
hf.pw.groups.f11$variables <- rownames(hf.pw.groups.f11)




fwrite(epr.tukey.groups, file = "Statistics/EPR_tukey_groups.txt")
fwrite(hf.tukey.groups, file = "Statistics/HF_tukey_groups.txt")

fwrite(epr.pw.groups, file = "Statistics/EPR_pw_groups.txt")
fwrite(hf.pw.groups, file = "Statistics/HF_pw_groups.txt")

fwrite(epr.pw.groups.f11, file = "Statistics/EPR_pw_groups_F11.txt")
fwrite(hf.pw.groups.f11, file = "Statistics/HF_pw_groups_F11.txt")


# pairwise comparisons with grouping variables
epr.pairwise.all <- tidy(pairwise.t.test(EPRtot$EPRtot, EPRtot$Generation:EPRtot$Treatment, p.adj = "none"))

epr.pairwise.all <- filter(epr.pairwise.all, epr.pairwise.all$p.value < 0.05)

hf.pairwise.all <- tidy(pairwise.t.test(EPRtot$HFtot, EPRtot$Generation:EPRtot$Treatment, p.adj = "none"))

hf.pairwise.all <- filter(hf.pairwise.all, hf.pairwise.all$p.value < 0.05)

fwrite(epr.pairwise.all, file = paste(epr.directory,"Statistics/EPR.pairwise.all.txt", sep = ""), sep = "\t")

fwrite(hf.pairwise.all, file = paste(epr.directory,"Statistics/hf.pairwise.all.txt", sep = ""), sep = "\t")




EPR.lm2 <- lm(EPRtot ~ Generation*Treatment, data = EPRtot)
emm_EPR <- emmeans(EPR.lm2, pairwise ~ Generation | Treatment)

p <- pairs(emm_EPR)
p
epr.groups <- CLD(emm_EPR)


HF.lm2 <- lm(HFtot ~ Generation*Treatment, data = EPRtot)

emm_HF <- emmeans(HF.lm2, pairwise ~ Generation | Treatment)


p2 <- pairs(emm_HF)
p2

hf.groups <- CLD(emm_HF)

fwrite(epr.groups, file = "Statistics/EPR_groups.txt")
fwrite(hf.groups, file = "Statistics/HF_groups.txt")



EPRtot$Temp <- as.factor(as.numeric(EPRtot$Temp))
EPRtot$pH <- as.factor(as.numeric(EPRtot$pH))


library(car)
## 3 way anova

epr.3.way <- aov(EPRtot~Generation.c*Temp*pH, data = EPRtot)
summary(epr.3.way)
Three.way.anova.epr <- Anova(epr.3.way)
Three.way.anova.epr$factors <- rownames(Three.way.anova.epr)
fwrite(Three.way.anova.epr, file = "Three_way_anova_epr.txt")

hf.3.way <- aov(HFtot~Generation*Temp*pH, data = EPRtot)
summary(hf.3.way)
Three.way.anova.hf <- Anova(hf.3.way)
Three.way.anova.hf$factors <- rownames(Three.way.anova.hf)
fwrite(Three.way.anova.hf, file = "Three_way_anova_hf.txt")

