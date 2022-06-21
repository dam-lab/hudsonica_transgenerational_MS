library(data.table)
library(dplyr)
library(ggplot2)
library(lme4)
library(mgcv)
library(ggpubr)
library(tidyr)
library(sjPlot)
library(emmeans)


bs.data <- fread("Body_size_data_complete.txt")
bs.data$Generation <- as.numeric(bs.data$Generation)

## scalars taken from Durbin and Durbin, 1978, Limnology and Oceanography 23(5), 958-969
adult.weight.scalar <-12.37
adult.weight.power <- 3.6276

CI.weight.scalar <- 13.185
CI.weight.power <- 3.1858

bs.data <- bs.data %>% 
  mutate(Weight = case_when(Stage == "C1" ~ CI.weight.scalar*Length^(CI.weight.power),
                            Stage == "C6F" ~ adult.weight.scalar*Length^(adult.weight.power),
                            Stage == "C6M" ~ adult.weight.scalar*Length^(adult.weight.power)))





bs.data.adult <- bs.data %>% 
  rename(Rep = Replicate) %>% 
  filter(!Stage == "C1")

bs.data.cope <- bs.data %>% 
  rename(Rep = Replicate) %>% 
  filter(Stage == "C1")


bs.data.joined <- inner_join(bs.data.adult, bs.data.cope,
                             by = c("Generation",
                                    "Treatment",
                                    "Rep"))


bs.data.joined$weight.diff <- bs.data.joined$Weight.x - bs.data.joined$Weight.y


Dev.time.diff <- fread("Dev_time_diff.txt")


bs.data.dev.joined <- inner_join(bs.data.joined, Dev.time.diff,
                                 by= c("Generation",
                                       "Treatment",
                                       "Rep"))

bs.data.dev.joined$Growth.Rate <- bs.data.dev.joined$weight.diff / bs.data.dev.joined$Dev.diff

fwrite(bs.data.dev.joined, file = "BS_GR_data_all_w_11.txt")
bs.data.dev.joined <- fread("BS_GR_data_all_w_11.txt")

Growth.Rate.sum <- bs.data.dev.joined %>% 
  group_by(Generation, Treatment) %>% 
  summarise(mean = mean(Growth.Rate, na.rm = T),
            sd = sd(Growth.Rate, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)

##### Growth rate plots ######

text.size = 32
line.size = 2
point.size = 5

Growth.rate.plot.f11 <- ggplot(data = Growth.Rate.sum, aes(x = Generation, color = factor(Treatment), group = factor(Treatment))) +
  
  geom_point(aes(y=mean),
            size = 2,
            position = position_dodge(width = 0.5)
            ) +
  
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0,
                position = position_dodge(width = 0.5)) +
  
 # geom_smooth(data = bs.data.dev.joined, aes(x = Generation, y = Growth.Rate),
  #            method = y ~ poly(x,2),
   #           size = line.size)+
  
  
  scale_color_manual(name = "Treatment",
                     labels = c("AM","OA", "OW", "OWA"),
                     values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
                     )) +
  
  #scale_shape_manual(name = "Sex",
   #                  values = c(15,16,17),
    #                 labels = c("Female", "Male")) +
  
  scale_x_continuous(breaks = c(0,2,4,11))+
  
  labs(y = "Estimated Growth Rate (??g/day)") +
  
  theme_bw() +
  
  theme(axis.title = element_text(size = text.size,
                                  colour = "black"), 
        axis.text = element_text(size = text.size,
                                 colour = "black"),
        legend.text = element_text(size = text.size,
                                   colour = "black"),
        legend.title = element_text(size = text.size,
                                    colour = "black"))

Growth.rate.plot.f11

ggsave(filename = "Growth_rate_plot_f11.pdf", path = bs.dir, plot = Growth.rate.plot.f11, height = 101, width = 180, units = "mm")



## plotting body size data

bs.data.sum <- bs.data %>% 
  group_by(Generation, Treatment, Stage) %>% 
  summarise(mean = mean(Length, na.rm = TRUE),
            sd = sd(Length, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)


bs.data.sum.short <- bs.data.sum %>% 
  filter(Generation < 5)

fwrite(bs.data.sum, file = "Body_size_data_sum.txt")

bs.plot <- ggplot(data = bs.data.sum.short, aes(x = Generation, color = factor(Treatment), group = factor(Treatment))) +
  
  geom_point(aes(y=mean, 
                 shape = factor(Stage)
  ),
  size = 2,
  position = position_dodge(width = 1)) +
  
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0,
                position = position_dodge(width = 1)) +
  
  scale_color_manual(name = "Treatment",
                     labels = c("AM","OA", "OW", "OWA"),
                     values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
                     )) +
  
  scale_shape_manual(name = "Stage",
                     values = c(15,16,17)) +
  
  scale_x_continuous(breaks = c(0,2,4,11))+
  
  labs(y = "Body Size (mm)") +
  
  theme_bw() +
  
  theme(axis.title.x = element_text(size = text.size), 
        axis.text.x = element_text(size = text.size),
        axis.title.y = element_text(size = text.size),
        axis.text.y = element_text(size = text.size),
        legend.text = element_text(size = text.size),
        legend.title = element_text(size = text.size))

bs.plot

ggsave(filename = "Body_size_plot.pdf", path = bs.dir, plot = bs.plot, height = 101, width = 180, units = "mm")


bs.plot.f11 <- ggplot(data = bs.data.sum, aes(x = Generation, color = factor(Treatment), group = factor(Treatment))) +
  
  geom_point(aes(y=mean, 
                 shape = factor(Stage)
  ),
  size = 2,
  position = position_dodge(width = 1)) +
  
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0,
                position = position_dodge(width = 1)) +
  
  scale_color_manual(name = "Treatment",
                     labels = c("AM","OA", "OW", "OWA"),
                     values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
                     )) +
  
  scale_shape_manual(name = "Stage",
                     values = c(15,16,17)) +
  
  scale_x_continuous(breaks = c(0,2,4,11))+
  
  labs(y = "Body Size (mm)") +
  
  theme_bw() +
  
  theme(axis.title.x = element_text(size = text.size), 
        axis.text.x = element_text(size = text.size),
        axis.title.y = element_text(size = text.size),
        axis.text.y = element_text(size = text.size),
        legend.text = element_text(size = text.size),
        legend.title = element_text(size = text.size))

bs.plot.f11

ggsave(filename = "Body_size_plot_F11.pdf", path = bs.dir, plot = bs.plot.f11, height = 101, width = 180, units = "mm")




###### Statistics ######
setwd("C:/Users/james/Documents/Grad_school/OA_hudsonica/Body_size_pictures/A.hudsonica/")

body.size <- fread("Body_size_data_complete.txt")

body.size.no.f11 <- body.size %>% 
  filter(Generation < 10) %>% 
  mutate(Temperature = if_else(Treatment < 3, 13, 15))


bs.data$Generation.c <- as.numeric(bs.data$Generation)
gamC1 <- gam(Length ~ s(Generation.c, by = factor(Treatment), k=3), data = subset(bs.data, Stage=="C1"))
summary(gamC1)

gamC6M <- gam(Length ~ s(Generation.c, by = factor(Treatment), k=3), data = subset(bs.data, Stage=="C6M"))
summary(gamC6M)

gamC6F <- gam(Length ~ s(Generation.c, by = factor(Treatment), k=3), data = subset(bs.data, Stage=="C6F"))
summary(gamC6F)

m1 <- lmer(Length ~ Generation.c*Treatment*Stage + (1|Replicate), data = bs.data)
tab_model(m1)



##### Growth rate stats #####
growth.rate <- fread("BS_GR_data_all_w_11.txt")

is.even <- function(x) x %% 2 == 0


growth.rate.no.f11 <- growth.rate %>% 
  filter(Generation < 10) %>% 
  mutate(Temperature = if_else(Treatment < 3, 13, 15),
         pH = if_else(is.even(Treatment), 7.8, 8.2))

growth.rate.no.f11 <- unite(growth.rate.no.f11,#the data frame
                  Treat.Rep, #the name of the new column
                  c(Treatment, Rep), #the existing columns you want to combine
                  remove = FALSE) 

growth.rate.no.f11$Generation.c <- as.numeric(growth.rate.no.f11$Generation)

gam.gr.C6M <- gam(Growth.Rate ~ s(Generation.c, by = factor(Treatment), k=3), data = subset(growth.rate.no.f11, Stage.x=="C6M"))
summary(gam.gr.C6M)

gam.gr.C6F <- gam(Growth.Rate ~ s(Generation.c, by = factor(Treatment), k=3), data = subset(growth.rate.no.f11, Stage.x=="C6F"))
summary(gam.gr.C6F)

#m2 <- lmer(Growth.rate ~ Generation.c*Treatment*Stage + (1|Rep), data = growth.rate.no.f11)
#tab_model(m2)

#m3 <- lmer(Growth.rate ~ Generation.c*Treatment*Stage + (1|Rep), data = growth.rate.no.f11)
#tab_model(m3)

m4.M <- lmer(Growth.Rate ~ Generation*factor(Treatment) + (1|Treat.Rep), data = subset(growth.rate.no.f11, Stage.x == "C6M"))
tab_model(m4.M)

# model for pairwise
m4.M.pw <- lmer(Growth.Rate ~ factor(Generation)*factor(Treatment) + (1|Treat.Rep), data = subset(growth.rate.no.f11, Stage.x == "C6M"))
male.emmeans <- emmeans(m4.M.pw, pairwise ~ factor(Generation) | factor(Treatment))
pairs(male.emmeans)


m4.F <- lmer(Growth.Rate ~ Generation*factor(Treatment) + (1|Treat.Rep), data = subset(growth.rate.no.f11, Stage.x == "C6F"))
tab_model(m4.F)
m4.F.pw <- lmer(Growth.Rate ~ factor(Generation)*factor(Treatment) + (1|Treat.Rep), data = subset(growth.rate.no.f11, Stage.x == "C6F"))

female.emmeans <- emmeans(m4.F.pw, pairwise ~ factor(Generation) | factor(Treatment))
pairs(female.emmeans)



##### with F11 AM and OWA #####


growth.rate.1.4 <- growth.rate %>% 
  filter(Treatment == c(1,4))


growth.rate.1.4 <- unite(growth.rate.1.4,#the data frame
                            Treat.Rep, #the name of the new column
                            c(Treatment, Rep), #the existing columns you want to combine
                            remove = FALSE) 
growth.rate.1.4$Generation.c <- as.numeric(growth.rate.1.4$Generation)


gam.gr.C6M.1.4 <- gam(Growth.Rate ~ s(Generation.c, by = factor(Treatment), k=3), data = subset(growth.rate.1.4, Stage.x=="C6M"))
summary(gam.gr.C6M.1.4)

gam.gr.C6F.1.4 <- gam(Growth.Rate ~ s(Generation.c, by = factor(Treatment), k=3), data = subset(growth.rate.1.4, Stage.x=="C6F"))
summary(gam.gr.C6F.1.4)



m5.M <- lmer(Growth.Rate ~ Generation*factor(Treatment) + (1|Treat.Rep), data = subset(growth.rate.1.4, Stage.x == "C6M"))
tab_model(m5.M)
m5.M.pw <- lmer(Growth.Rate ~ factor(Generation)*factor(Treatment) + (1|Treat.Rep), data = subset(growth.rate.1.4, Stage.x == "C6M"))

male.emmeans.F11 <- emmeans(m5.M.pw, pairwise ~ factor(Generation) | factor(Treatment))
pairs(male.emmeans.F11)

male.emmeans.F11 <- emmeans(m5.M.pw, pairwise ~ factor(Treatment) | factor(Generation))
pairs(male.emmeans.F11)


m5.F <- lmer(Growth.Rate ~ Generation*factor(Treatment) + (1|Treat.Rep), data = subset(growth.rate.1.4, Stage.x == "C6F"))
tab_model(m5.F)

m5.F.pw <- lmer(Growth.Rate ~ factor(Generation)*factor(Treatment) + (1|Treat.Rep), data = subset(growth.rate.1.4, Stage.x == "C6F"))
female.emmeans.F11 <- emmeans(m5.F.pw, pairwise ~ factor(Generation) | factor(Treatment))
pairs(female.emmeans.F11)

female.emmeans.F11 <- emmeans(m5.F.pw, pairwise ~ factor(Treatment) | factor(Generation))
pairs(female.emmeans.F11)

