
library(ggplot2)

library(data.table)

library(dplyr)

library(car)

library(emmeans)

#library(seacarb)

hudsonica.physical.data <- fread(file = "Physical_data_MS.txt")

Ahudsonica.pH.mean <- hudsonica.physical.data %>% 
  group_by(treatment) %>% 
  summarise(mean = mean(pH, na.rm = TRUE),
            sd = sd(pH, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)



Ahudsonica.temp.mean <- hudsonica.physical.data %>% 
  group_by(treatment) %>% 
  summarise(mean = mean(`Temperature (C)`, na.rm = TRUE),
            sd = sd(`Temperature (C)`, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)



fwrite(x = Ahudsonica.temp.mean, file = "A_hudsonica_temp_mean.txt", sep = "\t")
fwrite(x = Ahudsonica.pH.mean, file = "A_hudsonica_pH_mean.txt", sep = "\t")
#fwrite(x = Ahudsonica.pco2.mean, file = "A_hudsonica_pCO2_mean.txt", sep = "\t")


## filter out the cost treatments for transgen analysis

hudsonica.physical.data$treatment <- as.factor(hudsonica.physical.data$treatment)
hudsonica.physical.data$Date <- as.Date(hudsonica.physical.data$Date, "%m/%d/%Y")
#hudsonica.physical.data$pCO2 <- as.numeric(hudsonica.physical.data$pCO2)

treatments <- c("AM","OA","OW","OWA")

hudsonica.physical.data <- hudsonica.physical.data %>% 
  filter(treatment %in% treatments)


levels(hudsonica.physical.data$treatment)

## find the comparisons among different treatments for how different temperature, pH, and pco2 are 

temp.glm <- glm(Temperature ~ treatment, data = hudsonica.physical.data)

emm.temp <- emmeans(temp.glm, pairwise ~ treatment)

temp.contrasts <- as.data.frame(emm.temp$contrasts)

fwrite(temp.contrasts, file = "OA_temp_contrasts_w_cost.txt", sep = "\t")



pH.glm <- glm(pH ~ treatment, data = hudsonica.physical.data)

emm.pH <- emmeans(pH.glm, pairwise ~ treatment)

pH.contrasts <- as.data.frame(emm.pH$contrasts)

fwrite(pH.contrasts, file = "OA_pH_contrasts_w_cost.txt", sep = "\t")



Ahudsonica.temp.plot <- ggplot(data = hudsonica.physical.data, aes(Date, `Temperature (C)`, color = factor(treatment)))+
  geom_point(size = 1)+
  scale_y_continuous(name = "Temperature (°C)",
                     breaks = c(11,12,13,14,15,16,17), limits = c(10,17))+
  scale_x_date(date_labels = "%b/%Y")+
  scale_color_manual(name = NULL,
                     values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
                     ),
                     labels = c("Control",
                                expression("High C"*O[2]),
                                "High Temp",
                                expression("High Temp High C"*O[2])))+
  theme_light()+
  theme(legend.position = "bottom",
        legend.text.align = 0)
  
  
  
Ahudsonica.temp.plot



Ahudsonica.pH.plot <- ggplot(data = hudsonica.physical.data, aes(Date, pH, color = factor(treatment)))+
  geom_point(size = 1)+
  scale_y_continuous(name = "pH")+
  scale_x_date(date_labels = "%b/%Y")+
  scale_color_manual(name = NULL,
                     values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
                     ),
                     labels = c("Control",
                                expression("High C"*O[2]),
                                "High Temp",
                                expression("High Temp High C"*O[2])))+
  theme_light()+
  theme(legend.position = "bottom",
        legend.text.align = 0)

Ahudsonica.pH.plot


### Summarize the alkalinity data

Ahudsonica.pco2 <- fread("Alkalinity_data_total_leuker2000.txt")

pco2.mean.lueker2 <- Ahudsonica.pco2 %>% 
  group_by(Generation, `Target temp`, `Target pH`) %>% 
  summarise(mean = mean(`pCO2 in (matm)`, na.rm = TRUE),
            sd = sd(`pCO2 in (matm)`, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

pco2.mean.lueker3 <- Ahudsonica.pco2 %>% 
  group_by(`Target temp`, `Target pH`) %>% 
  summarise(mean = mean(`pCO2 in (matm)`, na.rm = TRUE),
            sd = sd(`pCO2 in (matm)`, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

TA.mean <- Ahudsonica.pco2 %>% 
  group_by(Generation, `Target temp`, `Target pH`) %>% 
  summarise(mean = mean(`TA in (mmol/kgSW)`, na.rm = TRUE),
            sd = sd(`TA in (mmol/kgSW)`, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

sal.mean<- Ahudsonica.pco2 %>% 
  group_by(Generation, `Target temp`, `Target pH`) %>% 
  summarise(mean = mean(Salinity, na.rm = TRUE),
            sd = sd(Salinity, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

fco2.mean <- Ahudsonica.pco2 %>% 
  group_by(Generation, `Target temp`, `Target pH`) %>% 
  summarise(mean = mean(`fCO2 in (matm)`, na.rm = TRUE),
            sd = sd(`fCO2 in (matm)`, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

omega.ca.mean <- Ahudsonica.pco2 %>% 
  group_by(Generation, `Target temp`, `Target pH`) %>% 
  summarise(mean = mean(`WCa in`, na.rm = TRUE),
            sd = sd(`WCa in`, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)


omega.ar.mean <- Ahudsonica.pco2 %>% 
  group_by(Generation, `Target temp`, `Target pH`) %>% 
  summarise(mean = mean(`WAr in`, na.rm = TRUE),
            sd = sd(`WAr in`, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

DIC.mean <- Ahudsonica.pco2 %>% 
  group_by(Generation, `Target temp`, `Target pH`) %>% 
  summarise(mean = mean(`TCO2 in (mmol/kgSW)`, na.rm = TRUE),
            sd = sd(`TCO2 in (mmol/kgSW)`, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)


Ahudsonica.pco2$Treatment <- as.factor(Ahudsonica.pco2$Treatment)

pco2.glm <- glm(`pCO2 in (matm)` ~ Treatment, data = Ahudsonica.pco2)

emm.pco2 <- emmeans(pco2.glm, pairwise ~ Treatment)

pco2.contrasts <- as.data.frame(emm.pco2$contrasts)

fwrite(pco2.contrasts, file = "OA_hudsonica_pco2_contrasts.txt", sep = "\t")



fwrite(pco2.mean.lueker3, file = "Ahudsonica_pCO2_mean_no_gen.txt", sep = "\t")
fwrite(pco2.mean.lueker2, file = "Ahudsonica_pCO2_mean.txt", sep = "\t")
fwrite(TA.mean, file = "Ahudsonica_TA_mean.txt", sep = "\t")
fwrite(fco2.mean, file = "Ahudsonica_fCO2_mean.txt", sep = "\t")
fwrite(omega.ca.mean, file = "Ahudsonica_omegaca_mean.txt", sep = "\t")
fwrite(omega.ar.mean, file = "Ahudsonica_omegaar_mean.txt", sep = "\t")
fwrite(DIC.mean, file = "Ahudsonica_DIC_mean.txt", sep = "\t")


Ahudsonica.pco2 <- Ahudsonica.pco2 %>% 
  mutate(group = case_when(treatment == 1 ~ 1,
                           treatment == 2 ~ 2,
                           treatment == 3 ~ 1,
                           treatment == 4 ~ 2))
  



co2.data.low <- Ahudsonica.pco2 %>% 
  filter(`Target pH` == 8.2)

co2.model.low <- lm(`pCO2 in (matm)` ~ Treatment, data = co2.data.low)
summary(co2.model.low)
Anova(co2.model.low)


co2.data.high <- Ahudsonica.pco2 %>% 
  filter(`Target pH` == 7.85)

co2.model.high <- lm(`pCO2 in (matm)` ~ Treatment, data = co2.data.high)
summary(co2.model.high)
Anova(co2.model.high)


Ahudsonica.pCO2.plot <- ggplot(data = hudsonica.physical.data, aes(Date, pCO2, color = factor(treatment)))+
  geom_point(size = 1)+
  scale_y_continuous(name = "pCO2")+
  scale_x_date(date_labels = "%b/%Y")+
  scale_color_manual(name = NULL,
                     values = c("blue", "forestgreen", "orange", "red"),
                     labels = c("Control",
                                expression("High C"*O[2]),
                                "High Temp",
                                expression("High Temp High C"*O[2])))+
  theme_light()+
  theme(legend.position = "bottom",
        legend.text.align = 0)

Ahudsonica.pCO2.plot
########################################################################################################################################
############################################ Reciprocal Transplant experiments ########################################################


hudsonica.physical.data <- fread(file = "Physical_data_MS.txt")
treatments <- c("AM->OWA","OWA->AM")

hudsonica.cost.physical.data <- hudsonica.physical.data %>% 
  filter(treatment %in% treatments)


Ahudsonica.pH.mean <- hudsonica.cost.physical.data %>% 
  group_by(treatment) %>% 
  summarise(mean = mean(pH, na.rm = TRUE),
            sd = sd(pH, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)


Ahudsonica.temp.mean <- hudsonica.cost.physical.data %>% 
  group_by(treatment) %>% 
  summarise(mean = mean(`Temperature (C)`, na.rm = TRUE),
            sd = sd(`Temperature (C)`, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)


fwrite(Ahudsonica.pH.mean, file = paste(cost.dir,"Cost_pH_mean.txt", sep = ""), sep = "\t")
fwrite(Ahudsonica.temp.mean, file = paste(cost.dir,"Cost_temp_mean.txt", sep = ""), sep = "\t")



co2.comp <- glm(pH ~ treatment, data = hudsonica.cost.physical.data)
summary(co2.comp)
Anova(co2.comp)

temp.comp <- glm(`Temperature (C)` ~ treatment, data = hudsonica.cost.physical.data)
summary(temp.comp)
Anova(temp.comp)
