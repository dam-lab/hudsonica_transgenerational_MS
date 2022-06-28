library(popbio)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(broom)
library(sjPlot)
library(car)
library(emmeans)



SurvData <- fread("Survival_data_total.txt")
# Create uniqe sorting column to group each technical and biological replicate

SurvData <- unite(SurvData,#the data frame
                  unique, #the name of the new column
                  c(Generation, Treatment, Rep, Beak), #the existing columns you want to combine
                  remove = FALSE) 

# change number of individuals to numeric format
SurvData$nx <- as.numeric(SurvData$nx)




SurvData$lx <- as.numeric(SurvData$lx)

## Extract the sex ratio data for scaling fecundity
SurvData$M.Ratio <- as.numeric(SurvData$M.Ratio)
SurvData$F.Ratio <- as.numeric(SurvData$F.Ratio)

SurvDataSex <- SurvData %>%
  filter(F.Ratio >= 0 | M.Ratio >= 0) %>%
  group_by(Generation, Treatment, Rep, Beak) %>%
  summarise(F.Ratio = last(F.Ratio))

SurvDataSex.Mean <- SurvDataSex %>%
  group_by(Generation, Treatment, Rep) %>%
  summarise(F.Ratio = mean(F.Ratio, na.rm = TRUE)) %>%
  as.data.frame(SurvDataSex.Mean)

SurvDataSex.Mean

SurvDataSex.Mean$Generation <- as.numeric(SurvDataSex.Mean$Generation)




# Create a vector formatted for day-specific survivorship as it fits along the off-diagonal of a Leslie Matrix
# See Caswell, H. 2001. Matrix Population Models for further details


SurvData1 <- SurvData %>%
  group_by(Treatment, Rep, Beak) %>%
  mutate(lx_diag = if_else(time == 0, 1, as.numeric(lx/lag(lx, default = first(lx))))) %>% ## always start with 100%
  mutate(lx_diag = if_else(lx_diag <= 1.000000, lx_diag, lag(lx_diag, n=1, default =last(lx_diag)))) %>%
  mutate(days = if_else(time == 0, 1, as.numeric(time-lag(time, default = first(time))))) # create a new column for the number of days spent at the respective survivorships


# Survivorship is reflective of prior day. Therefore, there can be no 0 survivorship in this vector.
# If all animals die, then the vector ends and the matrix is truncated at the appropriate time
SurvData1 <- filter(SurvData1, lx_diag > 0) 


# Check if there is any super survivorship. There can be no survivorship > 1.
# No individuals can be lost and reappear
if (any(SurvData1$lx_diag > 1) == TRUE ){
  SurvData1$lx_diag <- if_else(SurvData1$lx_diag > 1, lag(SurvData1$lx_diag), SurvData1$lx_diag)
}

any(SurvData1$lx_diag > 1)

SurvData1 <- as.data.frame(SurvData1)


# elongate the data frame to make it reflect actual days spent over the experiment. This essentially changes the matrix to a daily matrix.

SurvData1 <- SurvData1[rep(seq(nrow(SurvData1)), SurvData1$days),] 


SurvData1$Generation <- as.numeric(SurvData1$Generation)


#################################################################################################################################################################################
##### Development time data #####

Dev.time <- filter(SurvData, SurvData$Cdev > 0)
Dev.time <- Dev.time[rep(seq(nrow(Dev.time)), Dev.time$Cdev),]

Dev.time.sum <- Dev.time %>%
  group_by(Generation, Treatment, Rep, Beak) %>%
  summarise(mean = mean(time, na.rm = TRUE),
            sd = sd(time, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

Dev.time.sum <- unite(Dev.time.sum,#the data frame
                      unique, #the name of the new column
                      c(Generation, Treatment, Rep, Beak), #the existing columns you want to combine
                      remove = FALSE)


Dev.time.sum$Generation <- as.numeric(Dev.time.sum$Generation)

Dev.time.sum.short <- Dev.time.sum[,-c(1,7:11)]

## Add the dev.time to the survival table
SurvData1 <- inner_join(SurvData1,
                        Dev.time.sum.short,
                        by = c("Generation", "Treatment", "Rep", "Beak"))

SurvData1 <- SurvData1 %>% rename(dev.time = mean) # the new name of the column has to come first

SurvData1 <- unite(SurvData1,#the data frame
                  unique2, #the name of the new column
                  c(Generation, Treatment, Rep), #the existing columns you want to combine
                  remove = FALSE) 


Survival.list <- split(SurvData1, f = SurvData1$unique2) # don't create the list with a list of organizers. Then you create a complete matrix of samples which include some tables with no data



# lists within lists
Survival.list <- lapply(Survival.list, function (x) split(x, f = x$Beak))






#################################################################################################################################################################################
##### EPR data #####



EPR.data <- fread("Data_frames/EPR_HF_data_total.txt")

EPR.data$EPRtot[is.na(EPR.data$EPRtot)] = 0


EPR.data$Hftot[is.na(EPR.data$Hftot)] = 0

EPR.data$fecundity <- EPR.data$EPRtot*EPR.data$Hftot



EPR.data$Generation[EPR.data$Generation==10] <- 9

EPR.data$Treatment <- case_when(EPR.data$Treatment == "AA" ~ 1, #Ambient
                                EPR.data$Treatment == "AH" ~ 2, #OA
                                EPR.data$Treatment == "HA" ~ 3, #OW
                                EPR.data$Treatment == "HH" ~ 4) #OWA

EPR.data <- EPR.data %>% 
  mutate(Rep = case_when(Number < 13 ~ 1,
                          Number >=13 & Number < 25 ~ 2,
                          Number >= 25 ~ 3)
  )



EPR.data <- inner_join(EPR.data, #the data frame you want to modify
                       
                       SurvDataSex.Mean, # the data frame that has the column you want to add
                       
                       by = c("Generation", "Treatment", "Rep"))



# Caluclate per capita sex-specific fecundity

EPR.data$sex.spec.fecundity <- EPR.data$fecundity*EPR.data$F.Ratio

EPR.data <- unite(EPR.data,#the data frame
                  unique, #the name of the new column
                  c(Generation, Treatment, Rep), #the existing columns you want to combine
                  remove = FALSE) 

which(is.na(EPR.data$sex.spec.fecundity))

# Create a list of lists to use for matrix creations.


### NOTE: lists of survivorship and EPR data MUST have same names and same length of items

EPR.list <- split(EPR.data, f = EPR.data$unique)

# no fecundity data available
Survival.list[["3_1_4"]] <- NULL

names(Survival.list)
names(EPR.list)
names.surv <- names(Survival.list)
names.surv <- as.list(names.surv)
names.surv

#EPR.list[[1]]

# Create a dummy table to add data to

lambda.results <- data.frame(Variables = "dummy", lambda = 0, surv = 0, epr = 0)


for (i in 1:length(Survival.list)) { #list of names for each gen and treatment
  
  diag.ls <- Survival.list[[i]]
  
  epr.df <- EPR.list[[i]] # select the appropriate data frame to pull epr values from
  
  for (j in 1: length(diag.ls)) {
    
    diag.df <-  diag.ls[[j]] # temporary list used to index at each point
  
    diag.v <- diag.df$lx_diag # extract the vector of interest for the diagonal
    
    surv.value <- min(diag.df$lx) # find the corresponding survival value for plotting the fitness landscape
  
    leslie.matrix <- diag(diag.v) # create the matrix without the top row added
  
    dev.time.value <- as.integer(mean(diag.df$dev.time))

    zero <- matrix(0, nrow = 1, ncol = (dev.time.value))
  
  fecundity.vector <- epr.df$sex.spec.fecundity # select the vector with sex spec fecundity
  
  epr.vector <- epr.df$EPRtot
  
  epr.count <- as.integer(dim(leslie.matrix)[1]-dev.time.value)
  
  if (epr.count < 1) {
    
    epr.count <- 1
    
  } 
 
  
   for(k in 1:length(epr.vector)) {
    
    fecundity.value <- fecundity.vector[k] # use the sex specific fecundity values in sequence with the same survival matrix
    
    
    if (is.na(fecundity.value) == TRUE) {
      
      fecundity.value <- 0 ## use an arrow when assigning numbers, not a "=="
      
    } 
    
    epr.value <- epr.vector[k]
    
    fecundity.row <- t(c(zero, rep(fecundity.value, epr.count))) # combine the zero row and the epr.value for the matrix and transpose it to make it a row
    
    if (ncol(fecundity.row) > ncol(leslie.matrix)) { # if the dev.time is somehow greater than the survival matrix, then only make the last column reflective of epr
      
      delete <- ncol(fecundity.row)-ncol(leslie.matrix) # find the number of days that the dev time is greater than the survivorship
      
      fecundity.row <- fecundity.row[,-c(1:delete)] # delete those days from the fecundity row to make it the same size as the leslie matrix
      
    }
    
    leslie.matrix1 <- rbind(fecundity.row, leslie.matrix) # add the fecundity row to the matrix
    
    matrix.final <- leslie.matrix1[-nrow(leslie.matrix1),] # eliminate the last row to make it a square
    
    eigen.calcs <- eigen.analysis(matrix.final, zero = FALSE) # calculate the eigen values
    
    lambda.value <- eigen.calcs$lambda1 # extract the dominant eigen values
    
    lambda.row <- data.frame(Variables = names.surv[i], lambda = lambda.value, surv = surv.value, epr = epr.value) # create a 1x2 data frame to add to the end of the final data frame
    # data frame has to have the same colnames in order to rbind
    
    colnames(lambda.row) <- colnames(lambda.results) # make sure the data frames have the same names
    
    lambda.results <- bind_rows(lambda.results, lambda.row) # append the data frame with new results
    
    
  }
  
  
  
  }
}


lambda.results <- separate(lambda.results, "Variables", into = c("Generation", "Treatment", "Rep"))


# remove the first row as a last step
lambda.results <- lambda.results[-1,] 


fwrite(lambda.results, file = "lambda_results_devtime_surv_epr.txt")


# calculate malthusian parameter
lambda.results$Malthusian <- log10(lambda.results$lambda)







######################################################################################################################################################

####################### START STATISTICS WITH COMBINED LAMBDA RESULTS ################################################################################

######################################################################################################################################################





lambda.results <- fread("lambda_results_devtime_surv_epr_hf_sex_w_f11.txt")

lambda.results <- as.data.frame(lambda.results)

lambda.results <- lambda.results %>% 
  group_by(Generation, Treatment, Rep) %>% 
  mutate(lambda.stand = (lambda-mean(lambda))/sd(lambda),
         surv.stand = (surv-mean(surv))/sd(surv),
         epr.stand = (epr-mean(epr))/sd(epr),
         hf.stand = (hf-mean(hf))/sd(hf),
         sex.stand = (sex-mean(sex))/sd(sex),
         dev.stand = (dev.time-mean(dev.time))/sd(dev.time),
         lambda.rel = lambda/mean(lambda),
         surv.rel = surv/mean(surv),
         epr.rel = epr/mean(epr),
         hf.rel = hf/mean(hf),
         sex.rel = sex/mean(sex),
         dev.rel = dev.time/mean(dev.time))

lambda.results <- unite(lambda.results,#the data frame
                        unique, #the name of the new column
                        c(Generation, Treatment), #the existing columns you want to combine
                        remove = FALSE)

lambda.results <- unite(lambda.results,
                        Treat.Rep,
                        c(Treatment, Rep),
                        remove = FALSE)

lambda.results$Rep.c <- case_when(lambda.results$Treat.Rep == "1_1" ~ 1,
                                  lambda.results$Treat.Rep == "1_2" ~ 2,
                                  lambda.results$Treat.Rep == "1_3" ~ 3,
                                  lambda.results$Treat.Rep == "2_1" ~ 4,
                                  lambda.results$Treat.Rep == "2_2" ~ 5,
                                  lambda.results$Treat.Rep == "2_3" ~ 6,
                                  lambda.results$Treat.Rep == "3_1" ~ 7,
                                  lambda.results$Treat.Rep == "3_2" ~ 8,
                                  lambda.results$Treat.Rep == "3_3" ~ 9,
                                  lambda.results$Treat.Rep == "4_1" ~ 10,
                                  lambda.results$Treat.Rep == "4_2" ~ 11,
                                  lambda.results$Treat.Rep == "4_3" ~ 12)

# create continuous generation vector for anova and plotting
lambda.results$Generation.c <- as.numeric(as.character(lambda.results$Generation)) 
lambda.results$Generation <- as.factor(as.numeric(lambda.results$Generation))
lambda.results$Treatment <- as.factor(lambda.results$Treatment)
lambda.results$Rep.c <- as.numeric(lambda.results$Rep.c)

# check to see what the distribution looks like
hist(lambda.results$lambda)




# use a continuous generation model for the anova

malthusian.model <- lm(lambda~Generation.c*Treatment, data=lambda.results)


m <- Anova(malthusian.model)
m

m$factors <- rownames(m)
fwrite(m, file = "lambda_devtime_Anova.txt")

shapiro.test(malthusian.model$residuals)

hist(malthusian.model$residuals)


#leveneTest(malthusian.model)

foo <- plot_model(malthusian.model, 'int')

foo

# Pairwise comparison of malthusian parameters

# create the model to use for Tukey pairwise comparisons -- generation must be a factor, not a continuous variable
malthusian.model.2 <- lm(lambda~Generation*Treatment, data=lambda.results)

emmip(malthusian.model.2, Treatment ~ Generation)

emm_1 <- emmeans(malthusian.model.2, pairwise ~ Generation | Treatment)

p <- pairs(emm_1)
p

# model for writing the tukey results
malthusian.model.3 <- aov(lambda~unique, data = lambda.results)
p2 <- TukeyHSD(malthusian.model.3, "unique")
p2.tukey <- tidy(p2)

p3 <- tidy(pairwise.t.test(lambda.results$lambda, lambda.results$Generation:lambda.results$Treatment, p.adjust.method = "none"))



lambda.results$Temp <- case_when(lambda.results$Treatment == "1" ~ 13,
                                 lambda.results$Treatment == "2" ~ 13,
                                 lambda.results$Treatment == "3" ~ 15,
                                 lambda.results$Treatment == "4" ~ 15)


lambda.results$pH <- case_when(lambda.results$Treatment == "1" ~ 8.2,
                                 lambda.results$Treatment == "2" ~ 7.8,
                                 lambda.results$Treatment == "3" ~ 8.2,
                                 lambda.results$Treatment == "4" ~ 7.8)


malthusian.model.5 <- aov(lambda ~ Generation.c*as.factor(Temp)*as.factor(pH), data = lambda.results)

m5 <- Anova(malthusian.model.5)
m5
m5$factors <- rownames(m5)



fwrite(lambda.results, file = "lambda_results_devtime.txt")
fwrite(p2.tukey, file = "lambda_devtime_Tukey.txt")
fwrite(m5, file = "lambda_3way_anova.txt")
fwrite(p3, file = "lambda_pairwise.txt")




library(mgcv) # s() indicates a smooth function
gam1 <- gam(lambda ~ s(Generation.c, by = Treatment, k = 3), data = lambda.results)
summary(gam1)
#gam2 <- gam(lambda ~ s(Generation.c, by = Treatment, k = 4), data = lambda.results)
#gam3 <- gam(lambda ~ s(Generation.c, by = Treatment, k = 5), data = lambda.results)
#gam4 <- gam(lambda ~ s(Generation.c, by = Treatment, k = 6), data = lambda.results)
#gam5 <- gam(lambda ~ s(Generation.c, by = Treatment, k = 7), data = lambda.results)

AIC(gam1,gam2,gam3,gam4,gam5)

## try BIC also because that penalizes the more complex models more


library(itsadug)
library(sjPlot)

plot(gam1, pages = 1)
plot_smooth(gam1, view="Generation.c", plot_all="Treatment", rug=FALSE)
#plot_smooth(gam1, view="Generation.c", plot_all="Beak", rug=FALSE)
gam.stats <- summary(gam1)

gam.stats

# an HTML table of the results
tab_model(gam1)


## two-part mixed effects model

## use to identify how inflated the zeroes make the data

library(glmmTMB)

fit_zigauss <- glmmTMB(lambda~Generation.c*Treatment+(1|Rep.c),
                       data = lambda.results,
                       ziformula = ~Generation.c*Treatment+(1|Rep.c)) # this specifies the zero-inflated part of the model

## fixed effects results correspond to when response is >0, and zero-inflation results correspond to when lambda includes 0
fit_zigauss

summary(fit_zigauss)


## create a formatted output of model results with intra-class correlation (ICC)
tab_model(fit_zigauss)


fixef(fit_zigauss)$zi # the fixed-effects results of the zero-inflated model
ranef(fit_zigauss)$zi # the random effects intercepts of the model


## create linear mixed effects models that include zeroes and omit zeroes
library(lme4)
## create the predicted values graph
lambda.results1 <- lambda.results %>%
  filter(Generation.c <5) %>% 
  mutate(lambda_zero = if_else(lambda==0, 0, 1))



# model for predicting when lambda is either 0 or >0
gm1 <- glmer(lambda_zero ~ Generation.c*Treatment+(1|Rep.c),
             data = lambda.results1, family = binomial)
summary(gm1)
lambda.pred.prob.plot <- plot_model(gm1,type='int', colors = c("#0c7bdc", #AM
                                                               "#009e73", #OA
                                                               "#ffa500", #OW
                                                               "#d55e00" #OWA
)) + 
  theme_sjplot2()
lambda.pred.prob.plot

ggsave(filename = "Ahud_lambda_pred_prob_plot1.pdf", plot = lambda.pred.prob.plot, path = fit.directory, width = 6.5, height = 3.7, units = "in")

# model for predicting when lambda is not 0
lambda.nonzero = lambda.results1[lambda.results1$lambda>0,]

aa = lmer(lambda~Generation.c*Treatment+(1|Rep.c),
          data=lambda.nonzero)

set_theme(axis.textsize = 16)
lambda.pred.values.plot <- plot_model(aa,type='int', colors = c("#0c7bdc", #AM
                                                                "#009e73", #OA
                                                                "#ffa500", #OW
                                                                "#d55e00" #OWA
))+ 
  theme_sjplot2() 




lambda.pred.values.plot
ggsave(filename = "Ahud_lambda_pred_values_plot.pdf", plot = lambda.pred.values.plot, path = fit.directory, width = 16, height = 9, units = "in")

##### PATH analysis #####
mod <- 'lambda.rel ~ surv + epr + hf + dev.time + sex'

library(mnormt)
library(lavaan)
path.fit <- cfa(mod, data = lambda.results)

summary(path.fit)


## split by gen and treatment
mod2 <- 'lambda.rel ~ surv + epr + hf'

lambda.results <- unite(lambda.results,
                        unique,
                        c(Generation, Treatment),
                        remove = F)

# create a list of lambda data frames to create the SEMs
lambda.list <- split(lambda.results, f = lambda.results$unique)
path.list <- cfaList(mod2, lambda.list)
names(lambda.list)

# identify the coefficients for each factor within each SEM
coef(path.list)



##### multiple regression coefficients #####



library(lm.beta)



lambda.results <- unite(lambda.results,
                        unique,
                        c(Generation, Treatment),
                        remove = F)




lambda.list <- split(lambda.results, f = lambda.results$unique)

lm.list.2 <- lapply(lambda.list, function (x) lm(lambda.rel ~ surv + epr + hf, data = x))


beta.list <- lapply(lm.list.2, function (x) lm.beta(x))


summary(beta.list$`0_1`)

summary(beta.list$`0_2`)

summary(beta.list$`0_3`)

summary(beta.list$`0_4`)


summary(beta.list$`4_1`)
summary(beta.list$`4_2`)

summary(beta.list$`4_3`)
summary(beta.list$`4_4`)


summary(beta.list$`11_1`)

summary(beta.list$`11_4`)


Anova(beta.list$`0_1`,beta.list$`4_1`)
Anova(beta.list$`0_2`,beta.list$`4_2`)
Anova(beta.list$`0_3`,beta.list$`4_3`)
Anova(beta.list$`0_4`,beta.list$`4_4`)

Anova(beta.list$`0_1`,beta.list$`11_1`)
Anova(beta.list$`0_4`,beta.list$`11_4`)



###### lambda results with F11 ######

l11 <- lmer(lambda ~ as.factor(Generation) * as.factor(Treatment) + (1|Treat.Rep), data = lambda.results)

l11.emm <- emmeans(l11, pairwise ~ factor(Generation) | factor(Treatment))

l11.emm.Tukey <- pairs(l11.emm)
l11.emm.Tukey <- l11.emm.Tukey$emmeans
l11.emm.Tukey <- tidy(l11.emm.Tukey)

l11.emm.pw <- pairs(l11.emm, adjust = "none")
l11.emm.pw <- l11.emm.pw$emmeans
l11.emm.pw <- tidy(l11.emm.pw)


l11.emm2 <- emmeans(l11, pairwise ~ factor(Treatment) | factor(Generation))

l11.emm2.Tukey <- pairs(l11.emm2)
l11.emm2.Tukey <- l11.emm2.Tukey$emmeans
l11.emm2.Tukey <- tidy(l11.emm2.Tukey)

l11.emm2.pw <- pairs(l11.emm2, adjust = "none")
l11.emm2.pw <- l11.emm2.pw$emmeans
l11.emm2.pw <- tidy(l11.emm2.pw)



fwrite(l11.emm.Tukey, "Lambda_Tukey1.txt")
fwrite(l11.emm2.Tukey, "Lambda_Tukey2.txt")
fwrite(l11.emm.pw, "Lambda_contrasts1.txt")
fwrite(l11.emm2.pw, "Lambda_contrasts2.txt")



Lambda.11.sum <- lambda.results %>%
  group_by(Generation, Treatment) %>%
  summarise(mean = mean(lambda, na.rm = TRUE),
            sd = sd(lambda, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)


fwrite(Lambda.11.sum, "Lambda_mean_F11.txt")


################################################################################################################################################
##### Plots of lambda  #####

lambda.mean.c <- lambda.results %>%
  group_by(Generation.c, Treatment) %>%
  summarise(mean = mean(lambda, na.rm = TRUE),
            sd = sd(lambda, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)

pd=0.5
text.size = 12
lambdaPlotTotal <- ggplot()+
  geom_line(data = subset(lambda.mean.c, Generation.c < 5), aes(Generation.c, mean, color=factor(Treatment)),
            size=1,
            position = position_dodge(width = pd))+
  
  geom_point(data = lambda.mean.c, aes(Generation.c, mean, color=factor(Treatment)),
             size = 2, position = position_dodge(width = pd))+
  
  
  
  geom_errorbar(data = lambda.mean.c, aes(x=Generation.c,
                                          ymin=lower.ci, 
                                          ymax=upper.ci,
                                          color=factor(Treatment)), 
                size=0.3, width=0, position = position_dodge(width = pd))+
  
 
  
  geom_segment(aes(x=3.8, y=0.965, xend = 10.8, yend = 1.1393780), #AM dashed line
               color = "#0c7bdc", 
               linetype = "dashed",
               size = 1,
               alpha = 0.5)+
  
  geom_segment(aes(x=4.3, y=1.1428505, xend = 11.2, yend = 0.8867400), # OWA dashed line
               color = "#d55e00", 
               linetype = "dashed",
               size = 1,
               alpha = 0.5)+
  
  
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  ))+
  theme_classic()+
  labs(y=expression(Population~Fitness~(lambda)), x="Generation")+
  scale_x_continuous(breaks = c(0,2,4, 11))+
  scale_y_continuous(limits = c(0.4,1.2),
                     breaks = c(0.6,0.8,1.0,1.2))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = text.size), 
        axis.text.x = element_text(size = text.size),
        axis.title.y = element_text(size = text.size),
        axis.text.y = element_text(size = text.size))

lambdaPlotTotal

ggsave("lambda_plot_f11_dashed.pdf", plot = lambdaPlotTotal, path = fit.directory, width = 180, height = 101, units = "mm", dpi = 320)


## without F11
lambdaPlotTotal.no.f11 <- ggplot()+
  geom_line(data = subset(lambda.mean.c, Generation.c < 5), aes(Generation.c, mean, color=factor(Treatment)),
            size=1,
            position = position_dodge(width = pd))+
  
  geom_point(data = subset(lambda.mean.c, Generation.c < 5), aes(Generation.c, mean, color=factor(Treatment)),
             size = 2, position = position_dodge(width = pd))+
  
  
  
  geom_errorbar(data = subset(lambda.mean.c, Generation.c < 5), aes(x=Generation.c,
                                          ymin=lower.ci, 
                                          ymax=upper.ci,
                                          color=factor(Treatment)), 
                size=0.3, width=0, position = position_dodge(width = pd))+
  
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  ))+
  theme_classic()+
  labs(y=expression(Population~Fitness~(lambda)), x="Generation")+
  scale_x_continuous(breaks = c(0,2,4, 11))+
  scale_y_continuous(limits = c(0.4,1.2),
                     breaks = c(0.6,0.8,1.0,1.2))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title.x = element_text(size = text.size), 
        axis.text.x = element_text(size = text.size),
        axis.title.y = element_text(size = text.size),
        axis.text.y = element_text(size = text.size))

lambdaPlotTotal.no.f11

ggsave("lambda_plot.pdf", plot = lambdaPlotTotal.no.f11, path = fit.directory, width = 180, height = 101, units = "mm", dpi = 320)




fwrite(lambda.mean.c, file = "lambda_mean_devtime.txt")

# plot to see how data is distributed among treatments and reps
ggplot(data = lambda.results, aes(x = Generation.c,
                                  y = lambda,
                                  color = factor(Rep)))+
  geom_smooth(method = "lm", se = F)+
  
  geom_jitter(width = 0.7, alpha = 0.2)+
  
  geom_point()+
  
  facet_wrap(~Treatment)

# plot how many zeroes change across the experiment and how many positive lambda values change
lambda.results.zero <- lambda.results %>% 
  group_by(Generation.c, Treatment, Rep) %>% 
  filter(lambda == 0) %>% 
  summarise(zero_n = n())


ggplot(data = lambda.results.zero, aes(x = Generation.c,
                                       y = zero_n,
                                       color = factor(Rep)))+
  geom_smooth(method = "lm", se = F)+
  
  #geom_jitter(width = 0.7, alpha = 0.2)+
  
  geom_point()+
  
  facet_wrap(~Treatment)
  
## plot to see how many non-zero values for each replicate
lambda.results.pos <- lambda.results %>% 
  group_by(Generation.c, Treatment, Rep) %>% 
  filter(lambda > 0) %>% 
  summarise(zero_n = n())

ggplot(data = lambda.results.pos, aes(x = Generation.c,
                                       y = zero_n,
                                       color = factor(Rep)))+
  geom_smooth(method = "lm", se = F)+
  
  #geom_jitter(width = 0.7, alpha = 0.2)+
  
  geom_point()+
  
  facet_wrap(~Treatment)





## plots of traits vs. relative lambda


lambda.results.traits <- fread("lambda_results_devtime_surv_epr_hf_sex_w_f11.txt")

epr.lambda.plot <- ggplot(data = lambda.results.traits, aes(x = epr,
                                                            y = lambda
                                                            )) +
  geom_smooth(method = "lm", se = F) +
  
  #geom_point()+
  
  facet_wrap(~Treatment)
  
epr.lambda.plot


surv.lambda.plot <- ggplot(data = lambda.results.traits, aes(x = surv,
                                                             y = lambda,
                                                             color = factor(Generation))) +
  geom_smooth( se = F) +
  
  ylim(0,2)+
  #geom_point()+
  
  facet_wrap(~Treatment)

surv.lambda.plot




lambda.mean.c$Generation <- as.factor(lambda.mean.c$Generation.c)
lambdaBoxplot <- ggplot()+
  geom_boxplot(data = lambda.results, aes(Generation, lambda),
               lwd = 1.1, aes(fill = factor(Treatment)),
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  geom_point(data = lambda.mean.c, aes(Generation, mean, color = factor(Treatment)),
             size = 2)+
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
  labs(y="lambda", x="Generation")+
  scale_x_discrete(breaks = c(0,2,4,11))+
  #ylim(1.402, 1.405)+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


lambdaBoxplot

lambdaBoxplot <- ggplot(data = lambda.results, aes(x = Generation))+
  geom_boxplot(lwd = 1.1, aes(y = lambda, fill = factor(Treatment)),
               alpha = 0.3,
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  geom_point(data = lambda.mean.c, aes(Generation, mean, color=factor(Treatment)),
             size = 2, 
             shape = 18,
             position = position_dodge(width = 0.75)
  ) +
  
  
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
  labs(y="lambda", x="Generation")+  scale_x_discrete(breaks = c(0,2,4,11))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))


lambdaBoxplot
ggsave(filename = "Lambda_boxplot.pdf", plot = lambdaBoxplot, path = fit.directory, height = 101, width = 180, units = "mm")



