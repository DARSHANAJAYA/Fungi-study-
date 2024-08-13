# loading the libraries
rm(list = ls(all.names = TRUE))
library(tidyverse)
library(readxl)
library(glmmTMB)
library(emmeans)
library(gee)
library(gamlss)
library(gamlss.cens)
library(survival)
library(ggalt)
library(glmtoolbox)
library(DHARMa)
library(ggeffects) # version 1.2.0




#Loading the data 
fungi_data <- read_excel("/Share/localshare/djayakumari/Mixed marginal study/fungi.xlsx")


#Changing the type of the covariates 
fungi_data$Material <- as.factor(fungi_data$Material)
fungi_data$Revestimento <- as.factor(fungi_data$Revestimento)
fungi_data$Face <- as.factor(fungi_data$Face)
fungi_data$Rep <- as.factor(fungi_data$Rep)


#Transformation of the response variable by imputing the value 10 by 0.0001 and 0 by 0.95.
#this is done to avoid the boundary values of 0 and 1 for modelling using the beta distribution. 

fungi_data_tr <- fungi_data %>% mutate(percent_mold_growth= recode(coloniza, '10' = '0.0001',
                                                                   '9' = '0.05', '8'= '0.15',
                                                                   '7'= '0.25','6'= '0.35',
                                                                   '5'= '0.45','4'= '0.55',
                                                                   '3' = '0.65', '2'= '0.75',
                                                                   '1'= '0.85', '0'= '0.95'))
#transformation for interval response variable for the start of the interval and end of the interval
fungi_data_tr <- fungi_data_tr %>% mutate(start= recode(coloniza, '10' = '0.0001',
                                                        '9' = '.01', '8'= '.11',
                                                        '7'= '0.21','6'= '0.31',
                                                        '5'= '0.41','4'= '0.51',
                                                        '3' = '0.61', '2'= '0.71',
                                                        '1'= '0.81', '0'= '0.9'))
fungi_data_tr <- fungi_data_tr %>% mutate(end = recode(coloniza, '10' = '0.0001',
                                                       '9' = '0.1', '8'= '0.2',
                                                       '7'= '0.3','6'= '0.4',
                                                       '5'= '0.5','4'= '0.6',
                                                       '3' = '0.7', '2'= '0.8',
                                                       '1'= '0.9', '0'= '0.9999'))

##changing the name of the covariate Revestimento to english 

fungi_data_tr <- fungi_data_tr %>% mutate(Revestimento = recode(Revestimento, 'Com' = 'With Coating',
                                                                'Sem' = 'No Coating'))

##changing the name of the other covariates to english names 
names(fungi_data_tr)[2] <- "coating"
names(fungi_data_tr)[3] <- "replications"
names(fungi_data_tr)[4] <- "day"

#changing the response variable type to numeric types
percent_mold_growth <- as.numeric(fungi_data_tr$percent_mold_growth)
fungi_data_tr$start <- as.numeric(fungi_data_tr$start)
fungi_data_tr$end <- as.numeric(fungi_data_tr$end)



# changing the name of the response variable to the transformed response variable 
fungi_data_tr$percentage_mold_growth_tr <- percent_mold_growth




#Exploratory plot

fungi_data_tr %>%
  ggplot(aes(x = day, col = replications)) +
  theme_bw() +
  geom_linerange(aes(ymin = start * 100, ymax = end * 100, xmin = day, xmax = day),
                 position = position_dodge(2)) +
  geom_point(aes(y = end * 100), size = .5, position = position_dodge(2)) +
  geom_point(aes(y = start * 100), size = .5, position = position_dodge(2)) +
  facet_grid(Material~coating+Face)+
  ylab("% mold growth")+
  xlab("Week") +
  scale_x_continuous(breaks = c(7,14,21,28), labels = 1:4) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  scale_color_discrete(name = "Replicate")

# code to save the plot 
ggsave(filename =  "fungi_exploratory.jpeg", width =9, height = 5)


### Mixed modelling using 

fit_mix_gamlss <- gamlss(percentage_mold_growth_tr ~ day * coating * Material+random(replications:Face),
                         data=fungi_data_tr, family = BE)

###Diagnostic plot 
wp(fit_mix_gamlss)


##Other combinations of the covariates and checking the model selection using drop1 method

drop1(fit_mixed1, test = "Chisq")

fit_mixed2 <- glmmTMB(percentage_mold_growth_tr ~ (day + coating + Material)^2 + (1 | replications:Face),
                      family = "beta_family", 
                      data = fungi_data_tr)
drop1(fit_mixed2, test = "Chisq")
fit_mixed3 <- glmmTMB(percentage_mold_growth_tr ~ day + coating * Material + (1 | replications:Face),
                      family = "beta_family", 
                      data = fungi_data_tr)
fit_mixed4 <- glmmTMB(percentage_mold_growth_tr ~ day * coating + Material + (1 | replications:Face),
                      family = "beta_family", 
                      data = fungi_data_tr)
fit_mixed5 <- glmmTMB(percentage_mold_growth_tr ~ day * Material + coating + (1 | replications:Face),
                      family = "beta_family", 
                      data = fungi_data_tr)
fit_mixed6 <- glmmTMB(percentage_mold_growth_tr ~ day + coating + Material + (1 | replications:Face),
                      family = "beta_family", 
                      data = fungi_data_tr)
fit_mixed7 <- glmmTMB(percentage_mold_growth_tr ~ day + coating + (1 | replications:Face),
                      family = "beta_family", 
                      data = fungi_data_tr)
fit_mixed8 <- glmmTMB(percentage_mold_growth_tr ~ coating + Material + (1 | replications:Face),
                      family = "beta_family", 
                      data = fungi_data_tr)
fit_mixed9 <- glmmTMB(percentage_mold_growth_tr ~ day + Material + (1 | replications:Face),
                      family = "beta_family", 
                      data = fungi_data_tr)
fit_mixed10 <- glmmTMB(percentage_mold_growth_tr ~ day + (1 | replications:Face),
                       family = "beta_family", 
                       data = fungi_data_tr)
fit_mixed11 <- glmmTMB(percentage_mold_growth_tr ~ coating + (1 | replications:Face),
                       family = "beta_family", 
                       data = fungi_data_tr)
fit_mixed12 <- glmmTMB(percentage_mold_growth_tr ~ Material + (1 | replications:Face),
                       family = "beta_family", 
                       data = fungi_data_tr)
drop1(fit_mixed12, test = "Chisq")# material not significant 
drop1(fit_mixed11, test = "Chisq")#coating is significant 
drop1(fit_mixed10, test = "Chisq")#day is significant 
drop1(fit_mixed9, test = "Chisq")# day is significant 
drop1(fit_mixed8, test = "Chisq")# coating is significant, material is not 
drop1(fit_mixed7, test = "Chisq")##day and coating significant 
drop1(fit_mixed6, test = "Chisq")###check this 


# calculating the AIC and BIC values 
AIC_mixed <- c(AIC(fit_mixed1), AIC(fit_mixed2), AIC(fit_mixed3), AIC(fit_mixed4),
               AIC(fit_mixed5), AIC(fit_mixed6), AIC(fit_mixed7), AIC(fit_mixed8),
               AIC(fit_mixed9), AIC(fit_mixed10), AIC(fit_mixed11), AIC(fit_mixed12))
AIC_mixed

BIC_mixed <- c(BIC(fit_mixed1), BIC(fit_mixed2), BIC(fit_mixed3), BIC(fit_mixed4),
               BIC(fit_mixed5), BIC(fit_mixed6), BIC(fit_mixed7), BIC(fit_mixed8),
               BIC(fit_mixed9), BIC(fit_mixed10), BIC(fit_mixed11), BIC(fit_mixed12))
BIC_mixed



## Marginal modelling using gee 

fit_marginal_1 <-  glmgee(percentage_mold_growth_tr ~ day, data = fungi_data_tr,id = replications:Face, 
                          family = binomial, corstr = "Exchangeable")
fit_marginal_2 <-  glmgee(percentage_mold_growth_tr ~ coating, data = fungi_data_tr,id = replications:Face, 
                          family = binomial, corstr = "Exchangeable")
fit_marginal_3 <-  glmgee(percentage_mold_growth_tr ~ Material, data = fungi_data_tr,id = replications:Face, 
                          family = binomial, corstr = "Exchangeable")
fit_marginal_4 <-  glmgee(percentage_mold_growth_tr ~ day + coating , data = fungi_data_tr,id = replications:Face, 
                          family = binomial, corstr = "Exchangeable")
fit_marginal_5 <- glmgee(percentage_mold_growth_tr ~ coating + Material , data = fungi_data_tr,id = replications:Face, 
                         family = binomial, corstr = "Exchangeable") 
fit_marginal_6 <- glmgee(percentage_mold_growth_tr ~ day+ Material , data = fungi_data_tr,id = replications:Face, 
                         family = binomial, corstr = "Exchangeable") 
fit_marginal_7 <- glmgee(percentage_mold_growth_tr ~ day + coating + Material , data = fungi_data_tr,id = replications:Face, 
                         family = binomial, corstr = "Exchangeable") 
fit_marginal_8 <- glmgee(percentage_mold_growth_tr ~ (day+ coating + Material)^2 , data = fungi_data_tr,id = replications:Face, 
                         family = binomial, corstr = "Exchangeable") 
fit_marginal_9 <-  glmgee(percentage_mold_growth_tr ~ day*coating+Material , data = fungi_data_tr,id = replications:Face, 
                          family = binomial, corstr = "Exchangeable")

fit_marginal_10 <-  glmgee(percentage_mold_growth_tr ~ day + coating*Material , data = fungi_data_tr,id = replications:Face, 
                           family = binomial, corstr = "Exchangeable")
fit_marginal_11 <- glmgee(percentage_mold_growth_tr ~ day * coating * Material , data = fungi_data_tr,id = replications:Face, 
                          family = binomial, corstr = "Exchangeable")
residuals(fit_marginal_11, type="deviance", plot.it=TRUE, pch=16)


QIC_marginal <- c(QIC(fit_marginal_1), QIC(fit_marginal_2), QIC(fit_marginal_3), QIC(fit_marginal_4),
                  QIC(fit_marginal_5), QIC(fit_marginal_6), QIC(fit_marginal_7), QIC(fit_marginal_8),
                  QIC(fit_marginal_9), QIC(fit_marginal_10), QIC(fit_marginal_11))
#
##############--------#####################



residuals(, type="deviance", plot.it=TRUE, pch=16)







##### Interval censored 

fit_surv <- Surv(fungi_data_tr$start, fungi_data_tr$end, type = "interval2")

fit_intervalcensored1 <- gamlss(fit_surv ~ day * coating * Material + random(replications:Face), 
                                data = fungi_data_tr, family = cens(BE, type = "interval"))
wp(fit_intervalcensored1)

drop1(fit_intervalcensored1)

fit_intervalcensored2 <- gamlss(fit_surv ~ (day + coating + Material)^2 + random(replications:Face), 
                                data = fungi_data_tr, family = cens(BE, type = "interval"))
drop1(fit_intervalcensored2)
fit_intervalcensored3 <- gamlss(fit_surv ~ day + coating * Material + random(replications:Face), 
                                data = fungi_data_tr, family = cens(BE, type = "interval"))
drop1(fit_intervalcensored3)                     
fit_intervalcensored4 <- gamlss(fit_surv ~ day * coating + Material + random(replications:Face), 
                                data = fungi_data_tr, family = cens(BE, type = "interval"))
drop1(fit_intervalcensored4) #####this is the significant 

fit_intervalcensored5 <- gamlss(fit_surv ~ day * Material + coating + random(replications:Face), 
                                data = fungi_data_tr, family = cens(BE, type = "interval"))
drop1(fit_intervalcensored5)#####this is significant 
fit_intervalcensored6 <- gamlss(fit_surv ~ day + coating + Material + random(replications:Face), 
                                data = fungi_data_tr, family = cens(BE, type = "interval"))
drop1(fit_intervalcensored6)
fit_intervalcensored7 <- gamlss(fit_surv ~  day + coating + random(replications:Face), 
                                data = fungi_data_tr, family = cens(BE, type = "interval"))
drop1(fit_intervalcensored7)
fit_intervalcensored8 <- gamlss(fit_surv ~  coating + Material + random(replications:Face), 
                                data = fungi_data_tr, family = cens(BE, type = "interval"))
drop1(fit_intervalcensored8)

fit_intervalcensored9 <- gamlss(fit_surv ~  day + Material + random(replications:Face), 
                                data = fungi_data_tr, family = cens(BE, type = "interval"))

fit_intervalcensored10 <- gamlss(fit_surv ~  day + random(replications:Face), 
                                 data = fungi_data_tr, family = cens(BE, type = "interval"))
fit_intervalcensored11 <- gamlss(fit_surv ~  coating + random(replications:Face), 
                                 data = fungi_data_tr, family = cens(BE, type = "interval"))
fit_intervalcensored12 <- gamlss(fit_surv ~  Material + random(replications:Face), 
                                 data = fungi_data_tr, family = cens(BE, type = "interval"))
drop1(fit_intervalcensored2)


### continue here same as was done for fit_mixed then calculate information criteria
### ...
###

AIC_intervalcens <- c(AIC(fit_intervalcensored1),AIC(fit_intervalcensored2),AIC(fit_intervalcensored3),
                      AIC(fit_intervalcensored4),AIC(fit_intervalcensored5),AIC(fit_intervalcensored6),
                      AIC(fit_intervalcensored7),AIC(fit_intervalcensored8),AIC(fit_intervalcensored9),
                      AIC(fit_intervalcensored10),AIC(fit_intervalcensored11),AIC(fit_intervalcensored12))
AIC_intervalcens

BIC_intervalcens <- c(BIC(fit_intervalcensored1),BIC(fit_intervalcensored2),BIC(fit_intervalcensored3),
                      BIC(fit_intervalcensored4),BIC(fit_intervalcensored5),BIC(fit_intervalcensored6),
                      BIC(fit_intervalcensored7),BIC(fit_intervalcensored8),BIC(fit_intervalcensored9),
                      BIC(fit_intervalcensored10),BIC(fit_intervalcensored11),BIC(fit_intervalcensored12))
BIC_intervalcens





