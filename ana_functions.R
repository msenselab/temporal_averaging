# This code reads the raw experimental data from the study
# 'what you see depends on what you hear: temporal rate averaging and crossmodal integration'
# Authors: Lihan Chen, Xiaolin Zhou, Hermann J. Mueller, Zhuanghua Shi
# and generates the figures and statistical results for the paper
# Initial coded in 2016, modified in 2017, and made open-access in 2018
# Zhuanghua Shi, shi@lmu.de
# analyses functions 
# reading data, optimization etc. 

library(tidyverse)
library(broom)
library(multcomp)
library(nlme)
library(ggsignif)
library(quickpsy)
library(cowplot)
library(ez)

# read raw data function
readData <- function(fn) {
  raw <- read.csv(paste0(fn,'.csv'))
  # factorize mean auditory intervals (independent variable)
  raw$mIntv <- factor(raw$mIntv, labels=c('-70', '0','70'))
  raw$sub <- factor(raw$sub) # factorize subjects
  raw
}

# function to plot mean PSE bars
psebar <- function(data){
  data %>% group_by(mIntv) %>% summarise(m = mean(pse), se = sd(pse)/sqrt(nlevels(sub)-1)) %>%
    ggplot(.,aes(x=mIntv, y=m)) + 
    geom_bar(stat='identity', width=0.5, fill = 'grey', color = 'black') + 
    geom_errorbar(aes(ymin=m, ymax=m+se), stat='identity', width=0.3) + 
    coord_cartesian(ylim=c(100,175)) + xlab('Relative mean auditory interval (ms)') + ylab('PSEs (ms)') + 
    theme(legend.position="none")
}

# figure settings
theme_set(theme_bw(base_size = 12))
# control legends etc. theme
legend_pos <-function(l) {
  theme( plot.title = element_text(face="bold",size = 12),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = l,
         legend.title = element_blank(),
         #         legend.text = element_text(size=12),
         legend.key = element_blank() # switch off the rectangle around symbols
  )
}

# optimization function min.rss
# using two parameters model: auditory and Ternus intervals were integrated together
#       the integration weight depending on the discrepancy of auditory and Ternus interval
#       prior is guassian envelope based on interval difference (soa-mi)^2/s^2
#       see Roach et al. (2006), Ernst & Di Luca (2011)
#       or Koerding et al., 2007 - Causual inference p_common decreased as a gaussian function

min.rss2 <- function(data, par) {
  # input averaged data for each condition combination
  # par[1]: variability of ensemble auditory sequence: sig_a
  # par[2]: variability of integration prior (Inf -> full integration)
  # according to the MLE model, the weight of the auditory sequence in full integration mode
  # p = sig_v^2/(sig_v^2+sig_a^2)
  # i.e. p = data$sig^2/(data$sig^2 + par[1]^2)
  data %>% mutate(., p = sig^2/(sig^2+par[1]^2)) %>% # full MLE weight of Audition
    mutate(., w = p* exp(-(soa-mo)^2/par[2]/par[2]), # partial integration based on 
           pdur = (1-w)*soa + w*mo, # linear integration
           mresp = pnorm(a + b*pdur), #estimated mean response based on baseline curve
           rss = (mresp-m)^2) %>% #calcualte the residual
    summarise(.,rss = sum(rss)) # sum the residual and return
}

min.fullfusion <- function(data, par) {
  # full integration model
  # input averaged data for each condition combination
  # par[1]: variability of ensemble auditory sequence: sig_a
  # according to the MLE model, the weight of the auditory sequence in full integration mode
  # p = sig_v^2/(sig_v^2+sig_a^2)
  # i.e. p = data$sig^2/(data$sig^2 + par[1]^2)
  data %>% mutate(., w = sig^2/(sig^2+par[1]^2),
                  pdur = (1-w)*soa + w*mo, # linear integration
                  mresp = pnorm(a + b*pdur), #estimated mean response based on baseline curve
                  rss = (mresp-m)^2) %>% #calcualte the residual
    summarise(.,rss = sum(rss)) # sum the residual and return
}

# prediction based on crossmodal integration
predictResp <- function(data, ispartial = TRUE) {
  # inputs: data, baseline, and model selection: partial integraiton or full integration
  #  print(summary(avg.resp))
  # do optimization to estimate parameters
  # optimization is based on subject-level, not on condition-level
  if (ispartial == TRUE){
    targetfun  <-  min.rss2
    x0 <- c(0.5,100)
    cl <- c(1,0.001)
    cu <- c(300,Inf)
  } else{
    targetfun = min.fullfusion
    x0 <- c(0.5)
    cl <- c(1)
    cu <- c(300)
  }
  data %>% data.frame(.) %>%
    group_by(cond, sub, variance, mIntv) %>% # for each subject, do optimization
    do(tidy(optim(par=x0,targetfun,lower = cl, upper = cu,
                  method="L-BFGS-B", data = .))) %>%  # optimization for parameters
    spread(., parameter, value) -> avg.optim  # store estimated parameters
  
  # estimated mean responses based on baseline psychometric curve and estimated weights
  avg.estimate <- full_join(data, avg.optim, by=c('cond','sub','variance','mIntv')) # combined estimated parameters
  if (ispartial == TRUE) # two parameters
    avg.estimate <- avg.estimate %>% mutate( w=sig^2/(sig^2+parameter1^2)*exp(-(soa-mi)^2/parameter2/parameter2))
  else # one parameters
    avg.estimate <- avg.estimate %>% mutate( w=sig^2/(sig^2+parameter1^2))
  
  avg.estimate <- avg.estimate %>%  
    mutate(pdur = (1-w)*soa + w*mi, mresp = pnorm(a+b*pdur)) # and subjective duration based on weighted integration
  avg.estimate
}

# predicted PSEs
predictPSE <- function(avgResp) {
  # based on estiamted mean response (mresp), calculated estimated PSEs and sigmas
  pse_estimate <- avgResp %>% group_by(exp, sub, variance, mIntv) %>%
    do(tidy (glm(cbind(mresp,1-mresp) ~ soa, family = quasibinomial(probit), data=.))) %>% # pse estimation
    select(., one_of(c( 'term','estimate'))) %>% # presever the parameters only
    spread(.,term, estimate) %>% rename(., b=soa, a = `(Intercept)`) %>% # change name to a and b
    mutate(., ppse = -a/b, psig = abs(1/b))  # predicted pse and sigma

  pse_behavior <- avgResp %>% group_by(exp, sub, variance, mIntv) %>%
    do(tidy (glm(cbind(m,1-m) ~ soa, family = quasibinomial(probit), data=.))) %>% # pse estimation
    select(., one_of(c( 'term','estimate'))) %>% # presever the parameters only
    spread(.,term, estimate) %>% rename(., b=soa, a = `(Intercept)`) %>% # change name to a and b
    mutate(., pse = -a/b, sig = 1/b)  # predicted pse and sigma
  
  left_join(pse_estimate, pse_behavior, by=c('exp','sub','variance','mIntv')) # return
}

# new prediction function with quickpsy. The results were close to the previous methods
predictPSEnew <- function(avgResp) {
  # based on estiamted mean response (mresp), calculated estimated PSEs and sigmas
  pse_est <- avgResp %>% mutate(k0 = round(mresp*ntrials)) %>% 
    quickpsy(d =. , x = soa, k = k0, n = ntrials, 
             grouping = .(exp, sub),
             within = .(mIntv, variance),  
             lapses = TRUE, guess = TRUE, bootstrap = 'none', 
             parini = list(c(1,250),c(1,200),c(0,0.3),c(0,0.3))) 
  
  pse_behav <- avgResp %>% mutate(k0 = round(m*ntrials)) %>% 
    quickpsy(d =. , x = soa, k = k0, n = ntrials, 
             grouping = .(exp, sub),
             within = .(mIntv, variance),  
             lapses = TRUE, guess = TRUE, bootstrap = 'none', 
             parini = list(c(1,250),c(1,200),c(0,0.3),c(0,0.3))) 
  
  pses <- left_join(pse_est$thresholds, pse_behav$thresholds, 
                    by = c('exp','sub','variance','mIntv')) %>%
    rename(ppse = thre.x, pse = thre.y)
  pses
}

# plot predicted responses 
plotPredResponses <- function(data){
  data %>% 
    group_by(mIntv,variance, soa) %>%
    summarise(mr = mean(m), mse = mean(mse),
              mpred = mean(mresp)) %>% 
    ggplot(aes(soa, y = mr, color = mIntv, shape = mIntv, linetype = mIntv)) + geom_point(size = 2) + 
    geom_line(aes(y=mpred)) + 
    geom_errorbar(aes(x = soa, ymin = mr - mse, ymax = mr + mse, color = mIntv), width = 5) + 
    scale_shape(solid = FALSE) + xlab('SOA (ms)') + 
    ylab('Prop. of group motion')+ legend_pos(c(0.2,0.8))
}

# theme black-white
figTheme <- function(){
  theme_bw() +
    theme(strip.background = element_blank(), 
          legend.key = element_blank(), 
          legend.title = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) 
}

# ---- preprocess ----
# experiments: regular, irregular, ari_vs_geo_mean, variance, last_interval manipulations. 
exp_files <- c('exp1','exp2','exp3','exp4','exp5')
exp_names <- c('Regular','Irregular','Geometry','Variance','Last')
raw_data <- map(exp_files,readData) #read all data
names(raw_data) <- exp_names # add names for better readibility
# change levels for Exp3 and Exp5, because those are not -70,0, 70
levels(raw_data$Geometry$mIntv) <-  c('AriM','GeoM','Equal')
levels(raw_data$Last$mIntv) <-  c('Short','Median','Long')
# in variance manipulation experiment, we have additional factor
raw_data$Variance$var <- factor(raw_data$Variance$var, labels=c('CV:0.1','CV:0.3'))
