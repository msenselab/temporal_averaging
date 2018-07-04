
suppressWarnings(suppressMessages(suppressPackageStartupMessages({
  source("ana_functions.R")})))

# ------ Exp1_Regular -------
# due to lapsing and guessing, some participants may not have perfect psychometric curves
# so we estimate these two parameters as well. 
quickpsy(raw_data$Regular, x = soa, k = resp, grouping = .(sub), within = .(mIntv), 
         lapses = TRUE, guess = TRUE, bootstrap = 'none', 
         parini = list(c(1,250),c(1,200),c(0,0.3),c(0,0.3))) -> exp1

# select 1 participant for an example
quickpsy(raw_data$Regular %>% filter(sub == 21), x = soa, k = resp,  within = .(mIntv), 
         lapses = TRUE, guess = TRUE, bootstrap = 'none', 
         parini = list(c(1,250),c(1,200),c(0,0.3),c(0,0.3))) -> exp1sub1
#fig2a <- plotcurves(exp1sub1, threshold = FALSE) 
fig2a <- ggplot() + 
  geom_point(data = exp1sub1$averages, aes(x = soa, y = prob, shape = mIntv), size = 2) + 
  geom_line(data = exp1sub1$curves, aes(x, y , linetype = mIntv)) + 
  legend_pos(c(0.2,0.8)) + xlab ('ISI (ms)') + 
  ylab('Prop. of group motion') + 
  scale_shape_discrete(solid = FALSE)
fig2a

# mean PSE and anovas
exp1$thresholds  %>% group_by(mIntv) %>% summarise(mpse = mean(thre))
exp1$thresholds %>% ezANOVA(data =., dv = thre, wid = sub, within = mIntv)

# posthoc analysis
lme_exp1 = lme(thre ~ mIntv, data = exp1$thresholds, random = ~1|sub)
anova(lme_exp1)
summary(glht(lme_exp1, linfct = mcp(mIntv = 'Tukey')), test = adjusted(type = 'bonferroni'))
# posthoc suggest that significance between -70 and 0, -70 and 70, but not 70 and 0 conditions

# add signficance to the bar figure. 
fig2c <- exp1$thresholds %>% rename(pse = thre) %>% psebar(.)

fig2c <- fig2c + geom_signif(comparisons = list(c('-70','0')), annotation = '*', color = 'black') + 
  geom_signif(comparisons = list(c('-70','70')), annotation = '*', y_position = 160, color = 'black')
fig2c
# mean JND (sigma)
exp1$par %>% filter(parn == 'p2') %>% group_by(mIntv) %>% summarise(msig = mean(par), se = sd(par)/sqrt(nlevels(sub)))
lme_exp1sig = lme(par ~ mIntv, data =exp1$par %>% filter(parn == 'p2'), random = ~1|sub)
anova(lme_exp1sig)
summary(glht(lme_exp1sig, linfct = mcp(mIntv = 'Tukey')), test = adjusted(type = 'bonferroni'))

# ----- Exp2_Irregular ------
# irregular, individual pses
quickpsy(raw_data$Irregular, x = soa, k = resp, grouping = .(sub), within = .(mIntv), 
         lapses = TRUE, guess = TRUE, bootstrap = 'none', 
         parini = list(c(1,250),c(1,200),c(0,0.3),c(0,0.3))) -> exp2

# select 1 participant
quickpsy(raw_data$Irregular %>% filter(sub == 13), x = soa, k = resp,  within = .(mIntv), 
         lapses = TRUE, guess = TRUE, bootstrap = 'none', 
         parini = list(c(1,250),c(1,200),c(0,0.3),c(0,0.3))) -> exp2sub1
#fig2b <- plotcurves(exp2sub1, threshold = FALSE) + legend_pos(c(0.2,0.75)) + xlab ('ISI (ms') + 
#  ylab('Prop. of group motion')
fig2b <- ggplot() + 
  geom_point(data = exp2sub1$averages, aes(x = soa, y = prob, shape = mIntv), size = 2) + 
  geom_line(data = exp2sub1$curves, aes(x, y , linetype = mIntv)) + 
  legend_pos(c(0.2,0.8)) + xlab ('ISI (ms)') + 
  ylab('Prop. of group motion') + 
  scale_shape_discrete(solid = FALSE)
fig2b

# plot means and ANOVA analysis
exp2$thresholds  %>% group_by(mIntv) %>% summarise(mpse = mean(thre))
exp2$thresholds %>% ezANOVA(data =., dv = thre, wid = sub, within = mIntv)

# posthoc analysis
lme_exp2 = lme(thre ~ mIntv, data = exp2$thresholds, random = ~1|sub)
anova(lme_exp2)
summary(glht(lme_exp2, linfct = mcp(mIntv = 'Tukey')), test = adjusted(type = 'bonferroni'))
# same story: signicances between -70 and 0, -70 and 70, but not 70 and 0. 
# create average bar plots
fig2d <- exp2$thresholds %>% rename(pse = thre) %>% psebar(.)
fig2d <- fig2d + geom_signif(comparisons = list(c('-70','0')), annotation = '*', color = 'black') + 
  geom_signif(comparisons = list(c('-70','70')), annotation = '*', y_position = 170, color = 'black')
fig2d
# mean JND (sigma)
exp2$par %>% filter(parn == 'p2') %>% group_by(mIntv) %>% summarise(msig = mean(par), se = sd(par)/sqrt(nlevels(sub)))
lme_exp2sig = lme(par ~ mIntv, data =exp2$par %>% filter(parn == 'p2'), random = ~1|sub)
anova(lme_exp2sig)
summary(glht(lme_exp2sig, linfct = mcp(mIntv = 'Tukey')), test = adjusted(type = 'bonferroni'))

# combine subplots together into Figure 2
fig2 <- plot_grid(fig2a,fig2b, fig2c, fig2d, nrow = 2, labels = c('A','B','C','D'))
fig2
ggsave('figure2.pdf',fig2, width=8, height = 6)

# check laspe rate 
exp2$par %>% filter(parn == 'p4') %>% ezANOVA(data=., dv = par, wid = sub, within = mIntv)
# it is significant, which suggests short interval has a general bias toward the element motion.  

# save mean thresholds
if (FALSE) { # save once
  write.csv(exp1$thresholds,file = 'exp_reg_pses.csv')
  write.csv(exp2$thresholds,file = 'exp_irreg_pses.csv')
  # save adjusted mean, sigma, guess, lapsed parameters (p1, p2, p3, p4)
  write.csv(exp1$par, file = 'exp_reg_par.csv')
  write.csv(exp2$par, file = 'exp_irreg_par.csv')
}

#---- Exp3_variance ----
# variance manipulation experiment
quickpsy(raw_data$Variance, x = soa, k = resp, grouping = .(sub), within = .(var, mIntv), 
         lapses = TRUE, guess = TRUE, bootstrap = 'none', 
         parini = list(c(1,250),c(1,200),c(0,0.3),c(0,0.3))) -> exp3

# plot means and ANOVA analysis
exp3$thresholds %>%  group_by(var, mIntv) %>% summarise(mpse = mean(thre), se = sd(thre)/sqrt(nlevels(sub)-1)) -> mexp3
exp3$thresholds  %>% ezANOVA(data =., dv = thre, wid = sub, within = .(var, mIntv))

pd <- position_dodge(width = 0.1)
fig3 <- ggplot(mexp3, aes(mIntv, mpse, ymin = mpse - se, ymax = mpse + se, 
                          linetype = var, shape = var, group = var)) + 
  geom_line(position = pd) + geom_errorbar(width = 0.2, position = pd) + 
  geom_point(size = 2, position = pd) + 
  legend_pos(c(0.8,0.85)) + 
  scale_shape_discrete(solid = FALSE) +
  xlab('Relative mean auditory intervals (ms)') + ylab('PSEs (ms)')
fig3
ggsave(filename = 'figure3.pdf', fig3, width = 6, height = 4)
# posthoc analysis
lme_exp3 = lme(thre ~ mIntv * var, data = exp3$thresholds, random = ~1|sub)
anova(lme_exp3)
summary(glht(lme_exp3, linfct = mcp(mIntv = 'Tukey')), test = adjusted(type = 'bonferroni'))


# ---- Exp4_Geometric_manipulation ----
# arithmetic one
quickpsy(raw_data$Geometry, x = soa, k = resp, grouping = .(sub), within = .(mIntv), 
         lapses = TRUE, guess = TRUE, bootstrap = 'none', 
         parini = list(c(1,250),c(1,200),c(0,0.3),c(0,0.3))) -> exp4
# select one participant
quickpsy(raw_data$Geometry %>% filter(sub == 6), x = soa, k = resp, within = .(mIntv), 
         lapses = TRUE, guess = TRUE, bootstrap = 'none', 
         parini = list(c(1,250),c(1,200),c(0,0.3),c(0,0.3))) -> exp4sub1
#fig4a <- plotcurves(exp4sub1, threshold = FALSE) + legend_pos(c(0.2,0.75)) + xlab ('ISI (ms') + 
#  ylab('Prop. of group motion')
fig4a <- ggplot() + 
  geom_point(data = exp4sub1$averages, aes(x = soa, y = prob, shape = mIntv), size = 2) + 
  geom_line(data = exp4sub1$curves, aes(x, y , linetype = mIntv)) + 
  legend_pos(c(0.2,0.8)) + xlab ('ISI (ms)') + 
  ylab('Prop. of group motion') + 
  scale_shape_discrete(solid = FALSE)
# plot means and ANOVA analysis
exp4$thresholds  %>% group_by( mIntv) %>% summarise(mpse = mean(thre), se = sd(thre)/sqrt(nlevels(sub)-1)) -> mexp4
exp4$thresholds  %>% ezANOVA(data =., dv = thre, wid = sub, within = .( mIntv))

# plot bar figure with significances
fig4b <- exp4$thresholds %>% rename(pse = thre) %>% psebar(.)
fig4b <- fig4b + geom_signif(comparisons = list(c('GeoM','AriM')), annotation = '*', y_position = 165, color = 'black') + 
  geom_signif(comparisons = list(c('GeoM','Equal')), annotation = '*', y_position = 160, color = 'black')
fig4b

# combine two subplots together into Figure 4
fig4 <- plot_grid(fig4a,fig4b, nrow = 1, labels = c('A','B'))
fig4
ggsave(filename = 'figure4.pdf', fig4, width = 8, height = 4)
# posthoc analysis
lme_exp4 = lme(thre ~ mIntv, data = exp4$thresholds, random = ~1|sub)
anova(lme_exp4)
summary(glht(lme_exp4, linfct = mcp(mIntv = 'Tukey')), test = adjusted(type = 'bonferroni'))

# a control experiment

dat6 = read.csv('exp6.csv', header = FALSE)
names(dat6) <-  c('motion','Interval','pos','soa','resp','rt','dump','subno')
dat6$soa = dat6$soa*30 + 20
dat6$Interval = factor(dat6$Interval, labels = c("Short","Long"))

quickpsy(dat6, x = soa, k = resp, grouping = .(subno), within = .(Interval), 
         lapses = TRUE, guess = TRUE, bootstrap = 'none', 
         parini = list(c(1,250),c(1,200),c(0,0.3),c(0,0.3))) -> exp6

# select one participant
quickpsy(dat6 %>% filter(subno == 3), x = soa, k = resp, within = .(Interval), 
         lapses = TRUE, guess = TRUE, bootstrap = 'none', 
         parini = list(c(1,250),c(1,200),c(0,0.3),c(0,0.3))) -> exp6sub1

fig5 <- ggplot() + 
  geom_point(data = exp6sub1$averages, aes(x = soa, y = prob, shape = Interval), size = 2) + 
  geom_line(data = exp6sub1$curves, aes(x, y , linetype = Interval)) + 
  legend_pos(c(0.2,0.8)) + xlab ('ISI (ms)') + 
  ylab('Prop. of group motion') + 
  scale_shape_discrete(solid = FALSE)
fig5
ggsave(filename = 'figure5.pdf', fig5, width = 4, height = 4)

# Bayesian Modeling  
# read mean responses from formal experiments, together with the baseline conditions
meanResps <- readRDS('meanResps.rds')

# ---- integration_model ----
runOptim = FALSE
if (runOptim){
  resp.predicted <- predictResp(meanResps) # partial integration
  resp.fullfusion <- predictResp(meanResps, FALSE)# full integration
  saveRDS(resp.predicted,'partial_integration_new.rds')
  saveRDS(resp.fullfusion,'full_integration_new.rds')
}else{
  resp.predicted <- readRDS('partial_integration_new.rds')
  resp.fullfusion <- readRDS('full_integration_new.rds')
}

# sse, r^2, BIC
stat1 <- resp.predicted %>% group_by(exp) %>% 
  summarise(sse = sum((m-mresp)^2), 
            sst = sum((m- mean(resp.fullfusion$m))^2), 
            n = n()) %>% 
  mutate(bic_partial = n*log(sse/n) + 2*log(n), r2_partial = 1-sse/sst)

stat2 <- resp.fullfusion %>% group_by(exp) %>% 
  summarise(sse = sum((m-mresp)^2), 
            sst = sum((m- mean(resp.fullfusion$m))^2), 
            n = n()) %>% 
  mutate(bic_full = n*log(sse/n) + 2*log(n), r2_full = 1 - sse/sst)

stat_tab <- left_join(stat1, stat2, by = c('exp')) %>% 
  mutate(delta = bic_full - bic_partial) %>%
  dplyr::select(exp, bic_partial, r2_partial, bic_full, r2_full, delta) 
colnames(stat_tab) <- c("Exp.", "BIC P.I.", "R2 P.I.", 
                        "BIC F.I.", "R2 F.I.", "BIC Diff.")

knitr::kable(stat_tab, digits=2)

# plot predicted responses
fig6 <- resp.predicted %>%
  #resp.fullfusion%>%
  group_by(cond, mIntv,soa ) %>%
  summarise(mr = mean(m), mse = mean(mse),
            mpred = mean(mresp)) %>% 
  ggplot(aes(soa, y = mr,  shape = mIntv)) + 
  geom_point(size = 2) + 
  geom_line(aes(y=mpred, linetype = mIntv)) + 
  geom_errorbar(aes(x = soa, ymin = mr - mse, ymax = mr + mse), width = 5) + 
  scale_shape(solid = FALSE) + xlab('ISI (ms)') + 
  ylab('Prop. of group motion')+  facet_wrap(~cond, nrow = 2) + figTheme()

fig6
ggsave(filename = 'figure6.pdf', fig6, width = 8, height = 6)

# predict weights as function SOA
fig7 = resp.predicted %>% group_by(exp, mIntv, soa) %>% summarise(mw = mean(w)) %>% 
  ggplot(aes(x = soa, y = mw, shape = mIntv)) + geom_line(aes( linetype = mIntv)) + 
  geom_point(size = 2) +facet_wrap(~exp) + xlab('ISI (ms)') + ylab('Means auditory weights (wPam)') + 
  scale_shape(solid = FALSE) + figTheme()
fig7
ggsave(filename = 'figure7.pdf', fig7, width = 8, height = 4)

pse.predicted <- predictPSEnew(resp.predicted)
#pse.predicted <- predictPSE(resp.fullfusion)
# the following is to obtain slope 0.978 and R^2= 0.983
pse.lm <- pse.predicted %>% 
  lm( ppse ~ pse -1, data = .)
summary(pse.lm)

fig_pses <- pse.predicted %>% 
  group_by(exp, mIntv) %>% #summarise(pse = mean(pse), ppse = mean(ppse)) %>%
  ggplot(., aes(x=pse, y = abs(ppse), shape = exp)) + geom_point() +geom_abline(slope = 1, linetype = 'dashed') + xlab('Observed PSEs (ms)') + ylab('Predicted PSEs (ms)')+ 
  figTheme() 
fig_pses
ggsave(filename = 'figure8.pdf', fig_pses, width = 5, height = 3)
