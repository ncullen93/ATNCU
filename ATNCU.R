


library(ABA)

########################
## COGNITIVE OUTCOMES ##
########################
base.dir <- '~/desktop/studies/atn_cu/'
df <- read.csv(paste0(base.dir,'data/BIOFINDER_ATN_CU_LONG_PROCESSED.csv'))
df$PLASMA_ABETA_bl <- -1*df$PLASMA_ABETA_bl

# model parameters
outcomes <- c('PACC', 'MMSE')
predictors <- list('ATN' = c('PLASMA_ABETA_bl', 'PLASMA_PTAU217_bl', 'PLASMA_NFL_bl'),
                   'A' = 'PLASMA_ABETA_bl',
                   'T' = 'PLASMA_PTAU217_bl',
                   'N' = 'PLASMA_NFL_bl')
covariates <- c('AGE_bl', 'GENDER', 'EDUCAT')
control <- aba_control(coefs.ci=T,
                      metrics.ci=T,
                      aic.norm = T)

# PLASMA ATN for MMSE and PACC
res.lme1 <- aba_lme(data=df,
                    outcomes=outcomes,
                    predictors=predictors,
                    covariates=covariates,
                    id.factor='SUBJECT_ID',
                    time.factor='Years_bl',
                    std.beta=c(F,T))
s.lme1 <- aba_summary(res.lme1, control=control)


# APOE Sensitivity analysis for PACC
res.lme2 <- update(res.lme1,
                   covariates=c(covariates,'APOE'),
                   outcomes=c('PACC'))
s.lme2 <- aba_summary(res.lme2, control=control)

# CSF ATN for PACC
res.lme3 <- update(res.lme1,
                   predictors=list('ATN' = c('CSF_ABETA_bl', 'CSF_PTAU181_bl', 'CSF_NFL_bl'),
                                   'A' = 'CSF_ABETA_bl',
                                   'T' = 'CSF_PTAU181_bl',
                                   'N' = 'CSF_NFL_bl'),
                   outcomes=c('PACC'))
s.lme3 <- summary(res.lme3, control=control)

# Sensitivity analysis - adjustment for diagnostic status
res.lme4 <- update(res.lme1,
                   covariates=c(covariates, 'DX_bl'),
                   interaction.covariates=c('DX_bl'),
                   outcomes=c('PACC'))
s.lme4 <- aba_summary(res.lme4, control=control)


#######################
## CLINICAL OUTCOMES ##
#######################
df <- read.csv(paste0(base.dir,'data/BIOFINDER_ATN_CU_LONG_PROCESSED.csv'))
df <- subset(df, VISIT==0)
df$DX_bl <- factor(df$DX_bl)
df$PLASMA_ABETA_bl <- -1*df$PLASMA_ABETA_bl

# model parameters
outcomes <- c('ConvertedToAlzheimers', 'ConvertedToDementia')
predictors <- list('ATN' = c('PLASMA_ABETA_bl', 'PLASMA_PTAU217_bl', 'PLASMA_NFL_bl'),
                   'A' = 'PLASMA_ABETA_bl',
                   'T' = 'PLASMA_PTAU217_bl',
                   'N' = 'PLASMA_NFL_bl')
covariates <- c('AGE_bl', 'GENDER', 'EDUCAT')
control <- aba_control(coefs.ci=T,
                      metrics.ci=T,
                      aic.norm = T)

# PLASMA ATN for AD dementia and All-cause dementia
res.cox1 <- aba_cox(data=df,
                    outcomes=outcomes,
                    predictors=predictors,
                    covariates=covariates,
                    time.factor='TimeUnderRiskDementia',
                    std.beta=T)

plot.cox1 <- aba_coxplot(res.cox1,
                         outcomes='ConvertedToAlzheimers',
                         predictors=list('Plasma Ab42/Ab40'='A',
                                         'Plasma P-tau217'='T',
                                         'Plasma NfL'='N'),
                         time.vary=c(0,2,4,6),
                         predictor.vary=list('A'='PLASMA_ABETA_STATUS_bl',
                                             'T'='PLASMA_PTAU217_STATUS_bl',
                                             'N'='PLASMA_NFL_STATUS_bl'),
                         predictor.vary.labels=c('Biomarker-Negative',
                                                 'Biomarker-Positive'),
                         xlab='Years from baseline',
                         ylab='Risk of AD dementia (%)',
                         ylim=c(0,22.5))

##########################

s.cox1 <- aba_summary(res.cox1, control=control)

# APOE Sensitivity analysis for AD dementia
res.cox2 <- update(res.cox1,
                   covariates=c(covariates, 'APOE'),
                   outcomes=c('ConvertedToAlzheimers'))
s.cox2 <- aba_summary(res.cox2, control=control)

# Sensitivity analysis - adjustment for diagnostic status
res.cox3 <- update(res.cox1,
                   covariates=c(covariates, 'DX_bl'),
                   outcomes=c('ConvertedToAlzheimers'))
s.cox3 <- aba_summary(res.cox3, control=control)


##########################
## SAVE RESULTS TO FILE ##
##########################
save.all <- T
if (save.all) {
  aba_write(s.lme1, filebase=paste0(base.dir,'results/_LME_PLASMA'))
  aba_write(s.lme2, filebase=paste0(base.dir,'results/_LME_PLASMA_APOE'))
  aba_write(s.lme3, filebase=paste0(base.dir,'results/_LME_CSF'))
  aba_write(s.lme4, filebase=paste0(base.dir,'results/_LME_PLASMA_DX'))
  aba_write(s.cox1, filebase=paste0(base.dir,'results/_COX_PLASMA'))
  aba_write(s.cox2, filebase=paste0(base.dir,'results/_COX_PLASMA_APOE'))
  aba_write(s.cox3, filebase=paste0(base.dir,'results/_COX_PLASMA_DX'))
}


