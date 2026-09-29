library(lme4)

setwd("~/Desktop/MultilevelModeling/Project/scripts")

# load libraries
library(lme4)
library(lmerTest)  # Wald p-values (used for lmer; glmer uses z)
library(ordinal)   # clmm for cumulative link mixed models
library(ggplot2)
library(performance)
#library(DHARMa)
library(dplyr)
library(sjPlot)
#library(arm)       # binnedplot
#library(see)        # companion visualisation package for performance

########################################################################
########################################################################
########################################################################
# DATA SET UP
########################################################################
########################################################################
########################################################################

# load data
sch <- readRDS("../data/CY08MSP_SCH_QQQ.rdata")
stu <- readRDS("../data/CY08MSP_STU_QQQ.rdata")

# merge into one dataframe with relevant variables
sch_rel <- sch[, c("CNTSCHID", "W_SCHGRNRABWT", "SC211Q03JA")]

stu_rel <- stu[, c("MATHEFF", "CNTSCHID","W_FSTUWT",
                   # MATHEFF vars
                   "ST290Q01WA", "ST290Q02WA", "ST290Q03WA", "ST290Q04WA",
                   "ST290Q05WA", "ST290Q06WA", "ST290Q07WA", "ST290Q08WA", "ST290Q09WA",
                   #TEACHSUP vars
                   "ST270Q01JA", "ST270Q02JA", "ST270Q03JA", "ST270Q04JA",
                   # Predictors (ST004D01T = GENDER)
                   "TEACHSUP", "ESCS","ST004D01T",
                   # Predictor LANGMATCH
                   "LANGN", "LANGTEST_QQQ",
                   # Predictor ANXFAIL
                   "ST292Q06JA"
                   )]

pisa <- merge(
  stu_rel,
  sch_rel,
  by = "CNTSCHID",
  all.x = TRUE
)

########################################################################
# Calculated values and factorization
########################################################################
# calculate LANGMATCH

pisa$LANGN[pisa$LANGN == 999] <- NA # re code missing language as NA
#dummy code: if language spoken at home matches language of test
pisa$LANGMATCH <- factor(ifelse(pisa$LANGN == pisa$LANGTEST_QQQ, 1, 0))

#Anxious to Fail
pisa$ANXFAIL <- factor(pisa$ST292Q06JA, levels = 1:4,
                      labels = c("Strongly Agree", "Agree", "Disagree", "Strongly Disagree"),
                      ordered = TRUE)

#Gender: 1 = Female,   2 = Male
pisa$GENDER <- factor(
  pisa$ST004D01T,
  levels = c(1, 2),
  labels = c("Female", "Male")
)

pisa$CNTSCHID <- as.factor(pisa$CNTSCHID)

#calculate average ESCS for each school
schools <- unique(pisa$CNTSCHID)
school_escs <- numeric(length(schools))

for (i in seq_along(schools)) {

  schooldat <- pisa[pisa$CNTSCHID == schools[i], ]
  schooldat <- schooldat[!is.na(schooldat$ESCS), ]
  school_escs[i] <- sum(schooldat$ESCS * schooldat$W_FSTUWT) / sum(schooldat$W_FSTUWT)
}

# # rank by avg ESCS
#
schools_data <- data.frame(
  CNTSCHID = schools,
  SCH_ESCS = school_escs
)
# sorts descending: smallest to largest
# low ESCS --> more disadvantaged
schools_data <- schools_data[order(schools_data$SCH_ESCS), ]

# divide into quarters:
schools_data$quarter <- ceiling(
  4 * seq_len(nrow(schools_data)) / nrow(schools_data)
)
# 1 quarter = lowest ESCS = se_disadvantaged
# top quarter = se_advantaged
schools_data$SCH_ESCS_CAT <- ifelse(
  schools_data$quarter == 1,
  1,
  ifelse(
    schools_data$quarter == 4,
    3,
    2
  )
)

pisa <- merge(
  pisa,
  schools_data[, c("CNTSCHID", "SCH_ESCS", "SCH_ESCS_CAT")],
  by = "CNTSCHID",
  all.x = TRUE
)

pisa$SCH_ESCS_CAT <- factor(
                        pisa$SCH_ESCS_CAT,
                        levels = c(1, 2, 3),
                        labels = c("SES Disadvantaged", "Middle", "SES Advantaged")
                            )

########################################################################
# Null values 
########################################################################
npisa <- nrow(pisa) 
npisa # 6116 student observations
length(unique(pisa$CNTSCHID)) # 257 schools 

nc_var = "ANXFAIL"

vars <- c("MATHEFF", nc_var, "ESCS", "LANGMATCH", "TEACHSUP","SC211Q03JA")
vars_c <- c("MATHEFF","ESCS", "LANGMATCH", "TEACHSUP","SC211Q03JA")
vars_nc <- c(nc_var, "ESCS", "LANGMATCH", "TEACHSUP","SC211Q03JA")

# check how many are null
sapply(pisa[c(vars)], function(x) sum(is.na(x)))

# total null 
null_c <- sum(!complete.cases(pisa[c(vars_c)]))
null_nc <- sum(!complete.cases(pisa[c(vars_nc)]))

null_c/npisa #0.3736102
null_nc/npisa#0.4758012

# drop if null 
pisa_c <- pisa[complete.cases(pisa[, vars_c]), ]
pisa_nc <- pisa[complete.cases(pisa[, vars_nc]), ]
pisa_nn <- pisa[complete.cases(pisa[, vars]), ]

pisa_c_null <- pisa[!complete.cases(pisa[, vars_c]), ]
pisa_escs <- pisa_c_null[complete.cases(pisa[, "ESCS"]), ]

pisa_main_cont <- pisa_c

pisa_extend_ncont <- pisa_nc

nrow(pisa_main_cont)
length(unique(pisa_main_cont$CNTSCHID))
nrow(pisa_extend_ncont)
length(unique(pisa_extend_ncont$CNTSCHID))

####################################
# CENTERING, SCALING
####################################
# grand mean centering and scale
#pisa_main_cont
pisa_main_cont$TEACHSUP_C = pisa_main_cont$TEACHSUP - mean(pisa_main_cont$TEACHSUP, na.rm = TRUE)
pisa_main_cont$TEACHSUP_z = pisa_main_cont$TEACHSUP_C / sd(pisa_main_cont$TEACHSUP, na.rm = TRUE)

pisa_main_cont$ESCS_C = pisa_main_cont$ESCS - mean(pisa_main_cont$ESCS, na.rm = TRUE)
pisa_main_cont$ESCS_z = pisa_main_cont$ESCS_C / sd(pisa_main_cont$ESCS, na.rm = TRUE)

pisa_main_cont$SCH_PERCDIS_c = pisa_main_cont$SC211Q03JA - mean(pisa_main_cont$SC211Q03JA, na.rm = TRUE)
pisa_main_cont$SCH_PERCDIS_z = pisa_main_cont$SCH_PERCDIS_c / sd(pisa_main_cont$SC211Q03JA, na.rm = TRUE)

#pisa_extend_ncont
pisa_extend_ncont$TEACHSUP_C = pisa_extend_ncont$TEACHSUP - mean(pisa_extend_ncont$TEACHSUP, na.rm = TRUE)
pisa_extend_ncont$TEACHSUP_z = pisa_extend_ncont$TEACHSUP_C / sd(pisa_extend_ncont$TEACHSUP, na.rm = TRUE)

pisa_extend_ncont$ESCS_C = pisa_extend_ncont$ESCS - mean(pisa_extend_ncont$ESCS, na.rm = TRUE)
pisa_extend_ncont$ESCS_z = pisa_extend_ncont$ESCS_C / sd(pisa_extend_ncont$ESCS, na.rm = TRUE)

pisa_extend_ncont$SCH_PERCDIS_c = pisa_extend_ncont$SC211Q03JA - mean(pisa_extend_ncont$SC211Q03JA, na.rm = TRUE)
pisa_extend_ncont$SCH_PERCDIS_z = pisa_extend_ncont$SCH_PERCDIS_c / sd(pisa_extend_ncont$SC211Q03JA, na.rm = TRUE)

####################################
# WEIGHTS
####################################
# Re-scale the student weights to sum to sample size n (not population size)
# WHY IS THIS NECESSARY FOR LMER?
n <- nrow(pisa_main_cont) 
pisa_main_cont$W_HOUSEWT <- n * pisa_main_cont$W_FSTUWT / sum(pisa_main_cont$W_FSTUWT)

########################################################################
# Demographics
########################################################################

########################################################################
# MATHEFF
########################################################################
mean(pisa_main_cont$MATHEFF)
sd(pisa_main_cont$MATHEFF)

# Check MATHEFF items
MATHEFF_items <- c("ST290Q01WA", "ST290Q02WA", "ST290Q03WA", "ST290Q04WA",
                   "ST290Q05WA", "ST290Q06WA", "ST290Q07WA", "ST290Q08WA", "ST290Q09WA")
# variable to track total NA
pisa$MATHEFF_ITEM_NA_COUNT <- rowSums(is.na(pisa[MATHEFF_items]))
table(pisa$MATHEFF_ITEM_NA_COUNT)

# % skipped this item AND also missing 4-6 items total
pct_skipped_4to6 <- sapply(MATHEFF_items, function(x) {
  sum(is.na(pisa[[x]]) & pisa$MATHEFF_ITEM_NA_COUNT >= 4 & pisa$MATHEFF_ITEM_NA_COUNT <= 6) / nrow(pisa) * 100
})

# % skipped this item AND also missing 7+ items total
pct_skipped_7plus <- sapply(MATHEFF_items, function(x) {
  sum(is.na(pisa[[x]]) & pisa$MATHEFF_ITEM_NA_COUNT >= 7) / nrow(pisa) * 100
})

print(pct_skipped_4to6)
print(pct_skipped_7plus)

#check how this aligns with MATHEFF 
table(pisa$MATHEFF_ITEM_NA_COUNT, is.na(pisa$MATHEFF))
#MATHEFF null if 7+ skipped

colSums(is.na(pisa[MATHEFF_items]))
sum(is.na(pisa$MATHEFF))


########################################################################
# ANXFAIL
########################################################################
table(pisa_extend_ncont$ANXFAIL)
round(prop.table(table(pisa_extend_ncont$ANXFAIL)),5)

########################################################################
# Gender
########################################################################
table(pisa_main_cont$GENDER)
round(prop.table(table(pisa_main_cont$GENDER)),5)

table(pisa_extend_ncont$GENDER)
round(prop.table(table(pisa_extend_ncont$GENDER)),5)

########################################################################
# LANGMATCH 
########################################################################
table(pisa_main_cont$LANGMATCH)
round(prop.table(table(pisa_main_cont$LANGMATCH)),5)

table(pisa_extend_ncont$LANGMATCH)
round(prop.table(table(pisa_extend_ncont$LANGMATCH)),5)

########################################################################
# TEACHSUP
########################################################################

weighted_sd <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]
  w <- w[ok]
  
  sqrt(
    sum(w * (x - weighted.mean(x, w))^2) /
      (sum(w) - sum(w^2) / sum(w))
  )
}

TEACHSUP_items  <- c("ST270Q01JA", "ST270Q02JA", "ST270Q03JA", "ST270Q04JA")

# variable to track total NA
colSums(is.na(pisa[TEACHSUP_items]))

teachsup_dem <- data.frame(
  Statistic = c("Mean", "SD", "Min", "Max"),
  
  Main = c(
    weighted.mean(pisa_main_cont$TEACHSUP,
                  pisa_main_cont$W_HOUSEWT,
                  na.rm = TRUE),
    weighted_sd(pisa_main_cont$TEACHSUP,
                pisa_main_cont$W_HOUSEWT),
    min(pisa_main_cont$TEACHSUP, na.rm = TRUE),
    max(pisa_main_cont$TEACHSUP, na.rm = TRUE)
  ),
  
  Extended = c(
    mean(pisa_extend_ncont$TEACHSUP, na.rm = TRUE),
    sd(pisa_extend_ncont$TEACHSUP, na.rm = TRUE),
    min(pisa_extend_ncont$TEACHSUP, na.rm = TRUE),
    max(pisa_extend_ncont$TEACHSUP, na.rm = TRUE)
  )
)

teachsup_dem

########################################################################
# ESCS
########################################################################
escs_dem <- data.frame(
  Statistic = c("Mean", "SD", "Min", "Max"),
  
  Main = c(
    weighted.mean(pisa_main_cont$ESCS,
                  pisa_main_cont$W_HOUSEWT,
                  na.rm = TRUE),
    weighted_sd(pisa_main_cont$ESCS,
                pisa_main_cont$W_HOUSEWT),
    min(pisa_main_cont$ESCS, na.rm = TRUE),
    max(pisa_main_cont$ESCS, na.rm = TRUE)
  ),
  
  Extended = c(
    mean(pisa_extend_ncont$ESCS, na.rm = TRUE),
    sd(pisa_extend_ncont$ESCS, na.rm = TRUE),
    min(pisa_extend_ncont$ESCS, na.rm = TRUE),
    max(pisa_extend_ncont$ESCS, na.rm = TRUE)
  )
)

escs_dem

########################################################################
# SC211Q03JA (School Percent Disadvantaged)
########################################################################

pisa_main_schools <- pisa_main_cont[!duplicated(pisa_main_cont$CNTSCHID), ]
pisa_extend_schools <- pisa_extend_ncont[!duplicated(pisa_extend_ncont$CNTSCHID), ]
sch_perc_dis_dem <- data.frame(
  Statistic = c("Mean", "SD", "Min", "Max"),
  Main = c(
    mean(pisa_main_schools$SC211Q03JA, na.rm = TRUE),
    sd(pisa_main_schools$SC211Q03JA, na.rm = TRUE),
    min(pisa_main_schools$SC211Q03JA, na.rm = TRUE),
    max(pisa_main_schools$SC211Q03JA, na.rm = TRUE)
  ),
  Extended = c(
    mean(pisa_extend_schools$SC211Q03JA, na.rm = TRUE),
    sd(pisa_extend_schools$SC211Q03JA, na.rm = TRUE),
    min(pisa_extend_schools$SC211Q03JA, na.rm = TRUE),
    max(pisa_extend_schools$SC211Q03JA, na.rm = TRUE)
  )
)

sch_perc_dis_dem

########################################################################
# SCHOOL_ESCS for NULL
########################################################################
# MATHEFF subset: unique schools
sch_matheff <- pisa_escs %>%
  distinct(CNTSCHID, .keep_all = TRUE)

table(sch_matheff$SCH_ESCS_CAT)
round(prop.table(table(sch_matheff$SCH_ESCS_CAT)), 5)

schools_data <- merge(
  schools_data,
  unique(pisa_escs[, c("CNTSCHID", "MATHEFF")]),
  by = "CNTSCHID",
  all.x = TRUE
)

boxplot(
  MATHEFF ~ SCH_ESCS_CAT,
  data = schools_data,
  main = "Distribution of MATHEFF by School Type",
  xlab = "SCH_ESCS_CAT",
  ylab = "MATHEFF"
)

########################################################################
# Function for p-value corrections
########################################################################

LLtest <- function(LL,nrandom){
  #LL is a likelihood ratio value from the anova function
  #nrandom is the number of uncorrelated random effects tested simultaneously
  # returns a p value
  if(nrandom==1){
    # test k uncorrelated random effects vs. k+1 uncorrelated random effects
    1-(.5*pchisq(LL,df=0)+.5*pchisq(LL,df=1))
  }else if(nrandom==2){
    # test k uncorrelated random effects vs. k+2 uncorrelated random effects
    1-(.25*pchisq(LL,df=0)+.5*pchisq(LL,df=1)+.25*pchisq(LL,df=2))
  }else if(nrandom==3){
    # test k uncorrelated random effects vs. k+3 uncorrelated random effects
    1-(.125*pchisq(LL,df=0)+.375*pchisq(LL,df=1)+.375*pchisq(LL,df=2)+.125*pchisq(LL,df=3))
  }else if (nrandom==4){
    # test k uncorrelated random effects vs. k+4 uncorrelated random effects
    1-(.0625*pchisq(LL,df=0)+.25*pchisq(LL,df=1)+.375*pchisq(LL,df=2)+
         .25*pchisq(LL,df=3)+.0625*pchisq(LL,df=4))
  }else{NA}
}

########################################################################
########################################################################
########################################################################
# CONTINUOUS MODEL
########################################################################
########################################################################
########################################################################

####################################
# ICC
####################################

# random intercept model
m0 <- lmer(MATHEFF ~ 1 + (1|CNTSCHID), 
             data=pisa_main_cont, 
             weights = W_HOUSEWT,
             REML=TRUE)

summary(m0)
# Extract variances
covrand <- as.data.frame(VarCorr(m0))
# Between school variance
var_school <- covrand$vcov[covrand$grp == "CNTSCHID"]
# Within school (residual) variance
var_student <- covrand$vcov[covrand$grp == "Residual"]

ICC <- var_school / (var_school + var_student)
ICC # 0.1130553



####################################
# Test random effects
####################################
# saturated model
m1 <- lmer(
  MATHEFF ~ 
    ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z + 
    LANGMATCH:SCH_PERCDIS_z + 
    (1 + ESCS_z + TEACHSUP_z | CNTSCHID), # random slope and intercept
  data = pisa_main_cont, 
  weights = W_HOUSEWT,
  REML=TRUE)

AIC(m1)
BIC(m1)

# uncorrelated
m1a <- lmer(
  MATHEFF ~ 
    ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z + 
    LANGMATCH:SCH_PERCDIS_z + 
    (-1 + TEACHSUP_z| CNTSCHID) + 
    (-1 + ESCS_z| CNTSCHID) + 
    (1 | CNTSCHID), 
  data = pisa_main_cont, 
  weights = W_HOUSEWT,
  REML=TRUE)

anova(m1a,m1, refit=F)
#BIC lower for model m1a: AIC and BIC disagree, go with BIC

# model with no random slope for ESCS_z
m1b <- lmer(
  MATHEFF ~ 
    ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z + 
    LANGMATCH:SCH_PERCDIS_z + 
    (-1 + TEACHSUP_z| CNTSCHID) + 
    (1 | CNTSCHID), 
  data = pisa_main_cont, 
  weights = W_HOUSEWT,
  REML=TRUE)

anova(m1b,m1a, refit=F) 
test_b <- anova(m1b, m1a, refit=F)
test_b
test_b$Df[2]
test_b$Chisq[2]
LLtest(test_b$Chisq[2],1) # not significant at alpha = .05

# model with no random slopes
m1c <- lmer(
  MATHEFF ~ 
    ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z + 
    LANGMATCH:SCH_PERCDIS_z + 
    (1 | CNTSCHID), 
  data = pisa_main_cont, 
  weights = W_HOUSEWT,
  REML=TRUE)

anova(m1c,m1b, refit=F) 
test_c <- anova(m1c, m1a, refit=F)
test_c
test_c$Df[2]
test_c$Chisq[2]
LLtest(test_c$Chisq[2],2) # significant at alpha = .05

#final random effects structure
m2 <- lmer(
  MATHEFF ~ 
    ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z + 
    LANGMATCH:SCH_PERCDIS_z + 
    (-1 + TEACHSUP_z| CNTSCHID) + 
    (1 | CNTSCHID), 
  data = pisa_main_cont, 
  weights = W_HOUSEWT,
  REML=FALSE)

summary(m2)

####################################
# Test interaction effects
####################################
# no interaction effects
m2a <- lmer(
  MATHEFF ~ 
    ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z + 
    (-1 + TEACHSUP_z| CNTSCHID) + 
    (1 | CNTSCHID), 
  data = pisa_main_cont, 
  weights = W_HOUSEWT,
  REML=FALSE)

anova(m2a, m2) # significant at alpha = .05
anova(m2a, m2)$Chisq[2]
anova(m2a, m2)$Df[2]
# final interaction effects structure is m2.

####################################
# main effects
####################################

# remove TEACHSUP_z 
m3a <- lmer(
  MATHEFF ~ 
     GENDER + ESCS_z + SCH_PERCDIS_z + LANGMATCH + 
    LANGMATCH:SCH_PERCDIS_z + 
    (-1 + TEACHSUP_z| CNTSCHID) + 
    (1 | CNTSCHID),
  data = pisa_main_cont, 
  weights = W_HOUSEWT,
  REML=FALSE)

anova(m3a, m2) # significant, keep m2
anova(m3a, m2)$Chisq[2]
anova(m3a, m2)$Df[2]

# remove ESCS_z 
m3b <- lmer(
  MATHEFF ~ 
    TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z + 
    LANGMATCH:SCH_PERCDIS_z + 
    (-1 + TEACHSUP_z| CNTSCHID) + 
    (1 | CNTSCHID), 
  data = pisa_main_cont, 
  weights = W_HOUSEWT,
  REML=FALSE)

anova(m3b, m2) # significant
anova(m3b, m2)$Chisq[2]
anova(m3b, m2)$Df[2]

# final model
mFinal <- lmer(
  MATHEFF ~ 
    ESCS_z + TEACHSUP_z + GENDER + SCH_PERCDIS_z + LANGMATCH +  
    LANGMATCH:SCH_PERCDIS_z + 
    (-1 + TEACHSUP_z| CNTSCHID) + 
    (1 | CNTSCHID), 
  data = pisa_main_cont, 
  weights = W_HOUSEWT,
  REML=TRUE)

#############################################
# Main Findings
#############################################
summary(mFinal)
VarCorr(mFinal)

sjPlot::tab_model(mFinal)

#############################################
# Conduct a residual analysis for ϵ and u.
#############################################

# check for multicollinearity
check_collinearity(mFinal)

# Independence: The error terms should be independent.

# Residuals vs. Fitted Values Plot (Tests Equal Variance)
# Student-level residuals
eps <- resid(mFinal)
# Fitted values
yhat <- fitted(mFinal)

pdf("../figures/resid_vs_fitted.pdf", width = 6, height = 4)
#png("../figures/resid_vs_fitted.png", width = 6, height = 4, units = "in", res = 300)
plot(
  yhat, eps,
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "Residuals vs Fitted Values"
  )
abline(h = 0, lty = 2)
dev.off()

png("../figures/resid_vs_fitted_fan.png", width = 6, height = 4, units = "in", res = 300)
ggplot(data.frame(
  fitted = fitted(mFinal),
  residuals = residuals(mFinal)
), aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Fitted values",
    y = "Residuals",
    title = "Residuals vs Fitted Values"
  ) +
  theme_minimal()
dev.off()
# Linearity: The relationship between the predictors 
# and the outcome must be linear.
# no curve or funnel, but there are horizontal/diagonal bands

# Homoscedasticity (Equal Variance): The residual variance 
# should remain constant across all predicted values.
# no clear fanning on either side - 
# variance appears to be relatively equal


# Normality of Residuals: The individual-level errors (residuals) 
# should be normally distributed.

# Test normality qq plot
pdf("../figures/qq_residuals.pdf", width = 6, height = 4)
#png("../figures/qq_residuals.png", width = 6, height = 4, units = "in", res = 300)
qqnorm(resid(mFinal), main = "Normal Q-Q Plot: Residuals")
qqline(resid(mFinal), col = "black")
dev.off()
# heavy-tailed residual distribution QQ-plot
# indicates that the assumption of normality 
# is not fulfilled.
# Large sample size: CLM

# Normality of Random Effects: The random intercepts 
# for schools should follow a normal distribution 
# around the population average.
eb_estimates <- ranef(mFinal)
u <- ranef(mFinal)$CNTSCHID[[1]]

pdf("../figures/randeff_QQ.pdf", width = 6, height = 4)
#png("../figures/randeff_QQ.png", width = 6, height = 4, units = "in", res = 300)
qqnorm(u,
       main = "Normal Q-Q Plot: School Random Effects")
qqline(u)
dev.off()

# Extract Empirical Bayes estimates (random effects)
names(eb_estimates)

# Histogram to check for School-level outliers
ggplot(mapping = aes(x = u)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Distribution of School Random Effects", 
       x = "School Intercept Deviation")

########################################################################
########################################################################
########################################################################
# NON CONTINUOUS MODEL
########################################################################
########################################################################
########################################################################

mnc0 <- clmm(ANXFAIL ~ 1 + (1 | CNTSCHID),
             link = "logit", 
             data = pisa_extend_ncont)
summary(mnc0)

# ICC for ordinal logit: same formula as binary logit
# Level-1 variance is fixed to pi^2/3 (logistic distribution)
var_u0_logit <- as.data.frame(VarCorr(mnc0))[1, 1]
var_u0_logit
icc_logit    <- var_u0_logit / (var_u0_logit + pi^2 / 3)
cat("\nICC (logit null model):", round(icc_logit, 3), "\n")
#1.5% of variation is due to school

####################################
# Random effects
####################################

#random intercept only
mnc1 <- clmm(ANXFAIL ~ 
               ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z +
               LANGMATCH:SCH_PERCDIS_z + 
               (1 | CNTSCHID),
             link = "logit", 
             data = pisa_extend_ncont)

summary(mnc1)

# uncorrelated: test whether slope variance for TEACHSUP is needed 
mnc1a <- clmm(ANXFAIL ~ 
                 ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z +
                 LANGMATCH:SCH_PERCDIS_z + 
                 (1 | CNTSCHID) + (-1 + TEACHSUP_z | CNTSCHID),
               link = "logit", 
               data = pisa_extend_ncont)

# LRT: is the random slope for teachsup needed?
anova(mnc1, mnc1a) # significant at alpha = .05
# the effect of TEACHSUP_z varies across schools

# uncorrelated: test whether slope variance for ESCS is needed 
mnc1b <- clmm(ANXFAIL ~ 
                ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z +
                LANGMATCH:SCH_PERCDIS_z + 
                (1 | CNTSCHID) 
              + (-1 + ESCS_z | CNTSCHID)
              + (-1 + TEACHSUP_z | CNTSCHID),
              link = "logit", 
              data = pisa_extend_ncont)
anova(mnc1a, mnc1b)
VarCorr(mnc1b)

# correlated random intercept and slope: 
mnc1c <- clmm(ANXFAIL ~ 
                ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z +
                LANGMATCH:SCH_PERCDIS_z + 
                (1 + TEACHSUP_z | CNTSCHID),
              link = "logit", 
              data = pisa_extend_ncont)

anova(mnc1a, mnc1c) # methodology here is based on exercise 07

# check random intercept needed
mnc1d <- clmm(ANXFAIL ~
                ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z +
                LANGMATCH:SCH_PERCDIS_z +
                (-1 + TEACHSUP_z | CNTSCHID),
              link = "logit",
              data = pisa_extend_ncont)

anova(mnc1d, mnc1a) # significant at alpha = .05, random intercept needed

mnc2 <- mnc1a

####################################
# Interaction effects
####################################

# remove interaction effect for LANGMATCH
mnc2a <- clmm(ANXFAIL ~ 
                ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + SCH_PERCDIS_z +
                (1 | CNTSCHID) + (-1 + TEACHSUP_z | CNTSCHID),
              link = "logit", 
              data = pisa_extend_ncont)

anova(mnc2a, mnc2) # not significant at alpha = .05, drop interaction effect

mnc3 <- mnc2a

####################################
# Main effects
####################################
# drop LANGMATCH
mnc3a <- clmm(ANXFAIL ~ 
                ESCS_z + TEACHSUP_z + GENDER + SCH_PERCDIS_z + 
                (1 | CNTSCHID) + (-1 + TEACHSUP_z | CNTSCHID),
              link = "logit", 
              data = pisa_extend_ncont)

anova(mnc3a, mnc3) # significant at alpha = .05

# drop SCH_PERCDIS_z
mnc3b <- clmm(ANXFAIL ~ 
                ESCS_z + TEACHSUP_z + GENDER + LANGMATCH + 
                (1 | CNTSCHID) + (-1 + TEACHSUP_z | CNTSCHID),
              link = "logit", 
              data = pisa_extend_ncont)

anova(mnc3b, mnc3) # not significant

# drop TEACHSUP 
mnc3c <- clmm(ANXFAIL ~ 
                ESCS_z + GENDER + LANGMATCH +
                (1 | CNTSCHID) + (-1 + TEACHSUP_z | CNTSCHID),
              link = "logit", 
              data = pisa_extend_ncont)
anova(mnc3c, mnc3b) # significant at alpha = .05

# drop ESCS_z 
mnc3e <- clmm(ANXFAIL ~ 
                TEACHSUP_z + GENDER + LANGMATCH +
                (1 | CNTSCHID) + (-1 + TEACHSUP_z | CNTSCHID),
              link = "logit", 
              data = pisa_extend_ncont)

anova(mnc3e, mnc3b) # significant at alpha = .05

# final model
mncFinal <- clmm(ANXFAIL ~ 
                   ESCS_z + TEACHSUP_z + GENDER + LANGMATCH +
                   (1 | CNTSCHID) + (-1 + TEACHSUP_z | CNTSCHID),
                 link = "logit", 
                 data = pisa_extend_ncont)

summary(mncFinal)
####################################
# Residuals
####################################
sjPlot::tab_model(mncFinal)

# Convergence: ensure the model converged (no warnings about convergence failure)
# Complete separation: if any cell of the contingency table between the response 
#     and a predictor is empty, the model may produce unreliable thresholds
# Proportional odds assumption: 
# Fit C-1 = 3 binary logistic models and compare slope estimates.
# If the slopes agree, proportional odds is plausible.
pisa_extend_ncont$ANXFAIL_ge2 <- as.integer(pisa_extend_ncont$ANXFAIL >= "Agree")   # Y >=2
pisa_extend_ncont$ANXFAIL_ge3 <- as.integer(pisa_extend_ncont$ANXFAIL >= "Disagree")   # Y >=3
pisa_extend_ncont$ANXFAIL_ge4 <- as.integer(pisa_extend_ncont$ANXFAIL >= "Strongly Disagree") # Y >= 4
# Levels: Strongly Agree < Agree < Disagree < Strongly Disagree

lm_thresh1 <- glm(ANXFAIL_ge2 ~ ESCS_z + TEACHSUP_z + GENDER + LANGMATCH,
                  family = binomial(link = "logit"), data = pisa_extend_ncont)
lm_thresh2 <- glm(ANXFAIL_ge3 ~ ESCS_z + TEACHSUP_z + GENDER + LANGMATCH,
                  family = binomial(link = "logit"), data = pisa_extend_ncont)
lm_thresh3 <- glm(ANXFAIL_ge4 ~ ESCS_z + TEACHSUP_z + GENDER + LANGMATCH,
                  family = binomial(link = "logit"), data = pisa_extend_ncont)

cat("\nProportional odds check - slope estimates per threshold:\n")
comp <- rbind(
  threshold_1_2 = coef(lm_thresh1)[-1],
  threshold_2_3 = coef(lm_thresh2)[-1], 
  threshold_3_4 = coef(lm_thresh3)[-1]
)
print(round(comp, 3))
# Slopes should be similar across rows if proportional odds holds.
# Large discrepancies suggest the assumption may be violated.

# Transpose so predictors are rows 
comp_table <- as.data.frame(t(comp))
# Add absolute range across thresholds 
comp_table$Range <- apply(comp, 2, max) - apply(comp, 2, min) 

print(round(comp_table, 3))