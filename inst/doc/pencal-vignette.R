## ---- echo = F----------------------------------------------------------------
# CRAN limits number of available cores to 2 although it has 32
# in pencal this triggers warnings about the number of cores used
options(warn=-1)
# suppress warnings related to number of cores when rendering the vignette

## ---- echo = F, out.width = "650px"-------------------------------------------
knitr::include_graphics("PRC_diagram.jpg")

## ---- eval=FALSE, echo=TRUE, results='asis'-----------------------------------
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install('pencal')

## ---- eval=FALSE, echo=TRUE, results='asis'-----------------------------------
#  # step 1: install survcomp
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install('survcomp')
#  # step 2: install pencal and CRAN dependencies
#  install.packages('pencal')

## ---- eval=TRUE, echo=TRUE, results='asis'------------------------------------
library(pencal)

## ----simulate, cache = F------------------------------------------------------
set.seed(1234)
p = 10
simdata = simulate_prclmm_data(n = 100, p = p, p.relev = 5, 
                               lambda = 0.2, nu = 1.5,
                               seed = 1234, t.values = c(0, 0.2, 0.5, 1, 1.5, 2))
ls(simdata)

## ----view_longdata, fig.height=4, fig.width=5---------------------------------
# view the dataset in long format
head(simdata$long.data)
# visualize the trajectories for a randomly picked biomarker
# (the code in the if statement below relies on the ptmixed package)
if ('ptmixed' %in% rownames(installed.packages())) {
  library(ptmixed)
  ptmixed::make.spaghetti(x = age, y = marker5, 
               id = id, group = id,
               data = simdata$long.data, 
               margins = c(4, 4, 2, 2),
               legend.inset = - 1)
}


## ----view_survdata, fig.height=4.5, fig.width=4-------------------------------
# view the dataset with the survival data
head(simdata$surv.data)
# what is the proportion of censoring in this dataset?
simdata$censoring.prop
# visualize an estimate of the survival function
library(survival)
library(survminer)
surv.obj = survival::Surv(time = simdata$surv.data$time, 
                event = simdata$surv.data$event)
kaplan = survival::survfit(surv.obj ~ 1,  
                  type="kaplan-meier")
survminer::ggsurvplot(kaplan, data = simdata$surv.data)

## ---- eval=F------------------------------------------------------------------
#  n.cores = parallel::detectCores()

## -----------------------------------------------------------------------------
n.cores = 2

## ----step1, cache = F---------------------------------------------------------
y.names = paste('marker', 1:p, sep = '')
step1 = fit_lmms(y.names = y.names, 
                 fixefs = ~ age, ranefs = ~ age | id, 
                 long.data = simdata$long.data, 
                 surv.data = simdata$surv.data,
                 t.from.base = t.from.base,
                 n.boots = 10, n.cores = n.cores)

## -----------------------------------------------------------------------------
ls(step1)
# estimated LMM for marker1:
step1$lmm.fits.orig[1]

## -----------------------------------------------------------------------------
getlmm(step1, 'marker3', 'betas')
getlmm(step1, 'marker3', 'tTable')
getlmm(step1, 'marker3', 'variances')

## ----step2, cache = F---------------------------------------------------------
step2 = summarize_lmms(object = step1, n.cores = n.cores)

## -----------------------------------------------------------------------------
ls(step2)
# view predicted random effects for the first two markers
step2$ranef.orig[1:5, 1:4]

## ----step3, cache = F---------------------------------------------------------
step3 = fit_prclmm(object = step2, surv.data = simdata$surv.data,
                   baseline.covs = ~ baseline.age,
                   penalty = 'ridge', n.cores = n.cores)

## -----------------------------------------------------------------------------
print(step3)

## -----------------------------------------------------------------------------
ls(step3)
class(step3$pcox.orig)
library(glmnet)
t(as.matrix(coef(step3$pcox.orig)))

## -----------------------------------------------------------------------------
preds = survpred_prclmm(step1, step2, step3, times = c(1, 2, 3))
ls(preds)
head(preds$predicted_survival)

## ----cbocp, cache = F---------------------------------------------------------
cbocp = performance_prc(step2, step3, times = c(1, 2, 3), 
                   n.cores = n.cores)
# C index estimates:
cbocp$concordance
# time-dependent AUC estimates:
cbocp$tdAUC

