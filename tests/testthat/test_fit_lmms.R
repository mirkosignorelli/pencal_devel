# Tests for Step 1 LMM estimation

testthat::test_that("Step 1 LMM estimation max.ymissing error works", {
  set.seed(1234)
  # Number of outcomes
  p <- 4
  # Simulate data
  simdata <- simulate_prclmm_data(n = 100, p = p, p.relev = 2,
                                  seed = 123, t.values = c(0, 0.2, 0.5, 1, 1.5, 2))
  
  # Define names for biomarkers
  y.names <- paste0("marker", 1:p)
  
  max.ymissing <- 1.1
  
  testthat::expect_error(fit_lmms(y.names = y.names,
                                  fixefs = ~ age, ranefs = ~ age | id,
                                  long.data = simdata$long.data,
                                  surv.data = simdata$surv.data,
                                  t.from.base = t.from.base,
                                  n.boots = 0, n.cores = 1,
                                  max.ymissing = max.ymissing))
})

testthat::test_that("Step 1 LMM estimation result matches reference", {
  set.seed(1234)
  # Number of outcomes
  p <- 4
  # Simulate data
  simdata <- simulate_prclmm_data(n = 100, p = p, p.relev = 2,
                                 seed = 123, t.values = c(0, 0.2, 0.5, 1, 1.5, 2))
  
  # Define names for biomarkers
  y.names <- paste0("marker", 1:p)
  
  # Estimate the LMMs and test against snapshot
  testthat::expect_snapshot({
    fit_lmms(y.names = y.names,
             fixefs = ~ age, ranefs = ~ age | id,
             long.data = simdata$long.data,
             surv.data = simdata$surv.data,
             t.from.base = t.from.base,
             n.boots = 0, n.cores = 1)
    # We need to remove the identifier of the R environment from the call info 
    # because it differs for every run of the test
    }, transform = function(x) gsub("<.+>", replacement = "", x = x))
})
