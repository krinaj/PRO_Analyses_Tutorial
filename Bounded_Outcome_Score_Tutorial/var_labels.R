# Define questionnaire scales
  qs_scales <- tribble(
    ~code,                   ~label, ~n_items, ~range,        ~FLAG, ~is_func,
    "QL2",          "Global Health",        2,      6,     c(29, 30),   FALSE,
    "PF2",   "Physical Functioning",        5,      3,        c(1:5),    TRUE, 
    "RF2",       "Role Functioning",        2,      3,       c(6, 7),    TRUE,
     "EF",  "Emotional Functioning",        4,      3,      c(21:24),    TRUE, 
     "CF",  "Cognitive Functioning",        2,      3,     c(20, 25),    TRUE, 
     "SF",     "Social Functioning",        2,      3,     c(26, 27),    TRUE,
     "FA",       "Fatigue Symptoms",        3,      3, c(10, 12, 18),   FALSE,
     "NV",        "Nausea Symptoms",        2,      3,     c(14, 15),   FALSE, 
     "PA",          "Pain Symptoms",        2,      3,      c(9, 19),   FALSE, 
     "DY",      "Dyspnoea Symptoms",        1,      3,          c(8),   FALSE, 
     "SL",      "Insomnia Symptoms",        1,      3,         c(11),   FALSE, 
     "AP", "Appetite Loss Symptoms",        1,      3,         c(13),   FALSE, 
     "CO",  "Constipation Symptoms",        1,      3,         c(16),   FALSE, 
     "DI",     "Diarrhoea Symptoms",        1,      3,         c(17),   FALSE,
     "FI", "Financial Difficulties",        1,      3,         c(28),   FALSE
  )
  qs_scales$SCALE <- with(qs_scales, as.integer(factor(code, levels = code)))
  
# Define categorical covariates
  cov_definition <- list(
    STDIDf = c(
      "NCT00081796" = 2004
    ),
    SEXf = c(
      "Male" = 1,
      "Female" = 2
    ),
    RACEf = c(
      "White" = 1,
      "Black/African American" = 2,
      "Asian/Pacific Islander" = 3,
      "Latino/Hispanic" = 4,
      "Other" = 5
    ),
    MENOSf = c(
      "Premenopausal" = 1,
      "Other" = 2
    ), 
    ECOGf = c(
      "0" = 0,
      "1" = 1,
      "2-3" = 2,
      "Missing" = NA_real_
    ),
    DTHf = c(
      "Censored" = 0, 
      "Death" = 1,
      "Missing" = NA_real_
    ),
    CHEMOf = c(
      "No" = 0,
      "Yes" = 1,
      "Unknown" = NA_real_
    ),
    SURGERYf = c(
      "No" = 0,
      "Yes" = 1,
      "Unknown" = NA_real_
    ),
    RADIOf = c(
      "No" = 0,
      "Yes" = 1,
      "Unknown" = NA_real_
    ),
    MHDEPf = c(
      "No" = 0,
      "Yes" = 1,
      "Unknown" = NA_real_
    ),
    MDANXf = c(
      "No" = 0,
      "Yes" = 1,
      "Unknown" = NA_real_
    ),
    ERSf = c(
      "Negative" = 0,
      "Positive" = 1,
      "Unknown" = NA_real_
    ),
    PGRSf = c(
      "Negative" = 0,
      "Positive" = 1,
      "Unknown" = NA_real_
    )
  )

# Convert categorical covariates to dataframe with enframe
  catcov_labels <- cov_definition %>% 
    imap(enframe) %>%
    # Change name of value column from "value" to COV (COVf without the 'f')
    map(\(cov) rename(cov, !!str_remove(names(cov)[[1]], "f") := value)) %>%
    # For each covariate (i), change COVf column containing category strings (ii), to
    # type factor while maintaining the desired order of categories (iii).
    map(\(cov) {  # (i) for each covariate
      mutate(cov, across(ends_with("f"), \(catstring) {  # (ii) change COVf column
        factor(catstring, levels = unique(cov[[1]]))  # (iii) to type factor
      })
    )})
   