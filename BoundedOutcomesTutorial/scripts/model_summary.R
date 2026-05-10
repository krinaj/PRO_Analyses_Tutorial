# Initiate sourcing of script
  message("Sourcing model_summary.R")
if (!("tidyverse" %in% .packages())) { 
  message("Loading model_summary dependencies: tidyverse")
  library(tidyverse) 
}

#' @title  model_summary
#' @description  Standardised summary of NONMEM output.
#'
#' @param files a character vector detailing the file names to search for 
#'   when reading in NONMEM output. By default searches for file extensions
#'   .cov, .ext, .lst and .shk. Should be defined if more than one file with
#'   these extensions is present in the working directory.
#' @param output_format a character vector detailing which outputs file formats
#'   output should be saved to. By default saves to all three formats: .txt for
#'   human-readable plain text summary, .rds and .Rdata for full output.
#' @param path a character string of the path to the NONMEM output to summarise; 
#'   the default corresponds to the working directory, `getwd()`. 
#' @param print_rounded a logical value defining whether values output to R
#'   console should be rounded. Not rounded by default. Does not impact rounding
#'   of output formats.
#' @param condnum_sqrt a logical value defining whether the square root of the
#'   condition number should be provided in output.
#' @param param_labels a logical value defining whether parameter labels from
#'   the model control stream should be parsed. Parsing occurs according to 
#'   `theta_prefix`, `omega_prefix` and `sigma_prefix` arguments.
#' @param theta_prefix a character string defining the character pattern used
#'   in the model control stream prior to listing the name of a parameter in 
#'   the $THETA block. Must be different from parameters in other blocks.
#' @param omega_prefix a character string defining the character pattern used
#'   in the model control stream prior to listing the name of a parameter in 
#'   the $OMEGA block. Must be different from parameters in other blocks.
#' @param sigma_prefix a character string defining the character pattern used
#'   in the model control stream prior to listing the name of a parameter in 
#'   the $SIGMA block. Must be different from parameters in other blocks.
#' @param trans_labels a logical value defining whether parameter labels contain
#'   information about transformed parameter values. Parsing occurs according to
#'   `log_prefix` argument. 
#'
#' @rdname model_summary
#' @author Jessica Wojciechowski, \email{jessica.wojciechowski@@pfizer.com}
#' @author Jim Hughes, \email{jim.hughes@@pfizer.com}
#' 
#' @export 

model_summary <- function(files = c(".cov", ".ext", ".lst", ".shk"), 
  output_format = c(".txt", ".rds", ".Rdata"), path = ".", 
  print_rounded = FALSE, condnum_sqrt = TRUE, param_labels = FALSE,
  theta_prefix = "; TV", omega_prefix = "; PPV", sigma_prefix = "; ERR",
  trans_labels = FALSE, log_prefix = "LN") {
### Prepare workspace - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# Change working directory if required
  prev_path <- getwd()
  setwd(path)

### Check what files are available in the repository
  output_files <- files %>%
    paste0("$") %>%
  # Determine whether any files with extension "ext" are available
    map(str_subset, string = dir()) %>%
  # Return a message if a file is not in the directory
    walk2(files, function(file, ext) {
      if (length(file) == 0) message(paste("A file containing", ext, 
        "could not be found in the repository."))
    }) %>%
    unlist()
  
### Read in the available files into the workspace and extract values
# Suppress warnings here in case files with certain extensions are missing
# Dev note: regex pattern (?=.) is a lookahead, it makes sure to match text 
# with . directly after it (in case below this is ".ext")
  output_data <- list(
    run = na.omit(suppressWarnings(str_extract(output_files, ".*(?=\\.lst$)")))[[1]],
    cov = read_covdata(suppressWarnings(str_subset(output_files, ".cov$"))),
    ext = read_extdata(suppressWarnings(str_subset(output_files, ".ext$"))),
    lst = read_lstdata(suppressWarnings(str_subset(output_files, ".lst$")), 
      labels = param_labels, theta_prefix = theta_prefix, 
      omega_prefix = omega_prefix, sigma_prefix = sigma_prefix),
    shk = read_shkdata(suppressWarnings(str_subset(output_files, ".shk$"))),
    stepid = read_metadata(str_subset(dir(), "stepID\\.txt"))
  )
# Extract values of interest from .ext, .shk and .lst
  run_dims <- analysis_dimensions(output_data)
  ofv_value <- objective_function(output_data, npar = run_dims$npar)
  cond_number <- condition_number(output_data, square_root = condnum_sqrt)
  final_param <- final_parameters(output_data, trans = trans_labels, 
    log_prefix = log_prefix)
  covar_mat <- covariance_matrix(output_data, final_param)
  
### Generate list of final output
  results_summary <- list(
    "Run" = output_data$run, 
    "AnalysisStepID" = output_data$stepid, 
    "EstimationMethod" =  output_data$lst$est_method, 
    "Minimisation" = output_data$lst$est_status, 
    "CovarianceStep" = output_data$lst$cov_status, 
    "Dimensions" = run_dims,
    "ObjectiveFunctionValue" = ofv_value, 
    "ConditionNumber" = cond_number, 
    "Theta" = as.data.frame(final_param$theta), 
    "Omega" = as.data.frame(final_param$omega), 
    "Covariance" = as.data.frame(final_param$ocorr), 
    "Sigma" = as.data.frame(final_param$sigma),
    "SigmaCov" = as.data.frame(final_param$scorr),
    "CovarianceMatrix" = covar_mat
  )
 
### Create rounded versions of parameter data frames 
  if (".txt" %in% output_format | print_rounded) {
  # Round each set of parameter types
    rounded_param <- map(final_param, function(param) {
      mutate_if(param, names(param) %in% c("RSE", "Value", "Norm", "SD", "Corr", "Shrinkage"), 
      function(x) {
        round_which <- which(!str_detect(x, "^FIX$|^No\\s"))
        x[round_which] <- sprintf(as.double(x[round_which]), fmt ="%0.3g")
        return(x)
      })
    })
  # Define rounded results summary list
    rounded_summary <- list_modify(results_summary,
      Minimisation = str_pad(output_data$lst$est_status, width = 50, side = "right"),
      Theta = as.data.frame(rounded_param$theta), 
      Omega = as.data.frame(rounded_param$omega), 
      Covariance = as.data.frame(rounded_param$ocorr), 
      Sigma = as.data.frame(rounded_param$sigma),
      SigmaCov = as.data.frame(rounded_param$scorr),
      CovarianceMatrix = zap()
    ) # rounded_summary
  }
  
### Save results to output files
# Save results summary to .rds for use with mrgsolve and report tables
  results_file <- paste0(output_data$run, "_summary")
  if (".rds" %in% output_format) {
    saveRDS(results_summary, file = paste0(results_file, ".rds"))
  }
# Save rounded summary to .txt to view in improve
  if (".txt" %in% output_format) {
    capture.output(rounded_summary, file = paste0(results_file, ".txt"))
  }
# Save results summary to .Rdata for QC in R and compatability with old code
  if (".Rdata" %in% output_format) {
    results.l <- results_summary
    save(file = paste0(results_file, ".Rdata"), results.l, version = 2)
  }
  
### Clean up and provide output
# Return working directory to it's original location
  setwd(prev_path)
# Provide output but do not print to console
  if (print_rounded) {
    rounded_summary
  } else {
    results_summary
  }
}  # model_summary  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


#' @title  read_covdata
#' @description  Read a NONMEM .cov file into a tibble.
#'
#' @param file_name a character string defining the a path to a .cov file.
#'
#' @rdname read_covdata
#' @author Jessica Wojciechowski, \email{jessica.wojciechowski@@pfizer.com}
#' @author Jim Hughes, \email{jim.hughes@@pfizer.com}
#' 
#' @export 

read_covdata <- function(file_name) {
# Check number of files passed to function  - - - - - - - - - - - - - - - - - - 
  if (length(file_name) == 0) {
    return(NULL)
  } else if (length(file_name) > 1) {
    stop("More than 1 .cov file detected. Please specify file name with `files` argument")
  }
# Determine how many estimation methods are available in the file
# Extract values for final estimation method
  covdata <- file_name %>%
    readLines() %>%
    str_which("TABLE NO.") %>%
    {if_else(length(.) == 1, 1L, tail(., 1))} %>%
    {read_table(file_name, skip = ., col_types = cols(
      .default = col_double(),
        NAME = col_character())
    )} %>%
    select(-NAME)
# Remove parameters not estimated
# Remove columns that are completely zero
  keep_cols <- covdata %>%
    map(sum) %>%
    map_lgl(~ .x != 0)
# Remove rows that are completely zero
  clean_covdata <- covdata[keep_cols] %>%
    split(seq_len(nrow(.))) %>%
    map_dfr(function(row) mutate(row, rowsum = sum(row))) %>%
    filter(rowsum != 0) %>%
    select(-rowsum)
}  # read_covdata - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


#' @title  read_extdata
#' @description  Read a NONMEM .ext file into a tibble.
#'
#' @param file_name a character string defining the a path to a .ext file.
#'
#' @rdname read_extdata
#' @author Jessica Wojciechowski, \email{jessica.wojciechowski@@pfizer.com}
#' @author Jim Hughes, \email{jim.hughes@@pfizer.com}
#' 
#' @export 

read_extdata <- function(file_name) {
# Check number of files passed to function  - - - - - - - - - - - - - - - - - - 
  if (length(file_name) == 0) {
    stop
  } else if (length(file_name) > 1) {
    stop("More than 1 .ext file detected. Please specify file name with `files` argument")
  }
# Determine how many estimation methods are available in the file
# Extract values for final estimation method
  extdata <- file_name %>%
    readLines() %>%
    str_which("TABLE NO.") %>%
    {if_else(length(.) == 1, 1L, tail(., 1))} %>%
    {read_table(file_name, skip = ., col_types = cols(.default = col_double()))} %>%
    filter(ITERATION %in% c(
      -1000000000, # Final parameter estimates
      -1000000001, # Standard errors of final parameter estimates
      -1000000002, # Eigenvalues
      -1000000004, # Correlation matrix for random effects
      -1000000005, # Standard errors of correlation matrix for random effects
      -1000000006  # Identifier for fixed parameters
    ))
}  # read_extdata - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


#' @title  read_lstdata
#' @description  Read and parse a NONMEM .lst file.
#'
#' @param file_name a character string defining the a path to a .lst file.
#' @param labels a logical value defining whether parameter labels from
#'   the model control stream should be parsed. Parsing occurs according to 
#'   `theta_prefix`, `omega_prefix` and `sigma_prefix` arguments.
#' @param theta_prefix  a character string defining the character pattern used
#'   in the model control stream prior to listing the name of a parameter in 
#'   the $THETA block. Must be different from parameters in other blocks.
#' @param omega_prefix a character string defining the character pattern used
#'   in the model control stream prior to listing the name of a parameter in 
#'   the $OMEGA block. Must be different from parameters in other blocks.
#' @param sigma_prefix a character string defining the character pattern used
#'   in the model control stream prior to listing the name of a parameter in 
#'   the $SIGMA block. Must be different from parameters in other blocks.
#'
#' @rdname read_lstdata
#' @author Jessica Wojciechowski, \email{jessica.wojciechowski@@pfizer.com}
#' @author Jim Hughes, \email{jim.hughes@@pfizer.com}
#' 
#' @export 

read_lstdata <- function(file_name, labels = TRUE, theta_prefix = "; TV", 
  omega_prefix = "; PPV", sigma_prefix = "; ERR") {
# Check number of files passed to function
  if (length(file_name) == 0) {
    return(NULL)
  } else if (length(file_name) > 1) {
    stop("More than 1 .lst file detected. Please specify file name with `files` argument")
  }
# Read in file
  lstdata <- readLines(file_name)
# Extract information from file
  est_method <- lstdata %>%
    str_subset("#METH") %>%
    str_replace_all(" #METH: ","")
# Determine parameter estimation status
  est_status <- lstdata[str_which(lstdata, "#TERM") + 1] %>%
    str_trim(side = "left")
# Determine covariance estimation status
  cov_status <- lstdata %>%
    str_detect("STANDARD ERROR OF ESTIMATE") %>%
    any() %>%
    if_else("Passed", "Failed")
# Determine number of observations
  sample_messages <- c("TOT. NO. OF OBS RECS:", "TOT. NO. OF INDIVIDUALS:")
  sample_size <- sample_messages %>%
    map(str_subset, string = lstdata) %>%
    set_names("nobs", "nid") %>%
    map2(sample_messages, function(line, message) str_replace_all(line, message, "")) %>%
    map(str_trim, side = "left") %>%
    map_dfr(function(n) as.double(unique(n)))
# Return extracted information as a list
  lst_output <- list(
    est_method = est_method,
    est_status = est_status,
    cov_status = cov_status,
    sample_size = sample_size
  )
# Extract parameter labels
  if (labels) {
  # Extract theta labels
    lst_output$theta <- lstdata %>%
      str_subset(theta_prefix) %>%
      str_subset("^(\\$THETA)*\\s*\\(*[-]*\\d+") %>%  # ensure THETAs aren't commented out
      str_extract(paste0("(?<=", theta_prefix, ").*")) %>%
      # keep(str_detect(., "^[A-Za-z]")) # keep only values that start with letters
      str_trim()
  # Extract omega labels
    omega_labels <- lstdata %>%
      str_subset(omega_prefix) %>%
      str_subset("^(\\$OMEGA)*(\\s*BLOCK\\(\\d+\\))*\\s*(([-]*\\d+)|(SAME))") %>%  # ensure OMEGAs aren't commented out
      str_extract(paste0("(?<=", omega_prefix, ").*")) %>%
      str_trim()
  # Identify positions of diagonal and off-diagonal elements
    omega_index <- 1:length(omega_labels)*(1:length(omega_labels) + 1)/2
    omat_index <- seq_len(max(omega_index))
    ocorr_index <- omat_index[!omat_index %in% omega_index]
  # Generate labels for off-diagonal elements
    ocorrLabels <- map_chr(ocorr_index, function(o) {
      lower <- max(omega_index[omega_index < o])
      upper <- min(omega_index[omega_index > o])
      lowerlabel <- omega_labels[o - lower]
      upperlabel <- omega_labels[omega_index == upper]
      corrlabel <- paste0("R", lowerlabel, "-", upperlabel)
    }) # map_chr
  # Combine omega and correlation labels
    lst_output$omega <- omat_index
    lst_output$omega[omega_index] <- paste0("IIV", omega_labels)
    lst_output$omega[ocorr_index] <- ocorrLabels
  # Extract sigma labels
    sigma_labels <- lstdata %>%
      str_subset(sigma_prefix) %>%
      str_subset("^(\\$SIGMA)*(\\s*BLOCK\\(\\d+\\))*\\s*[-]*\\d+") %>%  # ensure SIGMAs aren't commented out
      str_extract(paste0("(?<=", sigma_prefix, ").*")) %>%
      str_trim()
  # If there are Sigma labels
    if (length(sigma_labels != 0)) {
    # Identify positions of diagonal and off-diagonal elements
      sigma_index <- 1:length(sigma_labels)*(1:length(sigma_labels) + 1)/2
      smat_index <- seq_len(max(sigma_index))
      scorr_index <- smat_index[!smat_index %in% sigma_index]
    # Generate labels for off-diagonal elements
      scorr_labels <- map_chr(scorr_index, function(o) {
        lower <- max(sigma_index[sigma_index < o])
        upper <- min(sigma_index[sigma_index > o])
        lowerlabel <- sigma_labels[o - lower]
        upperlabel <- sigma_labels[sigma_index == upper]
        corrlabel <- paste0("R", lowerlabel, "-", upperlabel)
      })
    # Combine sigma and correlation labels
      lst_output$sigma <- smat_index
      lst_output$sigma[sigma_index] <- sigma_labels
      lst_output$sigma[scorr_index] <- scorr_labels
  # If there are no Sigma labels
    } else {
    # Create default sigma label for the zero-sigma NONMEM output
      lst_output$sigma <- "ERRRES"
    } # if/else
  }
  return(lst_output)
}  # read_lstdata - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


#' @title  read_shkdata
#' @description  Read a NONMEM .shk file into a tibble.
#'
#' @param file_name a character string defining the a path to a .shk file.
#'
#' @rdname read_shkdata
#' @author Jessica Wojciechowski, \email{jessica.wojciechowski@@pfizer.com}
#' @author Jim Hughes, \email{jim.hughes@@pfizer.com}
#' 
#' @export 

read_shkdata <- function(file_name) {
# Check number of files passed to function
  if (length(file_name) == 0) {
    return(NULL)
  } else if (length(file_name) > 1) {
    stop("More than 1 .shk file detected. Please specify file name with `files` argument")
  }
# Determine how many estimation methods are available in the file
# Extract values for final estimation method
  shkdata <- readLines(file_name) %>%
    str_which("TABLE NO.") %>%
    {if_else(length(.) == 1, 1L, tail(., 1))} %>%
    {read_table(file_name, skip = ., col_types = cols(.default = col_double()))} %>%
    # select(-TYPE, -SUBPOP) %>%
    # rename(TYPE = X2, SUBPOP = X4) %>%
    filter(TYPE %in% c(4, 5)) # ETA shrinkage SD, EPS shrinkage SD
}


#' @title  read_metadata
#' @description  Read and parse improve metadata.
#'
#' @param file_name a character string defining the path to a metadata file.
#'
#' @rdname read_stdout
#' @author Jim Hughes, \email{jim.hughes@@pfizer.com}
#' 
#' @export 

read_metadata <- function(file_name) {
# Note: This function currently only handles simple metadata collection of the
# Analysis Step entity ID. Functionality may grow over time.
# Check number of files passed to function
  if (length(file_name) == 0) {
    return(NULL)
  }
# Read in file containing metadata
  readLines(file_name)
}  # read_metadata  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
#' @title analysis_dimensions
#' @description  Determines dimensions (including sample size) of the analysis
#' 
#' @param output_data an output object from the `model_summary` function.
#'
#' @rdname analysis_dimensions
#' @author Jim Hughes, \email{jim.hughes@@pfizer.com}
#' 

analysis_dimensions <- function(output_data) {
# Extract sample_size object derived from .lst
  dims_data <- output_data$lst$sample_size
# Determine number of parameters estimated from .cov file (if available)
  if (!is.null(output_data$cov)) {
    nparam <- ncol(output_data$cov)
  } else {
    nparam <- NA_real_
  }
# Add number of parameters to output and return dataframe
  dims_data$npar <- nparam
  return(dims_data)
}  # analysis_dimensions  - - - - - - - - - - - - - - - - - - - - - - - - - - - 


#' @title  objective_function
#' @description  Calculate the objective function value
#'
#' @param output_data an output object from the `model_summary` function.
#' @param npar number of parameters estimated during analysis as determined by
#'   the `analysis_dimensions` function.
#'
#' @rdname objective_function
#' @author Jessica Wojciechowski, \email{jessica.wojciechowski@@pfizer.com}
#' @author Jim Hughes, \email{jim.hughes@@pfizer.com}
#' 

objective_function <- function(output_data, npar) {
# Isolate sample size and number of parameters (if .cov available) from output
  nobs <- output_data$lst$sample_size$nobs
# Extract objective function value
  output_data$ext %>%
    filter(ITERATION == -1000000000) %>%
    select(OBJ) %>%
    rename(OFV = OBJ) %>%
    mutate(AIC = 2*npar + OFV) %>%
    mutate(BIC = npar*log(nobs) + OFV)
}  # objective_function - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  

#' @title  condition_number
#' @description  Calculate the condition number
#'
#' @param output_data an output object from the `model_summary` function.
#' @param square_root a logical value defining whether the square root of the
#'   condition number should be provided in output.
#'
#' @rdname condition_number
#' @author Jessica Wojciechowski, \email{jessica.wojciechowski@@pfizer.com}
#' @author Jim Hughes, \email{jim.hughes@@pfizer.com}
#' 

condition_number <- function(output_data, square_root = TRUE) {
# Pull eigenvalues from .ext output only if .cov is available
  if (!is.null(output_data$cov)) {
    eigen_values <- output_data$ext %>%
      filter(ITERATION == -1000000002) %>%
      select(-ITERATION,-OBJ) %>%
      c() %>%
      unname() %>%
      unlist() %>%
      `[`(. != 0)
  # Calculate condition number (ratio of the highest to lower eigenvalues)
  # Calculate square root of condition number if requested
    cond_number <- max(eigen_values)/min(eigen_values)
    if (square_root) cond_number <- sqrt(cond_number)
  } else {
    cond_number <- "Covariance Step Failed"
  } # if/else
# Return condition number
  return(cond_number)
}  # condition_number - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  

#' @title  final_parameters
#' @description  Generate table summaries of final parameter estimates.
#'
#' @param output_data an output object from the `model_summary` function.
#'
#' @rdname final_parameters
#' @author Jessica Wojciechowski, \email{jessica.wojciechowski@@pfizer.com}
#' @author Jim Hughes, \email{jim.hughes@@pfizer.com}
#' 

final_parameters <- function(output_data, trans = "TRUE", log_prefix = "LN") {
# Extract parameter values from .ext file
  all_param <- output_data$ext %>%
    filter(ITERATION == -1000000000) %>%
    select(-ITERATION, -OBJ) %>%
    pivot_longer(cols = everything(), names_to = "Number", values_to = "Value") %>%
    mutate(is_log = FALSE)
# Add parameter labels if available
  if (!is.null(output_data$lst$theta)) {
    param_labels <- with(output_data$lst, c(theta, sigma, omega))
    all_param <- mutate(all_param, Name = param_labels)
  # Perform transformations defined in the label if available
    if (trans) {
      log_regex <- paste0("^", log_prefix)
      all_param <- all_param %>% 
        mutate(is_log = {
          str_detect(Number, "^THETA") & str_detect(Name, log_regex)
        }) %>%
        mutate(Norm = if_else(is_log, exp(Value), Value))
    }
  }
# Add correlation values
  corr_param <- output_data$ext %>%
    filter(ITERATION == -1000000004) %>%
    select(-ITERATION, -OBJ) %>%
    pivot_longer(cols = everything(), names_to = "Number", values_to = "Corr") %>%
    left_join(all_param, by = "Number")
# Add relative standard error values (if .cov successful)
  param_rse <- filter(output_data$ext, ITERATION == -1000000001)
  if (nrow(param_rse) != 0) {
    final_param <- param_rse %>%
      select(-ITERATION, -OBJ) %>%
      pivot_longer(cols = everything(), names_to = "Number", values_to = "SE") %>%
      left_join(corr_param, by = "Number") %>%
      mutate(RSE = if_else(is_log, 100*sqrt(exp(SE^2) - 1), abs(SE/Value)*100)) %>% 
      mutate(RSE = ifelse(Number %in% colnames(output_data$cov), RSE, "FIX"))
  } else {
    final_param <- output_data$ext %>%
      filter(ITERATION == -1000000006) %>%
      select(-ITERATION, -OBJ) %>%
      pivot_longer(cols = everything(), names_to = "Number", values_to = "SE") %>%
      left_join(corr_param, by = "Number") %>%
      mutate(RSE = if_else(SE == 1, "FIX", NA_character_))
  }
# Isolate parameter sets and prepare relevant columns
# Prepare theta parameters
  theta_data <- final_param %>%
    filter(str_detect(Number, "^THETA")) %>%
    select(Number, matches("Name"), Value, matches("Norm"), matches("RSE"))
# Prepare omega variance parameters
  if (!is.null(output_data$shk)) {
  # Extract shrinkage values (if .shk available)
    eta_shrinkage <- output_data$shk %>%
      filter(TYPE == 4) %>%
      select(contains("ETA")) %>%
      select(seq_len(nrow(filter(final_param, str_detect(Number, "^OMEGA\\((\\d+),(\\1)\\)$"))))) %>%
      set_names(map_chr(seq_len(ncol(.)), function(n) paste0("OMEGA(", n, ",", n, ")"))) %>%
      pivot_longer(cols = everything(), names_to = "Number", values_to = "Shrinkage")
  } else {
  # Create empty shrinkage dataset (if .shk missing, i.e. MAXEVALS=0)
    eta_shrinkage <- final_param %>%
      filter(str_detect(Number, "^OMEGA\\((\\d+),(\\1)\\)$")) %>%
      mutate(Shrinkage = NA_real_) %>%
      select(Number, Shrinkage)
  }
# Use shrinkage values to isolate relevant omegas
  omega_data <- final_param %>%
    right_join(eta_shrinkage, by = "Number") %>%
    select(Number, matches("Name"), Value, SD = Corr, matches("RSE"), Shrinkage)
# Prepare omega covariance parameters
  if (nrow(omega_data) <= 1) {  # If matrix is 1x1 no off-diagonal elements 
    ocorr_data <- enframe("No off-diagonal elements", name = NULL, value = "Value")
  } else  {  # otherwise
  # Filter for omegas with no shrinkage values (i.e. covariance/correlations)
    ocorr_data <- final_param %>%
      filter(str_detect(Number, "^OMEGA")) %>%
      anti_join(eta_shrinkage, by = "Number") %>%  # filtering join
      select(Number, matches("Name"), Value, Corr, matches("RSE"))
  }
# Prepare sigma variance parameters
  if (!is.null(output_data$shk)) {
  # Extract shrinkage values (if .shk available)
    eps_shrinkage <- output_data$shk %>%
      filter(TYPE == 5) %>%
      select(contains("ETA")) %>%
      select(seq_len(nrow(filter(final_param, str_detect(Number, "^SIGMA\\((\\d+),(\\1)\\)$"))))) %>%
      set_names(map_chr(seq_len(ncol(.)), function(n) paste0("SIGMA(", n, ",", n, ")"))) %>%
      pivot_longer(cols = everything(), names_to = "Number", values_to = "Shrinkage")
  } else {
  # Create empty shrinkage dataset (if .shk missing, i.e. MAXEVALS=0)
    eps_shrinkage <- final_param %>%
      filter(str_detect(Number, "^SIGMA\\((\\d+),(\\1)\\)$")) %>%
      mutate(Shrinkage = NA_real_) %>%
      select(Number,  Shrinkage)
  }
# Use shrinkage values to isolate relevant sigmas
  sigma_data <- final_param %>%
    right_join(eps_shrinkage, by = "Number") %>%
    select(Number, matches("Name"), Value, SD = Corr, matches("RSE"), Shrinkage)
# Prepare sigma covariance parameters
  if (nrow(sigma_data) <= 1) {  # If matrix is 1x1 no off-diagonal elements 
    scorr_data <- enframe("No off-diagonal elements", name = NULL, value = "Value")
  } else  {  # otherwise
  # Filter for sigmas with no shrinkage values (i.e. covariance/correlations)
    scorr_data <- final_param %>%
      filter(str_detect(Number, "SIGMA")) %>%
      anti_join(eps_shrinkage, by = "Number") %>%  # filtering join
      select(Number, matches("Name"), Value, Corr, matches("RSE"))
  }
# Return list of parameter values
  list(
    theta = theta_data,
    omega = omega_data,
    ocorr = ocorr_data,
    sigma = sigma_data,
    scorr = scorr_data
  )
}  # final_parameters - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  

#' @title  covariance_matrix
#' @description  Adds parameter names to covariance matrix.
#'
#' @param output_data an output object from the `model_summary` function.
#' @param final_param an output object from the `final_parameters` function.
#'
#' @rdname covariance_matrix
#' @author Jim Hughes, \email{jim.hughes@@pfizer.com}
#' 

covariance_matrix <- function(output_data, final_param) {
# Extract unnamed covariance matrix
  unnamed_matrix <- output_data$cov
# Return NULL if input covariance matrix is NULL
  if (is.null(unnamed_matrix)) { return(NULL) }
# Extract parameter names from final parameters
  param_names <- final_param %>%
    map(when,
      is.character(.$Value) ~ NULL,  # drop omega and/or sigma covariances if empty
      TRUE ~ .
    ) %>%
    bind_rows()
# If parameter names weren't extracted...
  if (is.null(param_names$Name)) {
  # Return unnamed matrix
    return(unnamed_matrix)
  } else {
  # Otherwise, identify which parameters were estimated and are present in covariance matrix
    rawNames <- names(unnamed_matrix)
    matrix_names <- param_names %>% 
      filter(Number %in% rawNames) %>%
      pull(Name)
  # Name the covariance matrix and return as output
    named_matrix <- set_names(unnamed_matrix, matrix_names)
    return(named_matrix)
  }
}  # covariance_matrix - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

message("The following functions have been loaded:
 - model_summary
 - read_covdata
 - read_extdata
 - read_lstdata
 - read_shkdata
 - read_metadata
 - objective_function
 - condition_number
 - final_parameters
 - covariance_matrix
")
