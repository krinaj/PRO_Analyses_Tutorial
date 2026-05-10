# Initiate sourcing of script
message("Sourcing plot_vpc.R")
walk(c("tidyverse"), \(x) {
  if (!(x %in% .packages())) { 
    message("Loading plot_diagnostics dependencies: ", x)
    library(x, character.only = TRUE) 
  }
})

#' @title  bin_times
#' @description  Bin times for visual predictive check
#' 
#' @param data
#' @param time_var
#' @param breaks
#' @param labels
#' 
#' @rdname bin_times
#' @author Jim Hughes, \email{jim.hughes@@pfizer.com}
#' 
#' @export

bin_times <- function(data, time_var, breaks, labels) {
  bins <- cut(data[[time_var]], breaks, labels, right = FALSE)
  mutate(data, TIMEBIN = as.numeric(paste(bins)))
}

pred_correction <- function(data = NULL, obs_data = NULL, sim_data = NULL, 
  time_var = "TIME", timebin_var = "TIMEBIN", dv_var = "DV", 
  facet_cov = "PROJ", addl_facet = NULL, 
  id_var = "ID", sim_var = "SIM", pred_var = "PRED", log_dv = FALSE) {
### Perform prediction correct on VPC data  - - - - - - - - - - - - - - - - - - 
# TODO: Add logic to handle log-transformed DV
  if (log_dv) stop("Function does not currently handle log-transformed DV. Please transform prior to using this function.")
# Process input data
  if (!is.null(data)) { 
    if (length(data) != 2) stop("Please provide 2-element list including observed and simulated data when using `data`.")
    obs_input <- data[[1]]
    sim_input <- data[[2]]
  } else if (is.null(data)) { 
    if (is.null(obs_data) | is.null(sim_data)) stop("Please provide both `obs_data` and `sim_data`.")
    obs_input <- obs_data
    sim_input <- sim_data
  }
# Define summary groups and summarise data to obtain median value of PRED
  summaryGroup <- c(facet_cov, timebin_var, addl_facet)
  sim_output <- mutate(sim_input, PREDMED = median(get(pred_var)), .by = summaryGroup)
# Join PRED values with observed data (
# Dev note: Safer than col_bind workflow with small loss of efficiency
  obs_pred <- sim_output %>%
    filter(if_any(any_of(sim_var), \(x) x == 1)) %>%
    distinct(across(all_of(c(id_var, time_var, summaryGroup, pred_var, "PREDMED"))))
  obs_output <- left_join(obs_input, obs_pred, by = c(id_var, time_var, summaryGroup))
# Perform prediction-correction
  list_data <- list(obs = obs_output, sim = sim_output)
  map(list_data, mutate, pcDV = get(dv_var)*PREDMED/get(pred_var))
}

summary_fn <- function(ci) {
  p <- (1- 90/100)/2
  list(
    med = function(x) median(x, na.rm = TRUE), 
    qlo = function(x) quantile(x, probs = p, na.rm = TRUE), 
    qhi = function(x) quantile(x, probs = 1 - p, na.rm = TRUE)
  )
}

summarise_vpc <- function(data, facet_cov = "PROJ", time_var = "TIMEBIN", 
  dv_var = "DV", sim_var = "SIM", addl_facet = NULL, vpc_type = "continuous",
  line_summary = summary_fn(90), ribbon_summary = summary_fn(95), 
  interp_res = 1/60, haz_method = "local", haz_bound = "left") {
### Perform summary for VPC lines and ribbons - - - - - - - - - - - - - - - - - 
# Check whether data is observed or simulated
  if (is.null(data[[sim_var]])) sim_var <- NULL
# For continous and categorical VPCs...
  if (vpc_type %in% c("continuous", "categorical")) {
  # Define summary groups
    summaryGroup <- c(sim_var, facet_cov, time_var, addl_facet)
  # Summarise data based on line_summary function
    line_data <- data %>%
      group_by(across(all_of(summaryGroup))) %>%
      summarise(across(all_of(dv_var), line_summary, .names = "{.fn}"), .groups = "drop")
  # For categorical VPCs...
    if (vpc_type == "categorical") {
    # Pivot categories (summarised by custom line_summary function) into long format
      line_data <- line_data %>% 
        pivot_longer(all_of(names(line_summary)), names_to = "dv_cat", values_to = "med") %>%
        mutate(dv_cat = factor(dv_cat, levels = names(line_summary)))
    }  # if `vpc_type == "categorical"`
# For survival and hazard VPCs...
  } else if (vpc_type %in% c("survival", "hazard")) {
  # Define summary groups (time not included in summary) to nest/group data
    summaryGroup <- c(sim_var, facet_cov, addl_facet)
    nest_data <- nest(data, surv_data = !any_of(summaryGroup))
  # For survival VPCs...
    if (vpc_type == "survival") {
    # Define formula used for survival model
      survFormula <- as.formula(paste0("Surv(", time_var, ", ", dv_var, ") ~ 1"))
    # Calculate survival for each summary group
      line_data <- nest_data %>%
        mutate(surv_fit = map(surv_data, survfit, formula = survFormula)) %>%
        mutate(surv_out = map2(surv_fit, surv_data, surv_summary)) %>%
        select(-surv_data, -surv_fit) %>%
        unnest(surv_out)
    # If data was simulated...
      if (!is.null(sim_var)) {
      # Interpolate survival data across time bins for subsequent summarisation
        line_data <- line_data %>% 
          distinct(across(all_of(summaryGroup))) %>%
          expand_grid(time = seq(0, max(line_data$time), by = interp_res)) %>%
          left_join(line_data, by = c(summaryGroup, "time")) %>%
          group_by(across(all_of(summaryGroup))) %>%
          fill(surv, .direction = "down") %>%
          fill(n.risk, .direction = "up") %>%
          mutate(med = replace_na(surv, 1))
      }  # if `!is.null(sim_var)`
  # For hazard VPCs...
    } else if (vpc_type == "hazard") {
    # Estimate hazard using kernel-based methods
      line_data <- nest_data %>%
        mutate(surv_data = map(surv_data, mutate, hm = haz_method, hb = haz_bound)) %>%
        mutate(haz_fit = map(surv_data, \(x) muhaz(
          time = x[[time_var]], delta = x[[dv_var]], max.time = max(x[[time_var]]),
          bw.method = haz_method, b.cor = haz_bound, 
          n.est.grid = max(x[[time_var]])/interp_res + 1
        ))) %>%
        mutate(time = map(haz_fit, "est.grid")) %>%
        mutate(med = map(haz_fit, "haz.est") )%>%
        select(-surv_data, -haz_fit) %>%
        unnest(c(time, med))
    }  # if `vpc_type == "survival"` else if `vpc_type == "hazard"`
  }  # if `vpc_type %in% c("continuous", "categorical")` else if `vpc_type %in% c("survival", "hazard")`
# Perform summary for VPC ribbons (if simulation data)
  if (is.null(sim_var)) {
    return(line_data)
  } else if (!is.null(sim_var)) {
  # Define additional summary groups needed for simulated data
    sim_facet <- NULL
    if (vpc_type == "categorical") { sim_facet <- "dv_cat" }
    else if (vpc_type %in% c("survival", "hazard")) { sim_facet <- "time" }
  # Define variables to summarise for simulated data
    sim_vars <- "med"
    if (vpc_type == "continuous") { sim_vars <- names(line_summary) }
    else if (vpc_type == "survival") { sim_vars <- c(sim_vars, "n.risk") }
  # Summarise data based on ribbon_summary function
    ribbon_data <- line_data %>%
      group_by(across(all_of(c(summaryGroup[-1], sim_facet)))) %>%
      summarise(across(all_of(sim_vars), ribbon_summary), .groups = "drop")
    return(ribbon_data)
  }  # if `is.null(sim_var)` else if `!is.null(sim_var)`
}

plot_vpc <- function(obs_summary, sim_summary, obs_data = NULL, 
  facet_cov = NULL, vpc_type = "continuous", sim_var = "SIM", 
  time_var = "TIME", timebin_var = "TIMEBIN", dv_var = "DV", 
  time_scale = "linear", dv_scale = "linear", facet_scale = "fixed",
  time_label = waiver(), time_breaks = waiver(), time_brklab = waiver(), 
  dv_breaks = waiver(), dv_label = waiver(), dv_brklab = waiver(), 
  time_limits = c(NA, NA), dv_limits =  c(NA, NA),
  obs_point_colour = "#000067", obs_line_colour = "black", 
  sim_med_colour = "#D55E00", sim_int_colour = "#0095FF",
  sim_alpha = 0.3, line_size = 1, show_legend = TRUE, facet_ncol = NULL, 
  addl_facet = NULL, save_png = TRUE, png_name = "VPC_", png_width = 8.5, 
  png_height = 10, png_units = "in", png_res = 300) {
### Create plot  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# With the deprecation of aes_string in ggplot2 ver. 3.0.0, have opted to use
# the injection operator !! such that strings are converted to object names.
# aes(x = TIME) is equivalent to aes(x = !!sym(time_var)) where time_var == "TIME"
  p <- NULL
  p <- ggplot()
# Define individual observed elements
# Include observed data as points for continuous endpoints...
  if (vpc_type == "continuous") {
    # Include observed data as points and confidence intervals as dashed lines
    p <- p + geom_point(aes(x = !!sym(time_var), y = !!sym(dv_var)), 
      data = obs_data, colour = obs_point_colour, shape = 1)
  }
# Define simulated elements
# Always include median and confidence intervals of the median
  p <- p + geom_ribbon(aes(x = !!sym(timebin_var), ymin = med_qlo, ymax = med_qhi), 
    data = sim_summary, alpha = sim_alpha, fill = sim_med_colour)
  p <- p + geom_line(aes(x = !!sym(timebin_var), y = med_med), 
    data = sim_summary, colour = sim_med_colour, linewidth = line_size)
# For continuous endpoints...
  if (vpc_type == "continuous") {
  # Include median and confidence intervals of the Xth and Yth percentiles
    p <- p + geom_ribbon(aes(x = !!sym(timebin_var), ymin = qlo_qlo, ymax = qlo_qhi), 
      data = sim_summary, alpha = sim_alpha, fill = sim_int_colour)
    p <- p + geom_line(aes(x = !!sym(timebin_var), y = qlo_med), 
      data = sim_summary, colour = sim_int_colour, linetype = "dashed", 
      linewidth = line_size)
    p <- p + geom_ribbon(aes(x = !!sym(timebin_var), ymin = qhi_qlo, ymax = qhi_qhi), 
      data = sim_summary, alpha = sim_alpha, fill = sim_int_colour)
    p <- p + geom_line(aes(x = !!sym(timebin_var), y = qhi_med),  
      data = sim_summary, colour = sim_int_colour, linetype = "dashed", 
      linewidth = line_size)
  }  # if `vpc_type == "continuous"`
# Define summary observed elements
# For endpoints other than survival...
  if (vpc_type != "survival") {
  # Include median line (proportion for categorical, mean hazard for hazard)
    p <- p + geom_line(aes(x = !!sym(timebin_var), y = med), 
      data = obs_summary, colour = obs_line_colour, linewidth = line_size)
  } else if (vpc_type == "survival") {
  # Survival endpoints instead should have a step function
    p <- p + geom_step(aes(x = time, y = surv), 
      data = obs_summary, colour = obs_line_colour, linewidth = line_size)
  }  # if `vpc_type != "survival"` else if `vpc_type == "survival"`
# For continuous endpoints...
  if (vpc_type == "continuous") {
  # Include observed confidence intervals as dashed lines
    p <- p + geom_line(aes(x = !!sym(timebin_var), y = qlo), 
      data = obs_summary, colour = obs_line_colour, linetype = "dashed", 
      linewidth = line_size)
    p <- p + geom_line(aes(x = !!sym(timebin_var), y = qhi),
      data = obs_summary, colour = obs_line_colour, linetype = "dashed",
      linewidth = line_size)
  }  # if `vpc_type == "continuous"`
# Set labels, breaks and scale for axes
  p <- p + labs(x = time_label, y = dv_label)
  if (time_scale == "linear") {
    p <- p + scale_x_continuous(breaks = time_breaks, labels = time_brklab)
  } else if (time_scale == "log") {
    p <- p + scale_x_log10(breaks = time_breaks, labels = time_brklab)
  }  # if `time_scale == "linear"` else if `time_scale == "log"`
  if (dv_scale == "linear") {
    p <- p + scale_y_continuous(breaks = dv_breaks, labels = dv_brklab)
  } else if (dv_scale == "log") {
    p <- p + scale_y_log10(breaks = dv_breaks, labels = dv_brklab)
  }  # if.else
  p <- p + coord_cartesian(xlim = time_limits, ylim = dv_limits)
# Define plot facets
  if (vpc_type == "categorical") { addl_facet <- c(addl_facet, "dv_cat") }
  if (!is.null(facet_cov) & is.null(addl_facet)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_cov)), 
      ncol = facet_ncol, scales = facet_scale)
  } else if (!is.null(facet_cov) & !is.null(addl_facet)) {
    form_tail <- paste(addl_facet, collapse = "+")
    p <- p + facet_wrap(as.formula(paste(facet_cov, "~", form_tail)), 
      ncol = facet_ncol, scales = facet_scale)
  }  # if `(!is.null(facet_cov) & is.null(addl_facet))` else if `!is.null(facet_cov) & !is.null(addl_facet)`
# Save as separate .png files
  if (save_png) {
    plotFile <- paste0(png_name, facet_cov, ".png")
    png(plotFile, width = png_width, height = png_height, units = png_units, 
      res = png_res) 
    print(p)
    dev.off()
    p
  } else { # if `save_png`
    print(p)
  }
}
