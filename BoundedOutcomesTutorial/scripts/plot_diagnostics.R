# Initiate sourcing of script
message("Sourcing plot_diagnostics.R")
walk(c("tidyverse", "grid", "gridExtra", "GGally", "scales"), \(x) {
  if (!(x %in% .packages())) { 
    message("Loading plot_diagnostics dependencies: ", x)
    library(x, character.only = TRUE) 
  }
})

# Function for extracting legend from ggplot2 object
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
  extract_legend <- function(plotobj) {
    tmp <- ggplot_gtable(ggplot_build(plotobj))
    leg <- which(map_chr(tmp$grobs, "name") == "guide-box")
    if (length(leg) > 0) {
      return(tmp$grobs[[leg]])
    } else {
      return(NULL)
    }
  } # g_legend
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Diagnostic plots ------------------------------------------------------------
# Different colours for different values of covariates
# 1. Observed versus individual predicted
# 2. Observed versus population predicted
# 3. Residual variable versus time
# 4. Residual variable versus population predicted
  obs_gof_plot <- function(plot_data, pred_name, dv_name, cov_name, units = NULL, 
    pred_label = "Population Predicted", dv_label = "Observed Value", 
    legend = TRUE, log_scale = FALSE, range_names = NULL,
    log_breaks = c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000),
    palette = c("#000080", "#6343AB", "#008AEC", "#21B5A1", "#9D1B59")) {
    p <- NULL
    p <- ggplot(plot_data)
  # Line of identity
    p <- p + geom_abline(aes(intercept = 0, slope = 1), linewidth = 1)
  # Points and regression line (ipred vs. observed)
    p <- p + geom_point(aes(x = !!sym(pred_name), y = !!sym(dv_name), 
      colour = !!sym(cov_name)), shape = if_else(is.null(range_names), 1, 19), 
      size = 2, alpha = 0.7)
    if (is.null(range_names)) {
      p <- p + geom_smooth(aes(x = !!sym(pred_name), y = !!sym(dv_name), 
        colour = !!sym(cov_name)), linewidth = 1, method = lm, formula = "y ~ x", se = F)
    } else {
      p <- p + geom_linerange(aes(xmin = !!sym(range_names[[1]]), 
        xmax = !!sym(range_names[[2]]), y = !!sym(dv_name), colour = !!sym(cov_name)))
    }
  # Axes
    x_label <- pred_label
    y_label <- dv_label
    if (!is.null(units)) {
      x_label <- paste0(pred_label, " (", units, ")")
      y_label <- paste0(dv_label, " (", units, ")")
    }
    if (isTRUE(log_scale)) {
      p <- p + scale_x_log10(x_label, breaks = log_breaks, labels = log_breaks)
      p <- p + scale_y_log10(y_label, breaks = log_breaks, labels = log_breaks)
    } else {
      p <- p + scale_x_continuous(x_label, labels = comma)
      p <- p + scale_y_continuous(y_label, labels = comma)
    } # if/else
  # Legend
  # Only use custom colour palette if length of unique covariate values
  # is less than or equal to the length of unique colours provided
    if (length(unique(plot_data[[cov_name]])) <= length(unique(palette))) {
      p <- p + scale_colour_manual(cov_name, values = palette)
    } # if/else
  # Define legend if requested
    if (isTRUE(legend)) {
      p <- p + theme(legend.position = "bottom")
    } else {
      p <- p + theme(legend.position = "none")
    } # if/else
    p
  }
  
  res_gof_plot <- function(plot_data, res_name, x_name, cov_name, units = NULL, 
    x_label, res_label = "Conditional Weighted Residual", range_names = NULL,
    res_limit = 6, legend = TRUE, palette = c("#000080", "#6343AB", "#008AEC", "#21B5A1", "#9D1B59")) {
    p <- NULL
    p <- ggplot(plot_data)
  # Points of res variable versus time
    p <- p + geom_point(aes(x = !!sym(x_name), y = !!sym(res_name), 
      colour = !!sym(cov_name)), shape = if_else(is.null(range_names), 1, 19), size = 2)
    if (!is.null(range_names)) {
      p <- p + geom_linerange(aes(x = !!sym(x_name), ymin = !!sym(range_names[[1]]), 
        ymax = !!sym(range_names[[2]]), colour = !!sym(cov_name)))
    }
  # Lines for +/- 2 and residual limit as defined by res_limits
    p <- p + geom_abline(aes(intercept = 0, slope = 0))
    p <- p + geom_abline(aes(intercept = -2, slope = 0), linetype = "dashed")
    p <- p + geom_abline(aes(intercept = 2, slope = 0), linetype = "dashed")
    p <- p + geom_abline(aes(intercept = -res_limit, slope = 0), 
      linetype = "dashed", linewidth = 1)
    p <- p + geom_abline(aes(intercept = res_limit, slope = 0), 
      linetype = "dashed", linewidth = 1)
  # Regression line
    if (is.null(range_names)) {
      p <- p + geom_smooth(aes(x = !!sym(x_name), y = !!sym(res_name), 
        colour = !!sym(cov_name)), method = lm, formula = "y ~ x", se = F)
    }
  # Axes
    xax_label <- x_label
    if (!is.null(units)) xax_label <- paste0(x_label, " (", units, ")")
    p <- p + scale_y_continuous(res_label, labels = comma)
    p <- p + scale_x_continuous(xax_label, labels = comma)
  # Legend
  # Only use custom colour palette if length of unique covariate values
  # is less than or equal to the length of unique colours provided
    if (length(unique(plot_data[[cov_name]])) <= length(unique(palette))) {
      p <- p + scale_colour_manual(cov_name, values = palette)
    } # if/else
    if (isTRUE(legend)) {
      p <- p + theme(legend.position = "bottom")
  # Define legend if requested
    } else {
      p <- p + theme(legend.position = "none")
    } # if/else
    p
  }
  
  combine_gof_plot <- function(...) {
    all_plots <- list(...)
    nplots <- length(all_plots)
    # Define parameters for creation of custom layout_matrix for grid.arrange
    if (nplots == 4) {
      row_layout <- list(
        top = c(1, 1, 2, 2),
        mid = c(3, 3, 3, 3),
        bot = c(4, 4, 4, 4),
        leg = c(5, 5, 5, 5)
      )
      row_height <- list(top = 8, mid = 4, bot = 4, leg = 1)
    } else if (nplots == 5) {
      row_layout <- list(
        top = c(1, 1, 2, 2),
        tmd = c(3, 3, 3, 3),
        bmd = c(4, 4, 4, 4),
        bot = c(5, 5, 5, 5),
        leg = c(6, 6, 6, 6)
      )
      row_height <- list(top = 7, tmd = 4, bmd = 4, bot = 4, leg = 1)
    } else {
      stop("`combine_gof_plot` only supports 4-5 diagnostics plots.\n")
    }
    # If using "PROJ" covariate, remove legend from layout
    plot_legend <- extract_legend(all_plots[[1]])
    if (is.null(plot_legend)) {
      row_layout[["leg"]] <- NULL
      row_height[["leg"]] <- NULL
      plot_legend <- NULL
      diag_index <- seq_len(nplots)
      # Otherwise, extract then remove legend from ipredVsObs plot
    } else {
      plot_legend <- extract_legend(all_plots[[1]])
      diag_index <- seq_len(nplots + 1)
      all_plots <- map(all_plots, \(p) p + theme(legend.position = "none"))
    }  # if/else
    # Define layout 
    diag_layout <- row_layout %>%
      map2(row_height, rep.int) %>%
      map(matrix, nrow = 4) %>%
      map(t) %>%
      do.call(what = rbind)
    # Arrange plots using custom layout_matrix
    # Loop over for the number of res variables in resNames
      diag_plots <- unlist(list(all_plots, list(plot_legend)), recursive = FALSE)
      marrangeGrob(grobs = diag_plots[diag_index], layout_matrix = diag_layout, top = "")
  }
  
  hist_plot <- function(plot_data, sample_names, sample_label = NULL, palette = "#0093D0") {
    map(sample_names, function(sample_name) {
      # Create a bin.width appropriate for the data
      sample_values <- plot_data[[sample_name]]
      nval <- length(sample_values) # Number of ETA values in dataset
      nbins <- ceiling(sqrt(nval)) # Square of nval rounded up
      bin_width <- (max(sample_values) - min(sample_values))/nbins
      # Create plot including histogram and density line
      p <- NULL
      p <- ggplot()
      p <- p + geom_histogram(aes(x = sample_values, y = after_stat(density)), 
        fill = palette[1], colour = palette[1], alpha = 0.7, 
        binwidth = bin_width)
      p <- p + geom_density(aes(x = sample_values, y = after_stat(density)), size = 1)
      p <- p + geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1)
      p <- p + scale_y_continuous("Distribution Density")
      if (is.null(sample_label)) sample_label <- sample_name
      p <- p + scale_x_continuous(sample_label)
      return(p)
    }) # map
  }
  
  combine_dist_plot <- function(hist_plots, qq_plots, plots_per_page = 6) {
  # Arrange histograms and Q-Q plots into a grid and save
    dist_plots <- unlist(map2(hist_plots, qq_plots, list), recursive = FALSE)
  # Calculate number of plots to be arranged
    nplots <- length(dist_plots)
  # Number of pages required to plot all plots with nplots to a page
    npages <- ceiling(nplots/plots_per_page)
  # Assign a page number to each plot
    page_index <- seq_len(npages) %>%
      rep(times = plots_per_page) %>%
      sort()
    page_data <- tibble(plot_index = seq_len(nplots), page_index = page_index[plot_index])
  # Arrange and save
    plot_list <- map(seq_len(npages), function(page) {
      plot_index <- page_data$plot_index[page_data$page_index == page]
      comb_plot <- grid.arrange(grobs = dist_plots[plot_index], 
        ncol = 2, nrow = plots_per_page/2)
      print(comb_plot)
    }) # map
  }
  
  eta_cov_plot <- function(subject_level_data, eta_names, contcov_names, catcov_names, 
    palette = "#0093D0", char_per_level = 40, plots_per_page = 6) {
  # TODO: Currently prints to console, want to avoid that if possible
  # TODO: Assumes categorical covariates are factors.. could be more robust
    map(eta_names, function(eta) {
    # Plot eta versus continuous covariates
      eta_cont_plots <- map(contcov_names, function(cov) {
        plot_data <- subject_level_data %>%
          mutate(across(all_of(cov), \(x) na_if(x, -999)))
      # Define scatter plot with loess smooth
        p <- NULL
        p <- ggplot(data = plot_data)
        p <- p + geom_point(aes(x = !!sym(cov), y = !!sym(eta)), 
          colour = palette[1], size = 2, shape = 1)
        p <- p + geom_smooth(aes(x = !!sym(cov), y = !!sym(eta)), 
          colour = "black", size = 1, method = lm, formula = "y ~ x", se = T)
        p <- p + geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1)
        p <- p + scale_y_continuous(eta)
        p <- p + scale_x_continuous(cov)
      # Return plot object
        eta_cont_plot <- p
      }) # map
    # Plot eta versus categorical covariates
      eta_cat_plots <- map(catcov_names, function(cov) {
      # Determine number of categories and longest category by letters
        max_char <- max(nchar(levels(subject_level_data[[cov]])))
        nlevels <- length(levels(subject_level_data[[cov]]))
      # Define box plot
        p <- NULL
        p <- ggplot(data = subject_level_data)
        p <- p + geom_boxplot(aes(x = !!sym(cov), y = !!sym(eta)), 
          fill = palette[1], alpha = 0.7)
        p <- p + geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1)
        p <- p + scale_y_continuous(eta)
        if (max_char < char_per_level/nlevels) {
          p <- p + scale_x_discrete(cov)
        } else {
          p <- p + scale_x_discrete("")
          p <- p + theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust = 1.0))
        }  # if/else
      # Return plot object
        eta_cat_plot <- p
      })  # map
    # Arrange all plots into a grid
      plot_list <- unlist(list(eta_cont_plots, eta_cat_plots), recursive = FALSE)
    # Calculate number of plots to be arranged
      nplots <- length(plot_list)
    # Number of pages required to plot all plots with nplots to a page
      npages <- ceiling(nplots/plots_per_page)
    # Assign a page number to each plot
      page_index <- seq_len(npages) %>%
        rep(times = plots_per_page) %>%
        sort()
      page_data <- tibble(plot_index = seq_len(nplots), page_index = page_index[plot_index])
      eta_cov_plots <- map(seq_len(npages), function(page) {
        plot_index <- page_data$plot_index[page_data$page_index == page]
        comb_plot <- grid.arrange(grobs = plot_list[plot_index],
          ncol = 2, nrow = plots_per_page/2)
      })  # map
    })  # map
  }
  
  qq_plot <- function(plot_data, sample_names, use_dist = FALSE, palette = "#0093D0") {
    map(sample_names, function(sample_name) {
    # Extract sample values from data
      sample_values <- plot_data[[sample_name]]
    # Define theoretical normal distribution
      theor_dist <- list(mean = 0, sd = 1)
      if (use_dist) theor_dist <- list(mean = mean(sample_values), sd = sd(sample_values))
    # Q-Q plot of ETAs
      p <- NULL
      p <- ggplot(plot_data)
      p <- p + geom_abline(aes(intercept = 0, slope = 1), size = 1)
      p <- p + geom_qq(aes(sample = !!sym(sample_name)), dparams = theor_dist, 
        geom = "point", shape = 1, colour = palette[1], size = 2)
      p <- p + scale_y_continuous("Sample Quantile")
      p <- p + scale_x_continuous("Theoretical Quantile")
      return(p)
    }) # map
  }
  
  corr_plot <- function(subject_level_data, column_names, palette = "#0093D0") {
    ggpairs(subject_level_data,
      columns = column_names,
      upper = list(
        continuous = wrap("cor", colour = palette[1]),
        combo = wrap("box", fill = palette[1], alpha = 0.7),
        discrete = wrap("ratio", colour = palette[1])
      ),  # upper
      diag = list(
        continuous = wrap("densityDiag", fill = palette[1], alpha = 0.7),
        discrete = wrap("barDiag", fill = palette[1], alpha = 0.7)
      ),  # diag
      lower = list(
        continuous = wrap("smooth", colour = palette[1], shape = 1),
        combo = wrap("dot", colour = palette[1], shape = 1),
        discrete = wrap("facetbar", fill = palette[1])
      )  # lower
    )  # ggpairs
  }
  
  collin_plot <- function(cov_matrix, palette = "Spectral") {
  # Check whether covariance step was successful
    if (!is.null(cov_matrix)) {
    # Determine font sizes based on number of parameters in plot
      param_labels <- colnames(cov_matrix)
      nlabels <- length(param_labels)
      if (nlabels < 10) plot_font_size <- 5
      if (nlabels >= 10 & nlabels < 20) plot_font_size <- 4
      if (nlabels >= 20 & nlabels < 30) plot_font_size <- 3
      if (nlabels >= 30) plot_font_size <- 2
    # Plot
      p <- ggcorr(data = NULL, cor_matrix = cov_matrix,
        palette = palette,
        label = TRUE,
        label_color = "black",
        label_size = plot_font_size,
        label_round = 2,
        nbreaks = 10,
        hjust = 0.9,
        layout.exp = ceiling(length(param_labels)/20),
        legend.position = "none",
        check_overlap = TRUE,
        size = plot_font_size
      )  # ggcorr
      return(p)
    } else {
      message("Cannot find any files with .cov extension.", 
        "Will not generate plot to explore parameter colinearity.")
      return(NULL)
    }
  }
  
  individual_plots <- function(id_data, dv_name, ipred_name, pred_name, 
    time_name, id_name, flag_name, plots_per_page = 12, 
    x_label = NULL, y_label = NULL, x_limits = NULL, y_limits = NULL, 
    facet_nrow = NULL, facet_ncol = 3, 
    facet_scales = "free",
    palette = c("#000080", "#6343AB", "#008AEC", "#21B5A1", "#9D1B59")) {
    unique_id <- unique(id_data[[id_name]])
    unique_flag <- unique(id_data[[flag_name]])
    split_data <- id_data %>%
      mutate(facet_var = paste(get(id_name), get(flag_name), sep = "\nFLAG: ")) %>%
      mutate(plot_index = as.integer(factor(facet_var, levels = unique(facet_var)))) %>%
      mutate(page_index = ceiling(plot_index / plots_per_page)) %>%
      group_split(page_index)
    id_plots <- map(split_data, \(plot_data) {
      pal_index <- which(unique_flag %in% unique(plot_data[[flag_name]]))
      p <- ggplot(data = plot_data)
      p <- p + geom_point(aes(x = !!sym(time_name), y = !!sym(dv_name), colour = !!sym(flag_name)))
      p <- p + geom_point(aes(x = !!sym(time_name), y = !!sym(ipred_name), shape = "Individual"))
      p <- p + geom_point(aes(x = !!sym(time_name), y = !!sym(pred_name), shape = "Population"))
      if (nrow(plot_data) > 1) {
        p <- p + geom_line(aes(x = !!sym(time_name), y = !!sym(ipred_name), linetype = "Individual"))
        p <- p + geom_line(aes(x = !!sym(time_name), y = !!sym(pred_name), linetype = "Population"))
      }
      p <- p + labs(x = x_label, y = y_label)
      p <- p + coord_cartesian(xlim = x_limits, ylim = y_limits)
      p <- p + facet_wrap(~facet_var, scales = "free", nrow = facet_nrow, ncol = facet_ncol)
      p <- p + scale_colour_manual("Endpoint", values = palette[pal_index])
      p <- p + scale_shape_manual("Prediction", values = c(1, 2))
      p <- p + scale_linetype_manual("Prediction", values = c("solid", "dashed"))
      p <- p + theme(legend.position = "none", plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"),
        axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10))
    })
    plot_grid <- marrangeGrob(grobs = id_plots, ncol = 1, nrow = 1,
      npages = ceiling(length(id_plots)/plots_per_page))
  }

message("The following functions have been loaded:
 - obs_gof_plot
 - res_gof_plot
 - eta_cov_plot
 - combine_gof_plot
 - corr_plot
 - qq_plot
 - combine_eta_plot
 - collin_plot
")

  