library(sensobol)
library(ggplot2)
library(cowplot)
library(patchwork)
library(GGally)

source("scripts/sensitivity/functions.R")

theme_set(theme_cowplot())

sobol_neg_and_pos <- function() {
  # Functions
  # ---------
  load_trace_data <- function(pos = FALSE) {
    # Set the (global, hence <<-) pe variable, and the data dir, to the
    # appropriate value
    if (pos) {
      pe <<- 5.0
      res_dir_base <- "res_sensitivity_1000_pe_5_inhib_trace_only"
    } else {
      pe <<- -5.0
      res_dir_base <- "res_sensitivity_1000_pe_-5_inhib_trace_only"
    }

    # The names of the params that are being varied
    param_names <- c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag")

    # Read the parameters from file
    x <- readRDS(paste(res_dir_base, "param_sample.rds", sep = "/"))
    n_param_sample <- x$n_param_sample
    param_sample <- x$param_sample
    param_sample_dt <- x$param_sample_dt

    trace_data_full <- read_trace_data(param_sample_dt, res_dir_base)

    # Calculate the total cell outflux as the sum of the fluxes of the three different cell types
    trace_data_full[, cell_outflux := `-F_{phi_{C_u}}(x=0)` + `-F_{phi_{C_b}}(x=0)` + `-F_{phi_{C_s}}(x=0)`]

    trace_data <- trace_data_full[`t_{inf}` >= 0]

    list(
      n_param_sample = n_param_sample,
      param_sample = param_sample,
      param_sample_dt = param_sample_dt,
      trace_data = trace_data,
      param_names = param_names
    )
  }

  sobol_indices_cells_panel <- function(
    ind_neg,
    ind_pos
  ) {
    ind_neg$results[, pe := "neg"]
    ind_pos$results[, pe := "pos"]
    data <- rbind(ind_neg$results, ind_pos$results)
    data[, factor(
      parameters,
      levels = c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag")
    )]

    x_label <- "Sobol index"

    p <- ggplot(
      data,
      aes(
        x = original,
        y = parameters,
        xmin = low.ci,
        xmax = high.ci,
        shape = sensitivity,
        colour = pe
      )
    ) +
      geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
      geom_errorbar(width = 0.5, position = position_dodge(width = 0.6)) +
      scale_y_discrete(labels = param_labels_words_no_breaks, limits = rev) +
      scale_shape_discrete(
        labels = c("Si" = "First-order", "Ti" = "Total-order")
      ) +
      scale_colour_discrete(
        labels = c("pos" = expression(Pe == 5), "neg" = expression(Pe == -5))
      ) +
      coord_cartesian(xlim = c(0.0, 1.2)) +
      labs(
        x = x_label,
        y = NULL,
        shape = "Sobol index",
        colour = "Flow"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.text.align = 0 #,
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
      )

    p
  }

  cells_vs_params_long <- function(param_sample_dt, integrated_fluxes) {
    cells_vs_params <- param_sample_dt[
      ,
      .(j_phi_i_i_factor,
        m_i_factor,
        t_j_phi_i_lag,
        cells_in = integrated_fluxes[, cells_in],
        cells_out = integrated_fluxes[, cells_out],
        net_change = integrated_fluxes[, net_change]
      )
    ]

    cells_vs_params %>%
      melt(
        .,
        measure.vars = c(
          "j_phi_i_i_factor",
          "m_i_factor",
          "t_j_phi_i_lag"
        ),
        variable.name = "param",
        value.name = "param_value"
      ) %>%
      melt(
        .,
        measure.vars = c("cells_in", "cells_out", "net_change"),
        variable.name = "variable",
        value.name = "cells"
      )
  }

  cells_vs_params_panel <- function(cells_vs_params, variable_str, param_str) {
    p <- ggplot(cells_vs_params[variable == variable_str & param == param_str]) +
      geom_point(aes(
        x = param_value,
        y = cells,
        colour = factor(pe),
        group = factor(pe),
        shape = factor(pe)
      ), size = 0.2, alpha = 0.5) +
      labs(
        x = param_labels_words_no_breaks[param_str],
        y = variable_labels[variable_str],
        colour = NULL,
        shape = NULL
      ) +
      scale_colour_discrete(
        labels = c(
          "5" = expression(Pe == 5),
          "-5" = expression(Pe == -5)
        )
      ) +
      scale_shape_discrete(
        labels = c(
          "5" = expression(Pe == 5),
          "-5" = expression(Pe == -5)
        )
      ) +
      theme(
        #strip.background = element_blank(),
        #strip.placement = "outside",
        #legend.text.align = 0,
        # TODO re-enable the legend but find a better place to fit it in - it's
        # on the RHS at the moment and squeezes the plots a lot
        legend.position = "none"
      )

      if (variable_str == "net_change") {
        p <- p + geom_vline(xintercept = 0, linetype = "dashed")
      }

      p
  }

  # Plots
  # -----

  # Load trace data for neg Pe
  data_neg <- load_trace_data(pos = FALSE)
  integrated_fluxes_neg <- calculate_integrated_fluxes(data_neg$trace_data)
  integrated_fluxes_neg[, net_change := cells_in - cells_out]

  # Load trace data for pos Pe
  data_pos <- load_trace_data(pos = TRUE)
  integrated_fluxes_pos <- calculate_integrated_fluxes(data_pos$trace_data)
  integrated_fluxes_pos[, net_change := cells_in - cells_out]

  # Cells in
  ind_neg <- sobol_indices(
    Y = integrated_fluxes_neg$cells_in,
    N = data_neg$n_param_sample,
    params = data_neg$param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  ind_pos <- sobol_indices(
    Y = integrated_fluxes_pos$cells_in,
    N = data_pos$n_param_sample,
    params = data_pos$param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  p_sobol_indices_cells_in <- sobol_indices_cells_panel(
    ind_neg, ind_pos
  )

  # Cells out
  ind_neg <- sobol_indices(
    Y = integrated_fluxes_neg$cells_out,
    N = data_neg$n_param_sample,
    params = data_neg$param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  ind_pos <- sobol_indices(
    Y = integrated_fluxes_pos$cells_out,
    N = data_pos$n_param_sample,
    params = data_pos$param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  p_sobol_indices_cells_out <- sobol_indices_cells_panel(
    ind_neg, ind_pos
  )

  # Net change
  ind_neg <- sobol_indices(
    Y = integrated_fluxes_neg$net_change,
    N = data_neg$n_param_sample,
    params = data_neg$param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  ind_pos <- sobol_indices(
    Y = integrated_fluxes_pos$net_change,
    N = data_pos$n_param_sample,
    params = data_pos$param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  p_sobol_indices_net_change <- sobol_indices_cells_panel(
    ind_neg, ind_pos
  )

  # Cell vs params
  cells_vs_params_long_neg <- cells_vs_params_long(
    data_neg$param_sample_dt,
    integrated_fluxes_neg
  )
  cells_vs_params_long_pos <- cells_vs_params_long(
    data_pos$param_sample_dt,
    integrated_fluxes_pos
  )

  cells_vs_params_long_neg[, pe := -5]
  cells_vs_params_long_pos[, pe := 5]

  cells_vs_params_long_all <- rbind(
    cells_vs_params_long_neg,
    cells_vs_params_long_pos
  )

  p_cells_in_panels <- list()
  p_cells_out_panels <- list()
  p_net_change_panels <- list()

  param_strs <- cells_vs_params_long_all[, levels(param)]

  for (i in 1:length(param_strs)) {
    p_cells_in_panels[[i]] <- cells_vs_params_panel(
      cells_vs_params_long_all,
      "cells_in",
      param_strs[i]
    )

    p_cells_out_panels[[i]] <- cells_vs_params_panel(
      cells_vs_params_long_all,
      "cells_out",
      param_strs[i]
    )

    p_net_change_panels[[i]] <- cells_vs_params_panel(
      cells_vs_params_long_all,
      "net_change",
      param_strs[i]
    )
  }

  list(
    p_cells_in_panels,
    p_cells_out_panels,
    p_net_change_panels,
    p_sobol_indices_cells_in,
    p_sobol_indices_cells_out,
    p_sobol_indices_net_change
  )
}

plot_dir <- "plots_neg_and_pos_inhib"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

plots <- sobol_neg_and_pos()

p_cells_in_panels <- plots[[1]]
p_cells_out_panels <- plots[[2]]
p_net_change_panels <- plots[[3]]
p_sobol_indices_cells_in <- plots[[4]]
p_sobol_indices_cells_out <- plots[[5]]
p_sobol_indices_net_change <- plots[[6]]

p_sobol_and_cells_combined <- p_sobol_indices_cells_in +
  labs(tag = "A", title = "Cells in") +
  plot_spacer() +
  p_cells_in_panels[[1]] +
  p_cells_in_panels[[2]] +
  p_cells_in_panels[[3]] +
  p_sobol_indices_cells_out +
  labs(tag = "B", title = "Cells out") +
  plot_spacer() +
  p_cells_out_panels[[1]] +
  p_cells_out_panels[[2]] +
  p_cells_out_panels[[3]] +
  p_sobol_indices_net_change +
  labs(tag = "C", title = "Net change") +
  plot_spacer() +
  p_net_change_panels[[1]] +
  p_net_change_panels[[2]] +
  p_net_change_panels[[3]] +
  plot_layout(ncol = 5, nrow = 3, widths = c(1, 0.2, 1, 1, 1), guides = "collect", axis_title = "collect")

ggsave_with_defaults(
  plot = p_sobol_and_cells_combined,
  paste(plot_dir, "sobol_and_cells_combined.pdf", sep = "/"),
  width = 14,
  height = 7
)
