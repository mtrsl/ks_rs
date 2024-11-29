library(ggplot2)
library(cowplot)
library(patchwork)
library(ggh4x)

source("scripts/sensitivity/functions.R")

theme_set(theme_cowplot() + background_grid())

res_dir_base <- "res_special_inhib"
plot_dir_base <- paste(res_dir_base, "plots", sep = "/")

if (!dir.exists(plot_dir_base)) {
  dir.create(plot_dir_base, recursive = TRUE)
}

x <- readRDS(paste(res_dir_base, "param_sample.rds", sep = "/"))
params <- x$param_sample_dt

read_spatial_data_at_times <- function(times) {
  data <- lapply(
    times, function(t_inf) read_spatial_data(params, t_inf, res_dir_base)
  ) |> rbindlist()

  # Round the time to a few decimal places as a hack to ensure that all times
  # we consider "equal" when plotting (e.g. 30.000000 and 30.000020) are
  # treated as such by facet_grid etc.
  data[, time_inf := round(time_inf, 3)]
  data
}

read_spatial_data_at_times_long <- function(times) {
  read_spatial_data_at_times(times) |>
    melt(measure.vars = c("C_u", "C_b", "C_s", "phi_i", "phi_m", "phi_C_u", "phi_C_b", "phi_C_s", "J"))
}

custom_scientific <- function(x) {
  sapply(x, function(value) {
    # Check if scientific notation would be used by default
    if (grepl("e", format(value))) {
      # Remove leading zeros in exponent
      gsub("e\\+0?", "e", gsub("e-0?", "e-", format(value, scientific = TRUE)))
    } else {
      format(value, scientific = FALSE)
    }
  })
}

# Combined spatial plots
spatial_plot_subset_combined <- function(plot_times, pes, lag, c_var, phi_c_var) {
  data_subset_long <- read_spatial_data_at_times_long(plot_times, c_var, phi_c_var)

  # Manually set the factor labels so they appear correctly in the facet strips
  # without having to use a complex labeller
  data_subset_long[, j_phi_i_i_factor := factor(
    j_phi_i_i_factor,
    labels = c(
      "Ingress~ratio == 2",
      "Ingress~ratio == 1000"
    )
  )]

  data_subset_long[, m_i_factor := factor(
    m_i_factor,
    labels = c(
      "Maturation~ratio == 2",
      "Maturation~ratio == 1000"
    )
  )]

  data_subset_long[, time_inf := factor(
    time_inf,
    labels = sapply(plot_times, function(t) paste0("t[inf] == ", t / 10))
  )]

  data_subset_long[, variable := factor(
    variable,
    labels = c(
      variable_labels[c_var],
      variable_labels[phi_c_var]
    )
  )]

  ggplot(
    data_subset_long[t_j_phi_i_lag == lag & pe %in% pes],
    aes(
      x = x,
      y = value,
      group = rep,
      colour = factor(pe),
      linetype = factor(gamma)
    )
  ) +
    facet_nested(
      rows = vars(j_phi_i_i_factor, time_inf),
      cols = vars(m_i_factor, variable),
      scales = "free_y",
      independent = "y",
      labeller = label_parsed,
      switch = "y"
    ) +
    geom_line(linewidth = 1) +
    scale_linetype_discrete(limits = rev) +
    theme(
      strip.placement = "outside",
      strip.text = element_text(size = rel(1.25)),
      strip.background = element_rect(fill = "#eaeaea")
    ) +
    labs(
      x = expression(x),
      y = NULL,
      colour = param_labels_words["pe"],
      linetype = param_labels_words["gamma"]
    )
}

# Separate spatial plots
spatial_plot_subset <- function(variables, plot_time, pes, lag, j_phi_i_i_factor_val, m_i_factor_val) {
  data_subset <- read_spatial_data_at_times(plot_time) |>
    _[, .(x, rep, C_u, C_b, C_s, phi_C_u, phi_C_b, phi_C_s, j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma, pe)] |>
    _[pe %in% pes &
      t_j_phi_i_lag == lag &
      j_phi_i_i_factor == j_phi_i_i_factor_val &
      m_i_factor == m_i_factor_val
    ]

  plots <- list()

  for (i in 1:length(variables)) {
    variable <- variables[i]
    plots[[i]] <- ggplot(
      data_subset,
      aes(
        x = x,
        y = .data[[variable]],
        colour = factor(pe),
        linetype = factor(gamma),
        group = rep,
      )
    ) +
      geom_line(linewidth = 1) +
      scale_linetype_discrete(limits = rev, labels = function(x) sprintf("%.2f", as.numeric(x))) +
      scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1.0), labels = c(0, 0.25, 0.5, 0.75, 1)) +
      scale_y_continuous(labels = custom_scientific) +
      guides(color = guide_legend(nrow = 1)) +
      labs(
        title = all_labels[variable],
        x = "Spatial distance",
        y = NULL,
        #y = all_labels[variable],
        colour = param_labels_words_no_breaks["pe"],
        linetype = param_labels_words_no_breaks["gamma"]
      ) +
      theme(
        legend.direction = "horizontal",
      )
  }

  plots
}

# Separate spatial plots as fold change vs homeostasis
spatial_plot_subset_fold_change <- function(variables, plot_time, pes, lag, j_phi_i_i_factor_val, m_i_factor_val) {
  data_subset_hom <- read_spatial_data_at_times(0) |>
    _[, .(x, time_inf, rep, C_u, C_b, C_s, phi_C_u, phi_C_b, phi_C_s, j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma, pe)] |>
    _[pe %in% pes &
      t_j_phi_i_lag == lag &
      j_phi_i_i_factor == j_phi_i_i_factor_val &
      m_i_factor == m_i_factor_val
    ]

  data_subset_at_time <- read_spatial_data_at_times(plot_time) |>
    _[, .(x, time_inf, rep, C_u, C_b, C_s, phi_C_u, phi_C_b, phi_C_s, j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma, pe)] |>
    _[pe %in% pes &
      t_j_phi_i_lag == lag &
      j_phi_i_i_factor == j_phi_i_i_factor_val &
      m_i_factor == m_i_factor_val
    ]

  data_subset_wide <- rbindlist(list(data_subset_hom, data_subset_at_time)) |>
    _[, id := seq_len(.N), by = .(rep, time_inf)] |>
    dcast(... ~ time_inf, value.var = all_vars)

  plots <- list()

  for (i in 1:length(variables)) {
    variable <- variables[i]
    variable_hom <- paste0(variable, "_0")
    # divide by 10 here since we output every 0.1. bit of a hack
    # plot_time isn't actually a time - it's the number of the time step and
    # hence what appears in the name of the output file
    variable_at_time <- paste0(variable, "_", plot_time / 10)
    plots[[i]] <- ggplot(
      data_subset_wide,
      aes(
        x = x,
        y = .data[[variable_at_time]] / .data[[variable_hom]],
        colour = factor(pe),
        linetype = factor(gamma),
        group = rep,
      )
    ) +
      geom_line(linewidth = 1) +
      scale_linetype_discrete(limits = rev, labels = function(x) sprintf("%.2f", as.numeric(x))) +
      scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1.0), labels = c(0, 0.25, 0.5, 0.75, 1)) +
      scale_y_continuous(labels = custom_scientific) +
      guides(color = guide_legend(nrow = 1)) +
      labs(
        title = all_labels[variable],
        x = "Spatial distance",
        y = NULL,
        #y = all_labels[variable],
        colour = param_labels_words_no_breaks["pe"],
        linetype = param_labels_words_no_breaks["gamma"]
      ) +
      theme(
        legend.direction = "horizontal",
      )
  }

  plots
}

#pes <- c(-5, -3, -1, 1, 3, 5)
pes <- c(-5, 0, 5)

# Spatial plots for a range of time points

  #c_var = "C_u"
  #phi_c_var = "phi_C_u"

  #filename_suffix = paste0(c_var, "_", phi_c_var)

  ## No lag
  #plot_times <- c(0, 50, 100, 150, 250, 500)
  #p_spatial_lag_0 <- spatial_plot_subset_combined(plot_times, pes, 0, c_var, phi_c_var)

  #ggsave_with_defaults(
    #plot = p_spatial_lag_0,
    #paste(plot_dir_base, "spatial_lag_0.pdf", sep = "/"),
    #width = 14,
    #height = 20
  #)

  ## Lag 25
  #plot_times <- c(0, 250, 300, 350, 400, 500)
  #p_spatial_lag_25 <- spatial_plot_subset_combined(plot_times, pes, 25, c_var, phi_c_var)

  #ggsave_with_defaults(
    #plot = p_spatial_lag_25,
    #paste(plot_dir_base, "spatial_lag_25.pdf", sep = "/"),
    #width = 14,
    #height = 20
  #)

# Spatial plots for just homeostasis and one other time point

all_vars <- c("C_u", "C_b", "C_s", "phi_C_u", "phi_C_b", "phi_C_s")

time_label <- function(text, tag) {
  grid::textGrob(
    bquote(bold(.(tag)~~~.(text))),
    x = 0,
    y = 0.4,
    just = "left",
    gp = grid::gpar(fontsize = 14)
  )
}

group_label <- function(m_i_factor_val, j_phi_i_i_factor_val) {
  m_i_factor_val = format(m_i_factor_val, scientific = FALSE)
  j_phi_i_i_factor_val = format(j_phi_i_i_factor_val, scientific = FALSE)
  grid::textGrob(
    bquote(underline(list(
      Maturation~ratio == .(m_i_factor_val),
      Ingress~ratio == .(j_phi_i_i_factor_val)
    ))),
    x = 0.5,
    y = 0.3,
    just = "centre"
  )
}

options(scipen = 0)
p_hom <- spatial_plot_subset(all_vars, 0, pes, 0, 2, 2)
p_2_2 <- spatial_plot_subset(all_vars, 500, pes, 0, 2, 2)
p_2_1000 <- spatial_plot_subset(all_vars, 500, pes, 0, 2, 1000)
p_1000_2 <- spatial_plot_subset(all_vars, 500, pes, 0, 1000, 2)
p_1000_1000 <- spatial_plot_subset(all_vars, 500, pes, 0, 1000, 1000)
#p_2_2 <- spatial_plot_subset_fold_change(all_vars, 500, pes, 0, 2, 2)
#p_2_1000 <- spatial_plot_subset_fold_change(all_vars, 500, pes, 0, 2, 1000)
#p_1000_2 <- spatial_plot_subset_fold_change(all_vars, 500, pes, 0, 1000, 2)
#p_1000_1000 <- spatial_plot_subset_fold_change(all_vars, 500, pes, 0, 1000, 1000)

options(scipen = -2)
# Assemble the figure using patchwork
p_spatial_hom_vs_end_of_inf <-
  wrap_elements(full = time_label("Homeostatic steady state", "A")) +
    p_hom[[1]] + p_hom[[2]] + p_hom[[3]] +
    p_hom[[4]] + p_hom[[5]] + p_hom[[6]] +
    guide_area() +
    # gap
    plot_spacer() +
    # end of inflammation
    # j_phi_i_i_factor = 2, m_i_factor = 2; j_phi_i_i_factor = 2, m_i_factor = 1000
    wrap_elements(full = time_label("End of inflammation", "B")) +
    group_label(2, 2) + group_label(1000, 2) +
    p_2_2[[1]] + p_2_2[[2]] + p_2_2[[3]] + plot_spacer() + p_2_1000[[1]] + p_2_1000[[2]] + p_2_1000[[3]] +
    p_2_2[[4]] + p_2_2[[5]] + p_2_2[[6]] + plot_spacer() + p_2_1000[[4]] + p_2_1000[[5]] + p_2_1000[[6]] +
    # j_phi_i_i_factor = 1000, m_i_factor = 2; j_phi_i_i_factor = 1000, m_i_factor = 1000
    group_label(2, 1000) + group_label(1000, 1000) +
    p_1000_2[[1]] + p_1000_2[[2]] + p_1000_2[[3]] + plot_spacer() + p_1000_1000[[1]] + p_1000_1000[[2]] + p_1000_1000[[3]] +
    p_1000_2[[4]] + p_1000_2[[5]] + p_1000_2[[6]] + plot_spacer() + p_1000_1000[[4]] + p_1000_1000[[5]] + p_1000_1000[[6]] +
    plot_layout(
      guides = "collect",
      axes = "collect_x",
      design = c(
        # text
        area(1, 1, 1, 7),
        # homeostasis
        area(2, 1), area(2, 2), area(2, 3),
        area(3, 1), area(3, 2), area(3, 3),
        # guide area
        area(3, 5, 3, 7),
        # spacer
        area(4, 1, 4, 7),
        # text: "End of inflammation"
        area(5, 1, 5, 7),
        # 2 texts: "Maturation ratio = 2..."
        area(6, 1, 6, 3), area(6, 4, 6, 7),
        # plots: j_phi_i_i_factor = 2, m_i_factor = 2; spacer; plots: j_phi_i_i_factor = 2, m_i_factor = 1000
        area(7, 1), area(7, 2), area(7, 3), area(7, 4), area(7, 5), area(7, 6), area(7, 7),
        area(8, 1), area(8, 2), area(8, 3), area(8, 4), area(8, 5), area(8, 6), area(8, 7),
        # 2 texts: "Maturation ratio = 2..."
        area(9, 1, 9, 3), area(9, 4, 9, 7),
        # j_phi_i_i_factor = 1000, m_i_factor = 2; spacer; j_phi_i_i_factor = 1000, m_i_factor = 1000
        area(10, 1), area(10, 2), area(10, 3), area(10, 4), area(10, 5), area(10, 6), area(10, 7),
        area(11, 1), area(11, 2), area(11, 3), area(11, 4), area(11, 5), area(11, 6), area(11, 7)
      ),
      heights = c(
        # text
        0.15,
        # homeostasis
        1, 1,
        0.0,
        0.15,
        0.21,
        1, 1,
        0.21,
        1, 1
      ),
      widths = c(1, 1, 1, 0.2, 1, 1, 1)
    ) &
    theme(
      #axis.title.y = element_text(size = 10, angle = 0, vjust = 1, hjust = 0, margin = margin(0, 0, 0, 0)),
      #axis.title.y = element_text(size = 10, vjust = 0, margin = margin(0, 2, 0, 0)),
      plot.margin = margin(5, 0, 0, 0),
      plot.title = element_text(size = 11, margin = margin(0, 0, 0, 0), hjust = -0.05)
    )

ggsave_with_defaults(
  plot = p_spatial_hom_vs_end_of_inf,
  paste(plot_dir_base, "spatial_hom_vs_end_of_inf.pdf", sep = "/"),
  #paste(plot_dir_base, "spatial_hom_vs_end_of_inf_fold_change.pdf", sep = "/"),
  #paste(plot_dir_base, "spatial_hom_vs_end_of_inf_difference.pdf", sep = "/"),
  width = 16,
  height = 10
)
options(scipen = 0)
