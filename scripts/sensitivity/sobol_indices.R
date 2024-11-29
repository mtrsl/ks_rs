library(sensobol)
library(ggplot2)
library(cowplot)
library(patchwork)
library(GGally)
library(ggnewscale)
library(ggh4x)
library(ggbreak)

source("scripts/sensitivity/functions.R")

res_dir_base <- "res_sensitivity_1000_pe_5_traces_only"

# Set and create the plot directory if it doesn't exist
plot_dir <- paste(res_dir_base, "plots", sep = "/")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Read the param_sample object from disk. This contains the initial (Sobol
# sequence or uniformly random) parameter samples and the resulting "design
# matrix" (`x$X`)
x <- readRDS(paste(res_dir_base, "param_sample.rds", sep = "/"))
# The size of the original parameter samples, not the total number of runs
n_param_sample <- x$n_param_sample
param_sample <- x$param_sample
param_sample_dt <- x$param_sample_dt

trace_data_full <- read_trace_data(param_sample_dt, res_dir_base)
trace_data_full[, cell_outflux := `-F_{phi_{C_u}}(x=0)` + `-F_{phi_{C_b}}(x=0)` + `-F_{phi_{C_s}}(x=0)`]

trace_data <- trace_data_full[`t_{inf}` >= 0]

# add dimensional time so we can plot against it
trace_data[, t_inf_h := 0.6944444444 * `t_{inf}`]

# vertical line at the time of end of inflammation (in hours)
end_of_inf_vline_h <- geom_vline(xintercept = 50 * 0.6944444444, alpha = 0.5)

# Trace plots
#############

rep_sample <- sample(0:nrow(param_sample_dt), min(nrow(param_sample_dt), 1000)) |> sort()

# Just plot the total amounts of the variables first
# TODO this is hacky and involves copying a lot of code - make it better?

trace_data_longer <- trace_data[
  rep %in% rep_sample,
  .(
    rep, t, t_inf_h,
    `C_u^{tot}`, `C_b^{tot}`, `C_s^{tot}`, `phi_{C_u}^{tot}`, `phi_{C_b}^{tot}`, `phi_{C_s}^{tot}`, `J^{tot}`,
    j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag
  )
  ] |>
  melt(
    measure.vars = c(
      "C_u^{tot}", "C_b^{tot}", "C_s^{tot}", "phi_{C_u}^{tot}", "phi_{C_b}^{tot}", "phi_{C_s}^{tot}", "J^{tot}"
    )
  ) |>
  melt(
    measure.vars = c(
      "j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag"
    ),
    variable.name = "param",
    value.name = "param_value"
  )

p_trace_grid <- ggplot(
  trace_data_longer,
  aes(
    x = t_inf_h,
    y = value,
  )
)
for (param_name in trace_data_longer$param |> levels()) {
  p_trace_grid <- p_trace_grid +
    geom_line(
      aes(
        group = rep,
        colour = param_value
      ),
      alpha = 0.2,
      linewidth = 0.2,
      data = trace_data_longer[param == param_name]
    ) +
    geom_vline(xintercept = 50 * 0.6944444444, alpha = 0.5) +
    colour_scales[param_name] +
    labs(colour = param_labels_words_no_breaks[param_name]) +
    new_scale_colour()
}
p_trace_grid <- p_trace_grid +
  labs(x = "Time since inflammation (h)", y = NULL) +
  facet_grid(
    rows = vars(variable),
    cols = vars(param),
    scales = "free",
    switch = "y",
    labeller = as_labeller(function(v) all_labels[v], label_parsed)
  ) +
  theme_cowplot() +
  theme(
    plot.background = element_rect("white"),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.box.just = "bottom",
    legend.position = "top",
    legend.justification = "centre",
    # no idea why 1/24 seems to be correct here given that "npc" units mean
    # that surely it should be 0.25, but who cares
    legend.key.width = unit(1.0 / 24.0, "npc"),
    legend.key.height = unit(0.3, "cm")
  )

ggsave_with_defaults(
  paste(plot_dir, "trace_grid_totals.png", sep = "/"),
  plot = p_trace_grid
)

# Now plot just the fluxes
trace_data_longer <- trace_data[
  rep %in% rep_sample,
  .(
    rep, t, t_inf_h,
    `-F_{phi_{C_u}}(x=0)`, `-F_{phi_{C_b}}(x=0)`, `-F_{phi_{C_s}}(x=0)`, cell_outflux,
    j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag
  )
  ] |>
  melt(
    measure.vars = c(
      "-F_{phi_{C_u}}(x=0)", "-F_{phi_{C_b}}(x=0)", "-F_{phi_{C_s}}(x=0)", "cell_outflux"
    )
  ) |>
  melt(
    measure.vars = c(
      "j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag"
    ),
    variable.name = "param",
    value.name = "param_value"
  )

p_trace_grid <- ggplot(
  trace_data_longer,
  aes(
    x = t_inf_h,
    y = value,
  )
)
for (param_name in trace_data_longer$param |> levels()) {
  p_trace_grid <- p_trace_grid +
    geom_line(
      aes(
        group = rep,
        colour = param_value
      ),
      alpha = 0.2,
      linewidth = 0.2,
      data = trace_data_longer[param == param_name]
    ) +
    geom_vline(xintercept = 50 * 0.6944444444, alpha = 0.5) +
    colour_scales[param_name] +
    labs(colour = param_labels_words_no_breaks[param_name]) +
    new_scale_colour()
}
p_trace_grid <- p_trace_grid +
  labs(x = "Time since inflammation (h)", y = NULL) +
  facet_grid(
    rows = vars(variable),
    cols = vars(param),
    scales = "free",
    switch = "y",
    labeller = as_labeller(function(v) all_labels[v], label_parsed)
  ) +
  theme_cowplot() +
  theme(
    plot.background = element_rect("white"),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.box.just = "bottom",
    legend.position = "top",
    legend.justification = "centre",
    # no idea why 1/24 seems to be correct here given that "npc" units mean
    # that surely it should be 0.25, but who cares
    legend.key.width = unit(1.0 / 24.0, "npc"),
    legend.key.height = unit(0.3, "cm")
  )

ggsave_with_defaults(
  paste(plot_dir, "trace_grid_fluxes.png", sep = "/"),
  plot = p_trace_grid
)

# ************************************************************
# Plots also including pe and gamma - just for "special" runs!
# ************************************************************

# Just plot the total amounts of the variables first
trace_data_longer <- trace_data[
  rep %in% rep_sample,
  .(
    rep, t, t_inf_h,
    `C_u^{tot}`, `C_b^{tot}`, `C_s^{tot}`, `phi_{C_u}^{tot}`, `phi_{C_b}^{tot}`, `phi_{C_s}^{tot}`, `J^{tot}`,
    j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, pe, gamma
  )
  ] |>
  melt(
    measure.vars = c(
      "C_u^{tot}", "C_b^{tot}", "C_s^{tot}", "phi_{C_u}^{tot}", "phi_{C_b}^{tot}", "phi_{C_s}^{tot}", "J^{tot}"
    )
  ) |>
  melt(
    measure.vars = c(
      "j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "pe", "gamma"
    ),
    variable.name = "param",
    value.name = "param_value"
  )

p_trace_grid <- ggplot(
  trace_data_longer,
  aes(
    x = t_inf_h,
    y = value,
  )
)
for (param_name in trace_data_longer$param |> levels()) {
  p_trace_grid <- p_trace_grid +
    geom_line(
      aes(
        group = rep,
        colour = param_value
      ),
      alpha = 0.2,
      linewidth = 0.5,
      data = trace_data_longer[param == param_name]
    ) +
    geom_vline(xintercept = 50 * 0.6944444444, alpha = 0.5) +
    colour_scales[param_name] +
    labs(colour = param_labels_words_no_breaks[param_name]) +
    new_scale_colour()
}
p_trace_grid <- p_trace_grid +
  labs(x = "Time since inflammation (h)", y = NULL) +
  facet_grid(
    rows = vars(variable),
    cols = vars(param),
    scales = "free",
    switch = "y",
    labeller = as_labeller(function(v) all_labels[v], label_parsed)
  ) +
  theme_cowplot() +
  theme(
    plot.background = element_rect("white"),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.box.just = "bottom",
    legend.position = "top",
    legend.justification = "centre",
    # no idea why 1/24 seems to be correct here given that "npc" units mean
    # that surely it should be 0.25, but who cares
    legend.key.width = unit(1.0 / 35.0, "npc"),
    legend.key.height = unit(0.3, "cm")
  )

ggsave_with_defaults(
  paste(plot_dir, "trace_grid_totals_pe_gamma.png", sep = "/"),
  plot = p_trace_grid
)

# Now plot just the fluxes
trace_data_longer <- trace_data[
  rep %in% rep_sample,
  .(
    rep, t, t_inf_h,
    `-F_{phi_{C_u}}(x=0)`, `-F_{phi_{C_b}}(x=0)`, `-F_{phi_{C_s}}(x=0)`, cell_outflux,
    j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, pe, gamma
  )
  ] |>
  melt(
    measure.vars = c(
      "-F_{phi_{C_u}}(x=0)", "-F_{phi_{C_b}}(x=0)", "-F_{phi_{C_s}}(x=0)", "cell_outflux"
    )
  ) |>
  melt(
    measure.vars = c(
      "j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "pe", "gamma"
    ),
    variable.name = "param",
    value.name = "param_value"
  )

p_trace_grid <- ggplot(
  trace_data_longer,
  aes(
    x = t_inf_h,
    y = value,
  )
)
for (param_name in trace_data_longer$param |> levels()) {
  p_trace_grid <- p_trace_grid +
    geom_line(
      aes(
        group = rep,
        colour = param_value
      ),
      alpha = 0.5,
      linewidth = 0.5,
      data = trace_data_longer[param == param_name]
    ) +
    geom_vline(xintercept = 50 * 0.6944444444, alpha = 0.5) +
    colour_scales[param_name] +
    labs(colour = param_labels_words_no_breaks[param_name]) +
    new_scale_colour()
}
p_trace_grid <- p_trace_grid +
  labs(
    x = "Time since inflammation (h)",
    y = NULL
  ) +
  facet_grid(
    rows = vars(variable),
    cols = vars(param),
    scales = "free",
    switch = "y",
    labeller = as_labeller(function(v) all_labels[v], label_parsed)
  ) +
  theme_cowplot() +
  theme(
    plot.background = element_rect("white"),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.box.just = "bottom",
    legend.position = "top",
    legend.justification = "centre",
    # no idea why 1/24 seems to be correct here given that "npc" units mean
    # that surely it should be 0.25, but who cares
    legend.key.width = unit(1.0 / 35.0, "npc"),
    legend.key.height = unit(0.3, "cm")
  )

ggsave_with_defaults(
  paste(plot_dir, "trace_grid_fluxes_pe_gamma.png", sep = "/"),
  plot = p_trace_grid
)

# END of hacky trace plot bit


# Sobol indices as a function of time - general methods

t_inf_max <- trace_data |>
  _[, .(t_stop = max(`t_{inf}`)), by = rep] |>
  _[, min(t_stop)]

output_inf_max <- trace_data |>
  _[, .(output_stop = max(output_inf)), by = rep] |>
  _[, min(output_stop)]

# checks that the t_inf values for different reps are essentially equal
check_t_infs <- function() {
  get_t_inf <- function(i) {
    trace_data[rep == i & `t_{inf}` <= t_inf_max, .(`t_{inf}`)]
  }

  compare_t_infs <- function(i, j) {
    (abs(get_t_inf(i) - get_t_inf(j)) <= 1e-4) |> all()
  }

  for (i in 2:(nrow(param_sample_dt) - 1)) {
    print(i)
    if (!compare_t_infs(1, i)) {
      print(paste0("1 and ", i, " t_infs differ more than 1e-4"))
    }
  }
}

# check_t_infs()

sobol_at_time <- function(data, output, variable) {
  y <- data[output_inf == output, eval(substitute(variable))]
  ind <- sobol_indices(
    Y = y,
    N = n_param_sample,
    params = param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  ind$results
}

p_sobol_vs_time <- function(sobol_data, title, ylim) {
  ggplot(
    sobol_data,
    # I forgot to add the dimensional times when calculating sobol indices, so manually adjust here
    aes(x = `t_{inf}` * 0.6944444444, y = original, group = parameters, colour = parameters)
  ) +
    geom_ribbon(
      aes(
        ymin = low.ci,
        ymax = high.ci,
        fill = parameters,
        colour = NULL
      ),
      alpha = 0.2
    ) +
    geom_line() +
    geom_vline(xintercept = 50 * 0.6944444444, alpha = 0.5) +
    labs(
      x = "Time since inflammation (h)",
      y = title,
      colour = "Parameter",
      fill = "Parameter"
    ) +
    scale_colour_discrete(labels = param_labels_words_no_breaks) +
    scale_fill_discrete(labels = param_labels_words_no_breaks) +
    coord_cartesian(ylim = ylim) +
    theme_cowplot() +
    theme(legend.text.align = 0)
}

# Sobol indices of cell outflux as a function of time

flux_sobol_indices <- list()

if (!file.exists(paste(res_dir_base, "flux_sobol_indices.rds", sep = "/"))) {
  for (i in 1:output_inf_max) {
    print(i)
    results <- sobol_at_time(trace_data, i, cell_outflux)

    t_inf <- trace_data[output_inf == i & rep == 1, `t_{inf}`]
    results[, `t_{inf}` := t_inf]

    flux_sobol_indices[[i]] <- results
  }

  flux_sobol_indices <- rbindlist(flux_sobol_indices)

  saveRDS(
    flux_sobol_indices,
    paste(res_dir_base, "flux_sobol_indices.rds", sep = "/")
  )
} else {
  flux_sobol_indices <- readRDS(
    paste(res_dir_base, "flux_sobol_indices.rds", sep = "/")
  )
}

# First order

p_flux_first_order <- p_sobol_vs_time(
  flux_sobol_indices[sensitivity == "Si"],
  "First-order Sobol index",
  c(-0.1, 1.2)
)

ggsave_with_defaults(
  plot = p_flux_first_order,
  paste(plot_dir, "flux_first_order_sobol_indices.pdf", sep = "/"),
)

# Total order

p_flux_total_order <- p_sobol_vs_time(
  flux_sobol_indices[sensitivity == "Ti"],
  "Total-order Sobol index",
  c(0, 1.05)
)

ggsave_with_defaults(
  plot = p_flux_total_order,
  paste(plot_dir, "flux_total_order_sobol_indices.pdf", sep = "/")
)

# Calculate and plot the sums of the indices for each time
flux_sobol_sums <- flux_sobol_indices[
  ,
  .(sum = sum(original)),
  by = c("sensitivity", "t_{inf}")
]

ggplot(flux_sobol_sums, aes(x = `t_{inf}`, y = sum, colour = sensitivity)) +
  geom_point() +
  #scale_y_continuous(limits = c(-1, 2)) +
  theme_cowplot() +
  background_grid()

################

# Sobol indices of maximum gradient of C_b as a function of time

source("scripts/extract_max_cell_density.R")

# Sometimes starting at output_inf == 1 here doesn't work, presumably because
# the gradient is constant everywhere at this time making the Sobol index
# calculations fail? Start at output_inf == 2 to avoid this problem

# C_u

c_u_gradient_sobol_indices <- list()

if (!file.exists(paste(res_dir_base, "c_u_gradient_sobol_indices.rds", sep = "/"))) {
  for (i in 2:output_inf_max) {
    print(i)
    results <- sobol_at_time(max_dc_u_dx_all_and_params, i, `dC_u_dx`)

    t_inf <- trace_data[output_inf == i & rep == 1, `t_{inf}`]
    results[, `t_{inf}` := t_inf]

    c_u_gradient_sobol_indices[[i]] <- results
  }

  c_u_gradient_sobol_indices <- rbindlist(c_u_gradient_sobol_indices)

  saveRDS(
    c_u_gradient_sobol_indices,
    paste(res_dir_base, "c_u_gradient_sobol_indices.rds", sep = "/")
  )
} else {
  c_u_gradient_sobol_indices <- readRDS(
    paste(res_dir_base, "c_u_gradient_sobol_indices.rds", sep = "/")
  )
}

# C_b

c_b_gradient_sobol_indices <- list()

if (!file.exists(paste(res_dir_base, "c_b_gradient_sobol_indices.rds", sep = "/"))) {
  for (i in 2:output_inf_max) {
    print(i)
    results <- sobol_at_time(max_dc_b_dx_all_and_params, i, `dC_b_dx`)

    t_inf <- trace_data[output_inf == i & rep == 1, `t_{inf}`]
    results[, `t_{inf}` := t_inf]

    c_b_gradient_sobol_indices[[i]] <- results
  }

  c_b_gradient_sobol_indices <- rbindlist(c_b_gradient_sobol_indices)

  saveRDS(
    c_b_gradient_sobol_indices,
    paste(res_dir_base, "c_b_gradient_sobol_indices.rds", sep = "/")
  )
} else {
  c_b_gradient_sobol_indices <- readRDS(
    paste(res_dir_base, "c_b_gradient_sobol_indices.rds", sep = "/")
  )
}

# C_s

c_s_gradient_sobol_indices <- list()

if (!file.exists(paste(res_dir_base, "c_s_gradient_sobol_indices.rds", sep = "/"))) {
  for (i in 2:output_inf_max) {
    print(i)
    results <- sobol_at_time(max_dc_s_dx_all_and_params, i, `dC_s_dx`)

    t_inf <- trace_data[output_inf == i & rep == 1, `t_{inf}`]
    results[, `t_{inf}` := t_inf]

    c_s_gradient_sobol_indices[[i]] <- results
  }

  c_s_gradient_sobol_indices <- rbindlist(c_s_gradient_sobol_indices)

  saveRDS(
    c_s_gradient_sobol_indices,
    paste(res_dir_base, "c_s_gradient_sobol_indices.rds", sep = "/")
  )
} else {
  c_s_gradient_sobol_indices <- readRDS(
    paste(res_dir_base, "c_s_gradient_sobol_indices.rds", sep = "/")
  )
}

# First order

p_c_u_gradient_first_order <- p_sobol_vs_time(
  c_u_gradient_sobol_indices[sensitivity == "Si"],
  "First-order Sobol index",
  c(-0.1, 1.0)
)

p_c_b_gradient_first_order <- p_sobol_vs_time(
  c_b_gradient_sobol_indices[sensitivity == "Si"],
  "First-order Sobol index",
  c(-0.1, 1.0)
)

p_c_s_gradient_first_order <- p_sobol_vs_time(
  c_s_gradient_sobol_indices[sensitivity == "Si"],
  "First-order Sobol index",
  c(-0.1, 1.0)
)

ggsave_with_defaults(
  plot = p_c_u_gradient_first_order,
  paste(plot_dir, "c_u_gradient_first_order_sobol_indices.pdf", sep = "/")
)

ggsave_with_defaults(
  plot = p_c_b_gradient_first_order,
  paste(plot_dir, "c_b_gradient_first_order_sobol_indices.pdf", sep = "/")
)

ggsave_with_defaults(
  plot = p_c_s_gradient_first_order,
  paste(plot_dir, "c_s_gradient_first_order_sobol_indices.pdf", sep = "/")
)

# Total order

p_c_u_gradient_total_order <- p_sobol_vs_time(
  c_u_gradient_sobol_indices[sensitivity == "Ti"],
  "Total-order Sobol index",
  c(0.0, 1.05)
)

p_c_b_gradient_total_order <- p_sobol_vs_time(
  c_b_gradient_sobol_indices[sensitivity == "Ti"],
  "Total-order Sobol index",
  c(0.0, 1.05)
)

p_c_s_gradient_total_order <- p_sobol_vs_time(
  c_s_gradient_sobol_indices[sensitivity == "Ti"],
  "Total-order Sobol index",
  c(0.0, 1.05)
)

ggsave_with_defaults(
  plot = p_c_u_gradient_total_order,
  paste(plot_dir, "c_u_gradient_total_order_sobol_indices.pdf", sep = "/")
)

ggsave_with_defaults(
  plot = p_c_b_gradient_total_order,
  paste(plot_dir, "c_b_gradient_total_order_sobol_indices.pdf", sep = "/")
)

ggsave_with_defaults(
  plot = p_c_s_gradient_total_order,
  paste(plot_dir, "c_s_gradient_total_order_sobol_indices.pdf", sep = "/")
)

# Max value and location of max value plots for phi_c_b and dc_b_dx
# -----------------------------------------------------------------

# Some plotting styles and helper functions

p_location_panel <- function(data, colour_by, colour_style) {
  print(deparse(substitute(colour_by)))
  print(param_labels_words_no_breaks[deparse(substitute(colour_by))])

  ggplot(
    data,
    aes(time_inf, x, group = rep, colour = {{ colour_by }})
    ) +
  geom_path(alpha = 0.1) +
  coord_cartesian(ylim = c(0, 1)) +
  colour_style +
  labs(
    x = "Time since inflammation",
    y = expression(x),
    colour = param_labels_words_no_breaks[deparse(substitute(colour_by))]
  )
}

p_value_panel <- function(data, colour_by, colour_style, ylabel) {
  print(deparse(substitute(colour_by)))
  print(param_labels_words_no_breaks[deparse(substitute(colour_by))])

  ggplot(
    data,
    aes(time_inf, phi_C_b, group = rep, colour = {{ colour_by }})
    ) +
  geom_path(alpha = 0.1) +
  colour_style +
  labs(
    x = "Time since inflammation",
    y = ylabel,
    colour = param_labels_words_no_breaks[deparse(substitute(colour_by))]
  )
}

# Location of max phi_c_b

p_j_phi_i_i_factor <- p_location_panel(max_phi_c_b_all_and_params, j_phi_i_i_factor, blue)
p_m_i_factor <- p_location_panel(max_phi_c_b_all_and_params, m_i_factor, red)
p_t_j_phi_i_lag <- p_location_panel(max_phi_c_b_all_and_params, t_j_phi_i_lag, green)

p_phi_c_b_location <- p_j_phi_i_i_factor + p_m_i_factor + p_t_j_phi_i_lag +
  plot_annotation(
    title = expression(
      paste("Location of maxmimum ", phi[C[b]], " concentration")
    )
  ) &
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect("white")
  )

ggsave_with_defaults(
  plot = p_phi_c_b_location,
  paste(plot_dir, "max_phi_c_b_location.png", sep = "/")
)

# Value of max phi_c_b

p_phi_c_b_value_panel <- function(...) {
  p_value_panel(..., expression(phi[C[b]]))
}

p_j_phi_i_i_factor <- p_phi_c_b_value_panel(max_phi_c_b_all_and_params, j_phi_i_i_factor, blue)
p_m_i_factor <- p_phi_c_b_value_panel(max_phi_c_b_all_and_params, m_i_factor, red)
p_t_j_phi_i_lag <- p_phi_c_b_value_panel(max_phi_c_b_all_and_params, t_j_phi_i_lag, green)

p_phi_c_b_value <- p_j_phi_i_i_factor + p_m_i_factor + p_t_j_phi_i_lag +
  plot_annotation(
    title = expression(
      paste("Value of maxmimum ", phi[C[b]], " concentration")
    )
  ) &
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect("white")
  )

ggsave_with_defaults(
  plot = p_phi_c_b_value,
  paste(plot_dir, "max_phi_c_b_value.png", sep = "/")
)

# Location of max gradient of c_b

p_j_phi_i_i_factor <- p_location_panel(max_dc_b_dx_all_and_params, j_phi_i_i_factor, blue)
p_m_i_factor <- p_location_panel(max_dc_b_dx_all_and_params, m_i_factor, red)
p_t_j_phi_i_lag <- p_location_panel(max_dc_b_dx_all_and_params, t_j_phi_i_lag, green)

p_dc_b_dx_location <- p_j_phi_i_i_factor + p_m_i_factor + p_t_j_phi_i_lag +
  plot_annotation(
    title = expression(
      paste("Location of maxmimum ", C[b], " gradient")
    )
  ) &
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect("white")
  )

ggsave_with_defaults(
  plot = p_dc_b_dx_location,
  paste(plot_dir, "max_dc_b_dx_location.png", sep = "/")
)

# Value of max gradient of c_b

p_dc_b_dx_value_panel <- function(...) {
  p_value_panel(..., expression(dC[b]/dx))
}

p_j_phi_i_i_factor <- p_dc_b_dx_value_panel(max_dc_b_dx_all_and_params, j_phi_i_i_factor, blue)
p_m_i_factor <- p_dc_b_dx_value_panel(max_dc_b_dx_all_and_params, m_i_factor, red)
p_t_j_phi_i_lag <- p_dc_b_dx_value_panel(max_dc_b_dx_all_and_params, t_j_phi_i_lag, green)

p_dc_b_dx_value <- p_j_phi_i_i_factor + p_m_i_factor + p_t_j_phi_i_lag +
  plot_annotation(
    title = expression(
      paste("Value of maxmimum ", C[b], " gradient")
    )
  ) &
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect("white")
  )

ggsave_with_defaults(
  plot = p_dc_b_dx_value,
  paste(plot_dir, "max_dc_b_dx_value.png", sep = "/")
)

################

# Get the rows corresponding to the (global) maximum cell flux across lymphatic
# vessel
trace_max_flux <- trace_data[, .SD[which.max(cell_outflux)], by = rep]

trace_max_flux_long <- melt(
  trace_max_flux,
  measure.vars = c(
    "j_phi_i_i_factor",
    "m_i_factor",
    "t_j_phi_i_lag"),
  variable.name = "param",
  value.name = "param_value"
)

# Get the first local maximum assuming it's before t_inf = 50
trace_max_flux_first <-
  trace_data[`t_{inf}` < 49, .SD[which.max(cell_outflux)], by = rep]

trace_max_flux_first_long <- melt(
  trace_max_flux_first,
  measure.vars = c(
    "j_phi_i_i_factor",
    "m_i_factor",
    "t_j_phi_i_lag"),
  variable.name = "param",
  value.name = "param_value"
)

# Get the second local maximum assuming it's after t_inf = 50
trace_max_flux_second <-
  trace_data[`t_{inf}` >= 49, .SD[which.max(cell_outflux)], by = rep]

trace_max_flux_second_long <- melt(
  trace_max_flux_second,
  measure.vars = c(
    "j_phi_i_i_factor",
    "m_i_factor",
    "t_j_phi_i_lag"),
  variable.name = "param",
  value.name = "param_value"
)

plot_max_flux <- function(max_data) {
  max_flux_panel <- function(param) {
    ggplot(
      trace_max_flux,
      aes(x = `t_inf_h`, y = cell_outflux, colour = {{ param }})
    ) +
      geom_point(size = 2) +
      geom_vline(xintercept = 50 * 0.6944444444, alpha = 0.5) +
      #scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 5)) +
      #scale_x_break(c(10, 48)) +
      labs(
        x = "Time since inflammation (h)",
        y = variable_labels["cell_outflux"],
        colour = param_labels_words[deparse(substitute(param))]
      ) +
      theme_cowplot() +
      background_grid() +
      theme(legend.position = "top")
  }

  p_1 <- max_flux_panel(j_phi_i_i_factor) + blue
  p_2 <- max_flux_panel(m_i_factor) + red
  p_3 <- max_flux_panel(t_j_phi_i_lag) + green

  p_1 + p_2 + p_3 +
    plot_layout(nrow = 3, ncol = 1, guides = "collect", axes = "collect_x") &
    theme(
      legend.position = "top",
      legend.box.just = "bottom",
      legend.justification = "centre",
      legend.key.width = unit(1.0 / 24.0, "npc"),
      legend.key.height = unit(0.3, "cm")
    )
}

p_max_flux <- plot_max_flux(trace_max_flux)
#p_max_flux <- plot_max_flux(trace_max_flux_second)

ggsave_with_defaults(
  plot = p_max_flux,
  paste(plot_dir, "max_flux.pdf", sep = "/")
)

trace_max_flux_full_traj <-
  trace_data[rep %in% trace_max_flux[`t_{inf}` < 20, rep]]

ggplot(
  trace_max_flux_full_traj,
  aes(x = `t_{inf}`, y = cell_outflux, group = rep)
) +
  geom_line()

#################################

# Find all local maxima and minima of fluxes

# Helper functions
find_local_maxima_indices <- function(values) {
  max_ind <- NULL

  for (i in 2:(length(values) - 1)) {
    if (values[i] > values[i - 1] && values[i] > values[i + 1]) {
      max_ind <- c(max_ind, i)
    }
  }
  max_ind
}

find_local_minima_indices <- function(values) {
  min_ind <- NULL

  for (i in 2:(length(values) - 1)) {
    if (values[i] < values[i - 1] && values[i] < values[i + 1]) {
      min_ind <- c(min_ind, i)
    }
  }
  min_ind
}

flux_trace_plot_subset <- function(subset_rep_ids) {
  ggplot(
    trace_data[rep %in% subset_rep_ids],
    aes(
      x = t_inf_h,
      y = cell_outflux,
      group = rep
    )
  ) +
  geom_line() +
  theme_cowplot() +
  background_grid() +
  labs(
    x = "Time since inflammation (h)",
    y = "Cell flux at l.v."
  )
}

parameter_pairs_plot <- function(data) {
  ggpairs(
    data,
    columns = c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag"),
    upper = "blank",
    diag = list(continuous = "barDiag"),
    progress = FALSE,
    labeller = as_labeller(function(v) param_labels[v], label_parsed)
  ) +
  theme_cowplot() +
  background_grid()
}

# Find the local maxima of the flux for each rep
flux_local_maxima <- trace_data[,
  .SD[find_local_maxima_indices(cell_outflux)],
  by = rep
] |>
  _[, extrema_id := 1:.N, by = rep] |>
  _[, extrema_type := factor("maximum")]

# Find the local minima of the flux for each rep
flux_local_minima <- trace_data[,
  .SD[find_local_minima_indices(cell_outflux)],
  by = rep
] |>
  _[, extrema_id := 1:.N, by = rep] |>
  _[, extrema_type := factor("minimum")]

flux_local_extrema <- rbind(flux_local_maxima, flux_local_minima)

# Plot all extrema at once
ggplot(
  flux_local_extrema[extrema_id == 1],
  aes(
    x = t_inf_h,
    y = cell_outflux,
    colour = extrema_type,
    #shape = extrema_type,
    group = rep
  )
) +
  geom_point(size = 1) +
  geom_line(alpha = 0.2) +
  scale_x_continuous(limits = c(0, 39)) +
  scale_shape_manual(values = c(19, 6)) +
  theme_cowplot() +
  background_grid() +
  labs(
    x = "Time since inflammation (h)",
    y = "Cell flux at l.v.",
    colour = "Count",
    shape = "Type"
  ) #+
  #geom_point(data = flux_local_extrema[extrema_type == "minimum"])

# Plot all minima
ggplot(
  flux_local_minima,
  aes(x = `t_{inf}`, y = cell_outflux, colour = factor(extrema_id))
) +
  geom_point() +
  scale_x_continuous(limits = c(0, 55)) +
  theme_cowplot() +
  background_grid()

# Plot all maxima
ggplot(
  flux_local_maxima,
  aes(x = `t_{inf}`, y = cell_outflux, colour = factor(extrema_id))
) +
  geom_point() +
  geom_vline(xintercept = 8 - 2) +
  geom_vline(xintercept = 8 + 2) +
  scale_x_continuous(limits = c(0, 55)) +
  theme_cowplot() +
  background_grid()

# Find the reps that have only one local maximum
flux_single_max <- flux_local_maxima |>
  _[, .SD[which.max(extrema_id)], by = rep] |>
  _[extrema_id == 1]

# Find the reps that have only one local minimum
flux_single_min <- flux_local_minima |>
  _[, .SD[which.max(extrema_id)], by = rep] |>
  _[extrema_id == 1]

# Bindi's suggested max/min thresholds

# (i) Solutions that have only 1 max at early time, t1 (t1 = 20?)
flux_early_max_ids <- flux_single_max[`t_{inf}` < 20, rep]

# (ii) Solutions that have only 1 min at t2, where t1 < t2 < t3
# TODO what are sensible t1 and t3 here? in the previous case, the choice of t1
# makes no difference, as long as it's > 5...
flux_middle_min_ids <- flux_single_min[`t_{inf}` > 0, rep]

# (iii) Solutions that have only 1 max at the end of inflammation, t3 (t3 ~ 50)
flux_end_inf_max_ids <- flux_single_max[abs(`t_{inf}` - 50) < 1e-4, rep]

# (iv) Solutions that have only 1 max after the end of inflammation, "t4" (t_inf > 50 (+ a
# bit due to output times not strictly being evenly spaced)
flux_after_inf_max_ids <- flux_single_max[`t_{inf}` > 50.00002, rep]

flux_trace_plot_subset(flux_early_max_ids) +
  geom_point(
    data = flux_single_max[`t_{inf}` < 20],
    aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`),
    colour = "red",
    shape = "cross"
  )

flux_trace_plot_subset(flux_middle_min_ids) +
  geom_point(
    data = flux_single_min,
    aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`),
    colour = "red",
    shape = "cross"
  )

flux_trace_plot_subset(flux_end_inf_max_ids) +
  geom_point(
    data = flux_single_max[abs(`t_{inf}` - 50) < 1e-4],
    aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`),
    colour = "red",
    shape = "cross"
  )

flux_trace_plot_subset(flux_after_inf_max_ids) +
  geom_point(
    data = flux_single_max[`t_{inf}` > 50.00002],
    aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`),
    colour = "red",
    shape = "cross"
  )

parameter_pairs_plot(flux_single_min)
parameter_pairs_plot(trace_data[output_inf == 1])

#####

# Get the ids for the reps that have at least 2 maxima
reps_w_al_2_maxima <- flux_local_maxima[extrema_id == 2, rep]

# Get the actual rows from the flux table up to the 2nd maximum
flux_local_al_2_maxima <- flux_local_maxima[
  rep %in% reps_w_al_2_maxima & extrema_id <= 2
]

# Only keep the columns we need here. This makes transforming into wide form
# below easier because we don't need to account for columns that vary between
# rows for one rep...
flux_local_al_2_maxima <- flux_local_al_2_maxima[
  ,
  .(
    rep,
    `t_{inf}`,
    t_inf_h,
    j_phi_i_i_factor,
    m_i_factor,
    t_j_phi_i_lag,
    extrema_id,
    cell_outflux
  )
]

# dcast the two maxima into their own columns
flux_local_al_2_maxima_wide <- dcast(
  flux_local_al_2_maxima,
  ... ~ extrema_id,
  value.var = c("t_{inf}", "t_inf_h", "cell_outflux")
)

# calculate the ratio max_2/max_1
flux_local_al_2_maxima_wide[
  ,
  ratio_21 := cell_outflux_2 / cell_outflux_1,
  by = rep
]

# "join" the ratio_21 column with the trace_data
trace_data_with_max_21_ratio <-
  trace_data[rep %in% flux_local_al_2_maxima_wide[, rep]][
  flux_local_al_2_maxima_wide[, .(rep, ratio_21)],
  on = "rep"
]

## convert to long form
#trace_data_with_ratios_long <- melt(
  #trace_data_with_max_21_ratio,
  #measure.vars = c(
    #"C_u^{tot}",
    #"C_b^{tot}",
    #"C_s^{tot}",
    #"phi_i^{tot}",
    #"phi_m^{tot}",
    #"phi_{C_u}^{tot}",
    #"phi_{C_b}^{tot}",
    #"-F_{phi_i}(x=1)",
    #"-F_{phi_{C_b}}(x=0)"
  #)
#)

ggplot(
  flux_local_al_2_maxima_wide,
  aes(x = ratio_21, colour = j_phi_i_i_factor)
) +
  geom_histogram(bins = 50)
  #geom_density()
  #geom_point(position = position_jitter(seed = 1))

ggplot(
  trace_data_with_max_21_ratio,
  aes(x = t_inf_h, y = cell_outflux, colour = ratio_21, group = rep)
) +
  geom_line(alpha = 0.5, size = 0.1) +
  scale_colour_distiller(palette = "Spectral") +
  labs(
    x = "Time since inflammation (h)",
    y = "Cell flux at l.v.",
    colour = "Max 2 / Max 1"
  ) +
  theme_cowplot() +
  background_grid()

###########

# Find the ids of all reps that have no local minima

# First get the reps that do have minima
reps_w_some_minima <- flux_local_minima[, .(rep)] |> unique()

# Then take the complement. There's definitely a better way to do this but this
# works...
reps_w_no_minima <-
  data.table(rep = 0:trace_data[, max(rep)])[!reps_w_some_minima, rep, on = "rep"]

flux_trace_plot_subset(reps_w_no_minima)

# We want the subset of the above reps that have a maximum at (or maybe soon
# after?) t_inf = 50, regardless of the total number of maxima

reps_w_no_min_but_max_50 <-
  flux_local_maxima[rep %in% reps_w_no_minima & `t_{inf}` >= 50, rep]

flux_trace_plot_subset(reps_w_no_min_but_max_50)

parameter_pairs_plot(trace_data[rep %in% reps_w_no_min_but_max_50])

###########

# Work out, for reps that have a minimum before t_inf == 50 (i.e. somewhere
# during inflammation), the ratio of the minimum to the first maximum. The
# first condition finds all reps that have any kind of extrema before t = 50,
# including some reps that have say a maximum but no minimum. We filter these
# out using na.omit to leave just the reps we actually want.
flux_w_max_min_ratio <- flux_local_extrema[
  `t_{inf}` < 50 & extrema_id == 1,
  .(rep, cell_outflux, extrema_type)
] |>
  dcast(... ~ extrema_type, value.var = c("cell_outflux")) |>
  _[, ratio_max_min := maximum / minimum] |>
  na.omit()

ggplot(flux_w_max_min_ratio, aes(x = ratio_max_min)) +
  #geom_density() +
  geom_histogram(bins = 50) +
  coord_cartesian(xlim = c(0, NA))

reps_w_min_and_max_50 <- flux_w_max_min_ratio[, rep]

# "join" the ratio_max_min column with the trace_data
trace_data_with_max_min_ratio <-
  trace_data[rep %in% reps_w_min_and_max_50][
  flux_w_max_min_ratio[, .(rep, ratio_max_min)],
  on = "rep"
]

# Plot the flux vs time, coloured by the ratio of max to min
ggplot(
  trace_data_with_max_min_ratio,
  aes(
    x = `t_{inf}`,
    y = cell_outflux,
    group = rep,
    colour = ratio_max_min
  )
) +
  geom_line(alpha = 0.5, size = 1) +
  scale_colour_distiller(palette = "Spectral") + theme_cowplot() +
  background_grid()

# Plot the flux vs time for reps that have ratio below a threshold
ggplot(
  trace_data_with_max_min_ratio[ratio_max_min <= 1.1],
  aes(
    x = `t_{inf}`,
    y = cell_outflux,
    group = rep,
    colour = ratio_max_min
  )
) +
  geom_line(alpha = 0.5, size = 1) +
  scale_colour_distiller(palette = "Spectral") + theme_cowplot() +
  background_grid()

# UPDATE: these aren't actually the relevant reps at all (I misunderstood...)
# Combine the above subset of reps with those that have no minimum but still
# have a max at 50. These should be the relevant ones???
reps_relevant <- c(
  flux_w_max_min_ratio[ratio_max_min <= 1.1, rep],
  reps_w_no_min_but_max_50
)

ggplot(
  trace_data[rep %in% reps_relevant],
  aes(
    x = `t_{inf}`,
    y = cell_outflux,
    group = rep,
    #colour = ratio_max_min
  )
) +
  geom_line(alpha = 0.5, size = 1) +
  scale_colour_distiller(palette = "Spectral") + theme_cowplot() +
  background_grid()

parameter_pairs_plot(trace_data[rep %in% reps_relevant])

###########

# Find the actually relevant reps - those that have a first peak at around 6
# hours (need to decide size of range around 6h to include), then have a second
# peak with height not too much lower than the first.

first_max_range_mean <- 6

# The time either side of 6h we allow
first_max_range_width <- 2.78
ratio_21_threshold <- 0.75
relevant_reps <- flux_local_al_2_maxima_wide[
  first_max_range_mean - first_max_range_width <= `t_inf_h_1` &
    `t_inf_h_1` <= first_max_range_mean + first_max_range_width &
    ratio_21 >= ratio_21_threshold,
  rep
]
flux_trace_plot_subset(relevant_reps) +
  geom_rect(
    data = ~ head(.x, 1),
    xmin = first_max_range_mean - first_max_range_width,
    xmax = first_max_range_mean + first_max_range_width,
    ymin = -1,
    ymax = 1,
    colour = NA,
    fill = "black",
    alpha = 0.1
  ) +
  geom_point(
    data = flux_local_al_2_maxima_wide[rep %in% relevant_reps],
    aes(x = `t_inf_h_1`, y = cell_outflux_1),
    colour = "red"
  ) +
  geom_point(
    data = flux_local_al_2_maxima_wide[rep %in% relevant_reps],
    aes(x = `t_inf_h_2`, y = cell_outflux_2),
    colour = "blue"
  ) +
  geom_vline(xintercept = 50 * 0.6944444444, alpha = 0.5)

#ggsave(
  #plot = last_plot(),
  #"simulation_subset.png",
  #width = 9,
  #height = 5
#)

parameter_pairs_plot(
  flux_local_al_2_maxima_wide[rep %in% relevant_reps]
)

# TODO the next couple of bits were copied from earlier - refactor into
# functions if we keep it
trace_data_relevant_longer <- trace_data[
  rep %in% relevant_reps,
  .(
    rep, `t_{inf}`, t_inf_h,
    `C_u^{tot}`, `C_b^{tot}`, `C_s^{tot}`,
    `phi_{C_u}^{tot}`, `phi_{C_b}^{tot}`, `phi_{C_s}^{tot}`, 
    `cell_outflux`,
    j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag
  )
  ] |>
  melt(
    measure.vars = c(
      "C_u^{tot}", "C_b^{tot}", "C_s^{tot}",
      "phi_{C_u}^{tot}", "phi_{C_b}^{tot}", "phi_{C_s}^{tot}",
      "cell_outflux"
    )
  ) |>
  melt(
    measure.vars = c(
      "j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag"
    ),
    variable.name = "param",
    value.name = "param_value"
  )

p_relevant_trace_grid <- ggplot(
  trace_data_relevant_longer,
  aes(
    x = t_inf_h,
    y = value,
  )
)
for (param_name in trace_data_relevant_longer$param |> levels()) {
  p_relevant_trace_grid <- p_relevant_trace_grid +
    geom_line(
      aes(
        group = rep,
        colour = param_value
      ),
      alpha = 0.5,
      data = trace_data_relevant_longer[param == param_name]
    ) +
    colour_scales[param_name] +
    labs(colour = param_labels[param_name]) +
    new_scale_colour()
}
p_relevant_trace_grid <- p_relevant_trace_grid +
  labs(x = "Time since inflammation (h)", y = NULL) +
  facet_grid(
    rows = vars(variable),
    cols = vars(param),
    scales = "free",
    switch = "y",
    labeller = as_labeller(function(v) all_labels[v], label_parsed)
  ) +
  theme_cowplot() +
  theme(
    plot.background = element_rect("white"),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.box.just = "bottom",
    legend.position = "top",
    legend.justification = "centre",
    # no idea why 1/24 seems to be correct here given that "npc" units mean
    # that surely it should be 0.25, but who cares
    legend.key.width = unit(1.0 / 24.0, "npc"),
    legend.key.height = unit(0.3, "cm")
  )
p_relevant_trace_grid

# Comparing with all reps (with >= 2 max) with m_i_factor approx in the range
# indicated by the above subset

flux_trace_plot_subset(
  flux_local_al_2_maxima_wide[m_i_factor %between% c(50, 200), rep]
)

setdiff(
  relevant_reps,
  flux_local_al_2_maxima_wide[m_i_factor %between% c(50, 200), rep]
)

# Plot ratio_21 against parameters

flux_local_al_2_maxima_wide_long_params <- melt(
  flux_local_al_2_maxima_wide,
  measure.vars = c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag"),
  variable.name = "param",
  value.name = "param_value"
)

quantity_vs_param_plot <- function(quantity, y_label, colour_by) {
  p <- ggplot(
    # this join copies the un-melted parameter columns into the melted table so
    # we can use their values for e.g. colouring
    flux_local_al_2_maxima_wide_long_params[
      flux_local_al_2_maxima_wide[,
        .(rep, j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag)
      ],
      on = "rep"
    ],
    #][`t_{inf}_2` > 45], # & m_i_factor <= 200],
    aes(x = param_value, y = {{ quantity }}, colour = {{ colour_by }})
  ) +
    geom_point(size = 4, alpha = 0.6) +
    #geom_hline(yintercept = 1) +
    facet_wrap(
      vars(param),
      scales = "free",
      labeller = as_labeller(function(v) param_labels_words_no_breaks[v]),
      strip.position = "bottom"
    ) +
    scale_colour_distiller(palette = "Spectral") +
    scale_x_continuous(trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_continuous(trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    labs(
      x = NULL,
      y = y_label,
      colour = param_labels_words_no_breaks[deparse(substitute(colour_by))]
    ) +
    theme_cowplot() +
    background_grid() +
    theme(
      strip.placement = "outside",
      strip.background = element_blank()
    )

  q_str <- gsub("[{}]", "", deparse(substitute(quantity)))
  c_str <- gsub("[{}]", "", deparse(substitute(colour_by)))
  plot_name <- paste0(q_str, "_by_", c_str, ".pdf")

  q_vs_p_plot_dir <- paste(plot_dir, "quantity_vs_param", sep = "/")

  if (!dir.exists(q_vs_p_plot_dir)) {
    dir.create(q_vs_p_plot_dir, recursive = TRUE)
  }

  ggsave_with_defaults(
    plot = p,
    paste(plot_dir, "quantity_vs_param", plot_name, sep = "/")
  )

  p
}

quantity_vs_param_plot(t_inf_h_1, "Time of first peak (h)", j_phi_i_i_factor)
quantity_vs_param_plot(t_inf_h_1, "Time of first peak (h)", m_i_factor)
quantity_vs_param_plot(t_inf_h_1, "Time of first peak (h)", t_j_phi_i_lag)

quantity_vs_param_plot(t_inf_h_2, "Time of second peak (h)", j_phi_i_i_factor)
quantity_vs_param_plot(t_inf_h_2, "Time of second peak (h)", m_i_factor)
quantity_vs_param_plot(t_inf_h_2, "Time of second peak (h)", t_j_phi_i_lag)

quantity_vs_param_plot(cell_outflux_1, "Cell LV flux at second peak", j_phi_i_i_factor)
quantity_vs_param_plot(cell_outflux_1, "Cell LV flux at second peak", m_i_factor)
quantity_vs_param_plot(cell_outflux_1, "Cell LV flux at second peak", t_j_phi_i_lag)

quantity_vs_param_plot(cell_outflux_2, "Cell LV flux at second peak", j_phi_i_i_factor)
quantity_vs_param_plot(cell_outflux_2, "Cell LV flux at second peak", m_i_factor)
quantity_vs_param_plot(cell_outflux_2, "Cell LV flux at second peak", t_j_phi_i_lag)

#####################

# Basic flux plot
p_basic_flux <- ggplot(
  trace_data_with_max_21_ratio,
  aes(
    x = t_inf_h,
    y = cell_outflux,
    group = rep,
    colour = j_phi_i_i_factor
  )
) +
  geom_line(alpha = 0.4, linewidth = 0.6) +
  geom_vline(xintercept = 50 * 0.6944444444, alpha = 0.5, linewidth = 1) +
  scale_colour_distiller(palette = "Spectral") +
  theme_cowplot(font_size = 26) +
  theme(
    legend.position = c(0.80, 0.6),
    legend.title = element_text(size = 22)
  ) +
  labs(
    x = "Time since inflammation",
    y = "DC flux into l.v.",
    colour = "DC ingress\nratio",
    #tag = "B"
  )
p_basic_flux

ggsave(
  paste(plot_dir, "basic_flux.pdf", sep = "/"),
  plot = p_basic_flux,
  width = 9,
  height = 5,
  device = cairo_pdf
)

ggsave(
  paste(plot_dir, "basic_flux.png", sep = "/"),
  plot = p_basic_flux,
  width = 9,
  height = 5,
  dpi = 150,
  bg = "white"
)

################

# Basic flux at 2nd peak vs j_phi_i_i_factor coloured by gamma plot
p_basic_peak_vs_ingress_vs_gamma <- ggplot(
  flux_local_al_2_maxima[extrema_id == 2],
  aes(x = j_phi_i_i_factor, y = cell_outflux, colour = gamma)
) +
  geom_point(size = 4, alpha = 0.6) +
  scale_colour_distiller(palette = "Spectral", guide = guide_colorbar(title.position = "left", title.hjust = 1.0, barheight = 7)) +
  labs(
    x = "Ingress ratio",
    y = "DC flux at 2nd peak",
    colour = "Cleavage\nrate",
    #tag = "C"
  ) +
  theme_cowplot(font_size = 26) +
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22), legend.position = c(0.65, 0.23))
p_basic_peak_vs_ingress_vs_gamma
ggsave(
  paste(plot_dir, "basic_peak_vs_ingress_vs_gamma.pdf", sep = "/"),
  plot = p_basic_peak_vs_ingress_vs_gamma,
  width = 9,
  height = 5,
  device = cairo_pdf
)
