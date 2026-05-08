library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

## Import data -----------------------------------------------------------------
R <- readRDS("data/expected_loss/tabearly/post_mean/baytpgrpesreject1e+06.rds")

## Helpers ---------------------------------------------------------------------
stop_index_first <- function(x) {
  w <- which(x != 0L)
  if (length(w)) w[1] else NA_integer_
}

stop_sign_first <- function(x) {
  w <- which(x != 0L)
  if (length(w)) sign(x[w[1]]) else 0L
}

## Dimensions ------------------------------------------------------------------
n_sim   <- dim(R)[1]
K       <- dim(R)[2]
n_theta <- dim(R)[3]
n_prior <- dim(R)[4]

## Use first prior only --------------------------------------------------------
prior_idx <- 1L

## Apply helpers ---------------------------------------------------------------
first_looks_all <- apply(
  R,
  MARGIN = c(1, 3, 4),
  FUN = function(x) stop_index_first(x)
)

first_signs_all <- apply(
  R,
  MARGIN = c(1, 3, 4),
  FUN = function(x) stop_sign_first(x)
)

## Extract first prior only: matrix [sim x theta] ------------------------------
first_looks_p1 <- first_looks_all[, , prior_idx, drop = FALSE]
first_signs_p1 <- first_signs_all[, , prior_idx, drop = FALSE]

first_looks_p1 <- matrix(first_looks_p1, nrow = n_sim, ncol = n_theta)
first_signs_p1 <- matrix(first_signs_p1, nrow = n_sim, ncol = n_theta)

## Theta labels ----------------------------------------------------------------
## Edit this if your theta values differ
theta_labels <- c("-δ/2", "0", "δ/2", "δ", "3δ/2", "2δ")

if (length(theta_labels) != n_theta) {
  stop("theta_labels must have length n_theta")
}

## Build tidy simulation grid --------------------------------------------------
grid <- expand.grid(
  sim   = seq_len(n_sim),
  theta = seq_len(n_theta)
)

grid$stop_look <- as.integer(c(first_looks_p1))
grid$stop_type <- as.integer(c(first_signs_p1))

grid <- grid %>%
  mutate(
    Theta = factor(theta_labels[theta], levels = theta_labels),
    stop_look = as.integer(stop_look)
  )

## Generic summariser ----------------------------------------------------------
summarise_stops <- function(df, stop_code, label) {
  df %>%
    filter(stop_type == stop_code) %>%
    count(Theta, stop_look) %>%
    complete(
      Theta,
      stop_look = seq_len(K),
      fill = list(n = 0)
    ) %>%
    mutate(
      prop = n / n_sim,
      type = label
    )
}

## Combine efficacy and futility into one row per tile -------------------------
df_plot <- bind_rows(
  summarise_stops(grid,  1, "efficacy"),
  summarise_stops(grid, -1, "futility")
) %>%
  select(Theta, stop_look, type, prop) %>%
  pivot_wider(
    names_from  = type,
    values_from = prop,
    values_fill = 0
  ) %>%
  rename(
    eff_prop = efficacy,
    fut_prop = futility
  )

## Y positions -----------------------------------------------------------------
theta_levels <- levels(grid$Theta)

df_plot <- df_plot %>%
  mutate(
    y = match(as.character(Theta), theta_levels),
    eff_fill = alpha("green4", eff_prop),
    fut_fill = alpha("red3",   fut_prop)
  )

## Pretty y-axis labels --------------------------------------------------------
y_labs <- as.expression(c(
  bquote(-delta/2),
  bquote(0),
  bquote(delta/2),
  bquote(delta),
  bquote(3 * delta / 2),
  bquote(2 * delta)
))

## Helper to draw split triangles ----------------------------------------------
make_triangles <- function(df, which = c("upper", "lower"), fill_col) {
  which <- match.arg(which)
  
  if (which == "upper") {
    ## upper-right triangle
    x_off  <- c(-0.5,  0.5,  0.5)
    y_off  <- c( 0.5,  0.5, -0.5)
    prefix <- "u"
  } else {
    ## lower-left triangle
    x_off  <- c(-0.5, -0.5,  0.5)
    y_off  <- c( 0.5, -0.5, -0.5)
    prefix <- "l"
  }
  
  data.frame(
    group = rep(paste0(prefix, seq_len(nrow(df))), each = 3),
    x     = rep(df$stop_look, each = 3) + rep(x_off, times = nrow(df)),
    y     = rep(df$y,         each = 3) + rep(y_off, times = nrow(df)),
    fill  = rep(df[[fill_col]], each = 3)
  )
}

upper_tri <- make_triangles(df_plot, "upper", "eff_fill")
lower_tri <- make_triangles(df_plot, "lower", "fut_fill")

## Labels ----------------------------------------------------------------------
## These show all values, including 0.0%
eff_text <- df_plot %>%
  transmute(
    x = stop_look + 0.22,
    y = y + 0.22,
    label = percent(eff_prop, accuracy = 0.1)
  )

fut_text <- df_plot %>%
  transmute(
    x = stop_look - 0.22,
    y = y - 0.22,
    label = percent(fut_prop, accuracy = 0.1)
  )

## Dummy data for legend -------------------------------------------------------
legend_df <- data.frame(
  x    = c(1, 1),
  y    = c(1, 1),
  type = c("Efficacy", "Futility")
)

## Plot ------------------------------------------------------------------------
efffut_split <- ggplot() +
  geom_polygon(
    data = upper_tri,
    aes(x = x, y = y, group = group, fill = fill),
    color = "white",
    linewidth = 0.4
  ) +
  geom_polygon(
    data = lower_tri,
    aes(x = x, y = y, group = group, fill = fill),
    color = "white",
    linewidth = 0.4
  ) +
  geom_text(
    data = eff_text,
    aes(x = x, y = y, label = label),
    size = 4.8,
    fontface = "bold"
  ) +
  geom_text(
    data = fut_text,
    aes(x = x, y = y, label = label),
    size = 4.8,
    fontface = "bold"
  ) +
  geom_point(
    data = legend_df,
    aes(x = x, y = y, color = type),
    alpha = 0,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  scale_fill_identity() +
  scale_color_manual(
    name = "Stopping reason",
    values = c(
      "Efficacy" = "green4",
      "Futility" = "red3"
    ),
    guide = guide_legend(
      override.aes = list(
        alpha = 1,
        shape = 15,
        size  = 8
      )
    )
  ) +
  scale_x_continuous(
    breaks = seq_len(K),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(theta_levels),
    labels = y_labs,
    expand = c(0, 0)
  ) +
  coord_fixed(ratio = 1) +
  labs(
    x = "Analysis (k)",
    y = expression(theta)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "gray98", color = NA),
    axis.title.x  = element_text(size = 16),
    axis.title.y  = element_text(size = 16),
    axis.text.x   = element_text(size = 12),
    axis.text.y   = element_text(size = 14),
    legend.title  = element_text(size = 16),
    legend.text   = element_text(size = 14),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

efffut_split

## Save output -----------------------------------------------------------------
ggsave(
  filename = "output/bayesian/figures/tabthree/theta_prior/efffut_split.png",
  plot = efffut_split,
  width = 9,
  height = 9.6,
  units = "in",
  dpi = 800
)

## Optional: vector version for publication quality ----------------------------
# ggsave(
#   filename = "output/bayesian/figures/tabthree/theta_prior/efffut_split.pdf",
#   plot = efffut_split,
#   width = 9,
#   height = 9.6,
#   units = "in",
#   device = cairo_pdf
# )