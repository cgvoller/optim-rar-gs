df <- as.data.frame.table(fgtb3_results$rejecttb3, responseName = "reject") %>%
  setNames(c("sim", "look", "theta_idx", "reject")) %>%
  mutate(
    sim       = as.integer(sim),
    look      = as.integer(look),
    theta     = theta[theta_idx],
    reject    = as.integer(reject)
  )

# Compute the cumulative “stop for efficacy” up to each look
df_cum <- df %>%
  arrange(theta, sim, look) %>%
  group_by(theta, sim) %>%
  mutate(cum_eff   = cummax(reject == 1),
         cum_futil = cummax(reject == -1)) %>%
  ungroup()

# Summarise cumulative Type I error & Power by look, theta
summary_cum <- df_cum %>%
  group_by(theta, look) %>%
  summarise(
    cum_prob_reject = mean(cum_eff),
    se = sqrt((cum_prob_reject * (1 - cum_prob_reject)) / n()),
    .groups = "drop"
  ) %>%
  mutate(metric = if_else(theta == 0, "Type I error", "Power"))

summary_cum2 <- summary_cum %>%
  gather(key = "theta", value = "theta_value") %>%
  mutate(theta = as.factor(theta))

theta_map <- setNames(theta_fractions,
                      summary_cum %>% 
                        distinct(theta) %>% 
                        pull(theta))

summary_cum2 <- summary_cum %>%
  mutate(Theta_Fraction = theta_map[as.character(theta)])

summary_cum2$Theta_Fraction <- factor(summary_cum2$Theta_Fraction, levels=c(theta_fractions))

final_summary <- summary_cum2 %>%
  group_by(Theta_Fraction, metric) %>%
  filter(look == max(look)) %>%
  mutate(cum_prob_reject=paste(cum_prob_reject,"(",round2(se,5),")")) %>% 
  ungroup()

# Create the gt table
tab_group_oc <- final_summary %>%
  select(-c(theta,look,se)) %>% 
  gt(rowname_col = "Theta_Fraction") %>%
  tab_header(
    title = md("**Cumulative Rejection Probabilities: Type I Error and Power**"),
    subtitle = md("Final analysis (K=5) Results by $\\theta$")
  ) %>%
  fmt_number(
    columns = cum_prob_reject,
    decimals = 5
  ) %>%
  cols_label(
    metric = "Metric",
    cum_prob_reject = "Cumulative Probability"
  ) %>%
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    row_group.font.weight = "bold"
  )

save_outputs(table = tab_group_oc,
             name = "tab_group_oc",
             dir =  "output/frequentist/figures/tabthree")


