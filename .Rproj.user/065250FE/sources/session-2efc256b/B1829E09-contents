## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/BRAR
##
## Purpose of script: Function to create table output
##
## Author: Corey Voller
##
## Date Created: 21-02-2025
##
## QC'd by:
## QC date:
##
## ─────────────────────────────────────────────────────────────────────────────
##
## Notes:
##   
##
## ─────────────────────────────────────────────────────────────────────────────
##
## 
## Create tab ------------------------------------------------------------------

create_tab_fun <- function(col1,col2,col3,col4,col5,col6,col7,col8,col9,headertxt,footnotetxt){
  
# Create delta symbols
theta_fractions = paste(MASS::fractions(theta/delta), "\U1D6FF", sep = "")

# Create table with thetas using pooled theta
rarpooled <- data.frame(theta_fractions,
                        Design = c("RAR Using Pooled \U03B8"),
                        EN1 = col4,
                        EN2 = col5,
                        total=col6)
# Create table with thetas using group based theta
groupbased <- data.frame(theta_fractions,
                         Design = c("RAR Using Grouped \U03B8"),
                         EN1 = col7,
                         EN2 = col8,
                         total = col9)
# Create Theta and target columns
theta_tab <- data.frame(theta_fractions,
                        Design = c("Target"),
                        EN1 = col1,
                        EN2 = col2,
                        total = col3)
# Join data frames
tab_joined <- rbind(theta_tab, rarpooled,groupbased)

# Pivot
tab_wide <- tab_joined %>%
  pivot_wider(
    names_from = Design,
    values_from = c(EN1, EN2,total),
    names_glue = "{Design}_{.value}",
    names_vary = 'slowest'
  )
# Create pretty table
tab <- tab_wide %>%
  rename("\U03B8" = theta_fractions) %>%
  gt(tab_joined, rowname_col = "\U03B8") %>%
  tab_header(paste(headertxt)) %>%
  tab_stubhead(label = "\U1D703") %>%
  tab_spanner_delim(delim = "_") %>%
  #fmt_number(n_sigfig = 5) %>%
  tab_footnote(footnote = paste("Results are based on", nsims, "simulations;",footnotetxt)) %>%
  fmt_markdown(columns = everything()) %>%
  opt_table_font(
    font = list(gt::google_font("Lato")),
    weight = 400,
    size = px(11)  # Reduce font size
  )
return(tab)
}

## Create bayes/freq tab -------------------------------------------------------

create_Bayes_ESS <- function(comparison_df,
                                    theta_fractions,
                                    comparison_cols,
                                    col_labels,
                                    header_text,
                                    nsims,
                                    dig,
                                    priors_info) {
  comparison_df %>%
    mutate(theta = theta_fractions) %>%
    gt(rowname_col = "theta") %>%
    fmt_icon(
      columns = all_of(comparison_cols),
      fill_color = c(
        "arrow-up" = "green",
        "arrow-down" = "red",
        "equals" = "blue"
      )
    ) %>%
    tab_header(title = paste(header_text)) %>%
    cols_label(!!!col_labels) %>%
    tab_footnote(
      footnote = paste0(
        "Based on ",
        nsims,
        " Simulations; ",
        "mean (standard error) to ",
        dig,
        " decimal place"
      )
    ) %>%
    tab_footnote(footnote = priors_info) %>%
    opt_table_font(
      font = list(gt::google_font("Lato")),
      weight = 400,
      size = px(10)
    )
}

