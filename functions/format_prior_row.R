## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/BRAR_gsd
##
## Purpose of script: Format Rows to output priors in footnotes
##
## Author: Corey Voller
##
## Date Created: 13-03-2025
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
## Format Row ------------------------------------------------------------------

# Take a dataframe of priors and print them in a string to be used in footnote
format_row <- function(row_name, row_data) {
  formatted_values <- sapply(row_data, function(val) {
    if (is.numeric(val)) {
      # Format numeric values to 1 decimal places
      return(as.character(val))
    } else {
      # If not numeric, just return the value as is (character or factor)
      return(as.character(val))
    }
  })
  
  paste(row_name, ":",
        paste(names(row_data), ": ", formatted_values, sep = "", collapse = ", "),
        sep = "")
}
