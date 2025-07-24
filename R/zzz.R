# zzz.R - Initialization code for vacalibration package

# Declare global variables to avoid R CMD check NOTES
# These variables are used in non-standard evaluation (e.g., data.table, dplyr, etc.)

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "cause",
    "Var2",
    "Var3",
    "value",
    "calib_type",
    "llim",
    "ulim",
    "Mmat_champs"
  ))
}
