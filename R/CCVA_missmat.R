#' CCVA Misclassification Matrix Inventory
#'
#' This is the inventory of misclassification matrix estimates for \href{https://www.jogh.org/documents/issue201601/jogh-06-010601.pdf}{EAVA}, \href{https://www.tandfonline.com/doi/abs/10.1080/01621459.2016.1152191}{InSilicoVA}, and InterVA (\doi{10.3402/gha.v5i0.19281}) algorithms.
#' The estimates are derived using the misclassification matrix modeling framework from \href{https://projecteuclid.org/journals/annals-of-applied-statistics/volume-19/issue-2/Modeling-structure-and-country-specific-heterogeneity-in-misclassification-matrices-of/10.1214/24-AOAS2006.short}{Pramanik et al. (2025)}.
#' and paired CHAMPS–VA cause-of-death data from the Child Health and Mortality Prevention Surveillance (\href{https://champshealth.org/}{CHAMPS}) project.
#' Please refer to Pramanik et al. (2026; \doi{10.1136/bmjgh-2025-021747}) for details on analysis.
#' The package interpret CHAMPS and VA causes as true and estimated causes.
#'
#' @format Nested list.
#' \describe{
#'   \item{age_group}{\code{"neonate"} for 0-27 days, and \code{"child"} for 1-59 months}
#'   \item{va_algo}{\code{"eava"}, \code{"insilicova"}, and \code{"interva"}}
#'   \item{estimate types}{\code{"postsumm"} contains posterior summaries, \code{"postmean"} contains the posterior means, and \code{"asDirich"} contains Dirichlet approximation for each CHAMPS cause and country.}
#'   \item{country}{Seven specific countries: \code{"Bangladesh"}, \code{"Ethiopia"}, \code{"Kenya"}, \code{"Mali"}, \code{"Mozambique"}, \code{"Sierra Leone"}, and \code{"South Africa"}. For all other countries, use \code{"other"}.}
#'   \item{version}{Character. Date stamp (yyyymmdd) for version control Only for package maintainers.}
#' }
#'
#' @details
#' Format: \code{CCVA_missmat[[age_group]][[va_algo]][[estimate types]][[country]]}.
#'
#' \code{CCVA_missmat[[age_group]][[va_algo]][["postsumm"]][[country]]} contains posterior summaries of misclassification matrices for a given \code{age_group}, \code{va_algo}, and \code{country}.
#' It is an array arranged as the number of posterior summaries \out{×} CHAMPS cause \out{×} VA cause.
#'
#' Neonatal causes include \code{"congenital_malformation"}, \code{"pneumonia"}, \code{"sepsis_meningitis_inf"}, \code{"ipre"}, \code{"other"}, and \code{"prematurity"}.
#'
#' Child causes encompass \code{"malaria"}, \code{"pneumonia"}, \code{"diarrhea"}, \code{"severe_malnutrition"}, \code{"hiv"}, \code{"injury"}, \code{"other"}, \code{"other_infections"}, and \code{"nn_causes"}.
#'
#' For example, for \code{"neonate"} age group, \code{"eava"} algorithm in \code{"Mozambique"},
#' \itemize{
#'    \item \code{CCVA_missmat$neonate$eava$postsumm$Mozambique[,"pneumonia","pneumonia"]} are posterior summaries of the sensitivity for "pneumonia".
#'    \item \code{CCVA_missmat$neonate$eava$postsumm$Mozambique[,"pneumonia","ipre"]} are posterior summaries of the false negative rate for CHAMPS cause "pneumonia" and VA cause "ipre".
#' }
#'
#' \code{CCVA_missmat[[age_group]][[va_algo]][["postmean"]][[country]]} contains posterior means of misclassification matrices for a given \code{age_group}, \code{va_algo}, and \code{country}.
#' It is a matrix arranged as CHAMPS cause \out{×} VA cause.
#'
#' For example, for \code{"neonate"} age group, \code{"eava"} algorithm in \code{"Mozambique"},
#' \itemize{
#'    \item \code{CCVA_missmat$neonate$eava$postmean$Mozambique["pneumonia","pneumonia"]} is the posterior mean of the sensitivity for "pneumonia".
#'    \item \code{CCVA_missmat$neonate$eava$postmean$Mozambique["pneumonia","ipre"]} is the posterior mean of the false negative rate for CHAMPS cause "pneumonia" and VA cause "ipre".
#' }
#'
#' \code{CCVA_missmat[[age_group]][[va_algo]][["asDirich"]][[country]]} contains Dirichlet approximations of misclassification matrices for a given \code{age_group}, \code{va_algo}, and \code{country}.
#' It is a matrix arranged as CHAMPS cause \out{×} VA cause.
#' Each row contains Dirichlet scale parameters that best approximates the marginal posterior of misclassification for each CHAMPS cause (rows), \code{age_group}, \code{va_algo}, and \code{country}.
#'
#' For example, for \code{"neonate"} age group, \code{"eava"} algorithm in \code{"Mozambique"},
#' the Dirichlet distribution with scale parameters \code{CCVA_missmat$neonate$eava$asDirich$Mozambique["pneumonia",]} best approximates the marginal posterior of misclassification rates for CHAMPS cause \code{"pneumonia"}.
#'
#' Specific estimates are available for seven countries: \code{"Bangladesh"}, \code{"Ethiopia"}, \code{"Kenya"}, \code{"Mali"}, \code{"Mozambique"}, \code{"Sierra Leone"}, and \code{"South Africa"}. For all other countries, the package uses the estimate for \code{"other"}. This estimate is centered at the misclassification matrix pooled across countries, and its uncertainty reflects the degree of cross-country heterogeneity observed across the seven CHAMPS countries.
#'
#' Due to file size limit, the posterior samples corresponding to this inventory are available at \href{https://github.com/sandy-pramanik/CCVA-Misclassification-Matrices}{CCVA-Misclassification-Matrices} GitHub repository.
#'
#' For example, \code{CCVA_missmat$neonate$eava$postsamples$Mozambique} contains misclassification matrix samples for \code{eava} among \code{neonate} in \code{Mozambique}.
#'
#' The .rda file is available under the \href{https://github.com/sandy-pramanik/CCVA-Misclassification-Matrices/releases/tag/20241004}{release}.
#'
#' @references
#' Pramanik, S, et al. (2026)
#' Country-Specific Estimates of Misclassification Rates of Computer-Coded Verbal Autopsy Algorithms
#' BMJ Global Health
#' \doi{10.1136/bmjgh-2025-021747}
#'
#' Pramanik, S, et al. (2025)
#' Modeling structure and country-specific heterogeneity in misclassification matrices of verbal autopsy-based cause of death classifiers
#' Annals of Applied Statistics
#' \href{https://projecteuclid.org/journals/annals-of-applied-statistics/volume-19/issue-2/Modeling-structure-and-country-specific-heterogeneity-in-misclassification-matrices-of/10.1214/24-AOAS2006.short}{Link}
#'
#' Wilson E, et al. (2025)
#' EAVA: Deterministic Verbal Autopsy Coding with Expert Algorithm Verbal Autopsy
#' \href{https://CRAN.R-project.org/package=EAVA}{Link}
#'
#' Zehang Richard Li, et al. (2024)
#' openVA: Automated Method for Verbal Autopsy
#' R package version 1.1.2.
#' \href{https://CRAN.R-project.org/package=openVA}{Link}
#'
#' Zehang Richard Li, et al. (2023)
#' The openVA Toolkit for Verbal Autopsies
#' The R Journal \href{https://journal.r-project.org/articles/RJ-2023-020/}{Link}
#'
#' Kalter, H., et al. (2016)
#' Validating hierarchical verbal autopsy expert algorithms in a large data set with known causes of death.
#' J Glob Health
#' \href{https://www.jogh.org/documents/issue201601/jogh-06-010601.pdf}{Link}
#'
#' McCormick, Tyler H., et al. (2016)
#' Probabilistic Cause-of-Death Assignment Using Verbal Autopsies
#' Journal of the American Statistical Association
#' \href{https://www.tandfonline.com/doi/abs/10.1080/01621459.2016.1152191}{Link}
#'
#' Byass, Peter, et al. (2012)
#' Strengthening standardised interpretation of verbal autopsy data: the new InterVA-4 tool
#' Global Health Action
#' \doi{10.3402/gha.v5i0.19281}
#'
"CCVA_missmat"
