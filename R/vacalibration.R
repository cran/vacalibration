#' VA-calibration function
#'
#' @param va_data A named list. Algorithm-specific unlabeled VA-only data.
#'
#' For example, \code{list("algo1" = algo1_output, "algo2" = algo2_output, ...)}.
#'
#' Algorithm names (\code{"algo1"}, \code{"algo2"}, ...) can be "eava", "insilicova", or "interva".
#'
#' Data (\code{algo1_output}, \code{algo2_output}, ...) can be specific causes (output from \code{codEAVA()} function in \code{EAVA} and \code{crossVA()} function in \code{openVA}), or broad causes (output from the \code{cause_map()} function in this package), or broad-cause-specific death counts (integer vector).
#'
#' Can be different for different algorithms.
#'
#' Total number of deaths for different algorithms can be different.
#'
#'
#' @param age_group Character. Age-group of interest.
#'
#' \code{"neonate"} or \code{"child"}.
#'
#' \code{"neonate"} ages between 0-27 days, or \code{"child"} ages between 1-59 months.
#'
#'
#' @param country Character. The country \code{va_data} is from.
#'
#' Country-specific calibration is possible for "Bangladesh", "Ethiopia", "Kenya", "Mali", "Mozambique", "Sierra Leone", "South Africa".
#'
#' Any other country is matched with "other".
#'
#'
#' @param calibmodel.type Character. How to utilize misclassification estimates.
#'
#' \code{"Mmatprior"} (default). Propagates uncertainty in the misclassification matrix estimates.
#'
#' \code{"Mmatfixed"}. Uses fixed (default: posterior mean) misclassification matrix estimates.
#'
#'
#' @param Mmat.asDirich A named list. Similarly structured as \code{va_data}.
#'
#' Needed only if \code{calibmodel.type = "Mmatprior"} (propagates uncertainty).
#'
#' For example, \code{list("algo1" = Mmat.asDirich_algo1, "algo2" = Mmat.asDirich_algo2, ...)}.
#'
#' List of algorithm-specific Dirichlet prior on misclassification matrix to be used for calibration.
#'
#' Names and length must be identical to \code{va_data}.
#'
#' If algorithm names (\code{"algo1"}, \code{"algo2"}, ...) are \code{"eava"}, \code{"insilicova"} or \code{"interva"}, and \code{Mmat.asDirich} is missing, it by default uses the CHAMPS-based estimates (Dirichlet approximation of posterior) stored in \code{Mmat_champs} in this package.
#'
#' See \code{Mmat_champs} for details.
#'
#' If \code{Mmat.asDirich} is not missing, whatever provided is used.
#'
#' If any algorithm name (\code{"algo1"}, \code{"algo2"}, ...) is different from \code{"eava"}, \code{"insilicova"} or \code{"interva"}, \code{Mmat.asDirich} must be provided.
#'
#' \code{Mmat.asDirich_algo1} is a matrix of dimension CHAMPS ("gold standard") cause X VA cause.
#'
#' \code{Dirichlet(Mmat.asDirich_algo1[i,])} is used as informative prior on classification rates for CHAMPS cause \code{i}.
#'
#'
#' @param Mmat.fixed A named list. Similarly structured as \code{va_data} or \code{Mmat.asDirich}.
#'
#' Needed only if \code{calibmodel.type = "Mmatfixed"} (no uncertainty propagation).
#'
#' For example, \code{list("algo1" = Mmat.fixed_algo1, "algo2" = Mmat.fixed_algo2, ...)}
#'
#' List of algorithm-specific fixed misclassification matrix to be used for calibration.
#'
#' Names and length must be identical to \code{va_data}.
#'
#' If algorithm names (\code{"algo1"}, \code{"algo2"}, ...) are \code{"eava"}, \code{"insilicova"}, or \code{"interva"} and \code{Mmat.fixed} is missing, it by default uses the CHAMPS-based estimates (posterior mean) stored in \code{Mmat_champs} in this package.
#'
#' See \code{Mmat_champs} for details.
#'
#' If \code{Mmat.fixed} is not missing, whatever provided is used.
#'
#' If any algorithm name (\code{"algo1"}, \code{"algo2"}, ...) is different from \code{"eava"}, \code{"insilicova"} or \code{"interva"}, \code{Mmat.fixed} must be provided. \code{Mmat.fixed_algo1} is a matrix of dimension CHAMPS cause X VA cause. \code{Mmat.fixed_algo1[i,]} are the classification rates for CHAMPS cause \code{i}.
#'
#'
#' @param donotcalib A named list. Similarly structured as \code{va_data}, \code{Mmat.asDirich}, or \code{Mmat.fixed}.
#'
#' List of broad causes for each CCVA algorithm that we do not want to calibrate
#'
#' Default: \code{list("eava"="other", "insilicova"="other", "interva"="other")}. That is, \code{"other"} cause is not calibrated.
#'
#' For neonates, the broad causes are \code{"congenital_malformation"}, \code{"pneumonia"}, \code{"sepsis_meningitis_inf"}, \code{"ipre"}, \code{"other"}, or \code{"prematurity"}.
#'
#' For children, the broad causes are \code{"malaria"}, \code{"pneumonia"}, \code{"diarrhea"}, \code{"severe_malnutrition"}, \code{"hiv"}, \code{"injury"}, \code{"other"}, \code{"other_infections"}, \code{"nn_causes"} (neonatal causes).
#'
#' Set \code{list("eava" = NULL, "insilicova" = NULL, "interva" = NULL)} if you want to calibrate all causes.
#'
#'
#' @param donot.calib_type Character. \code{"fixed"} or \code{"learn"} (default).
#'
#' For \code{"fixed"}, only broad causes that are provided in \code{"donotcalib"} are not calibrated.
#'
#' For \code{"learn"}, it learns from \code{"Mmat.fixed"} or \code{"Mmat.asDirich"} if any other causes cannot be calibrated.
#'
#' For \code{"learn"}, it identifies VA causes for which the misclassification rates do not vary across CHAMPS causes.
#'
#' In that case, the calibration equation becomes ill-conditioned (see the footnote below Section 3.8 in Pramanik et al. (2025)). Currently, we address this by not calibrating VA causes for which the misclassification rates are similar along the rows (CHAMPS causes). VA causes (Columns) for which the rates along the rows (CHAMPS causes) do not vary more that \code{"nocalib.threshold"} are not calibrated. \code{"donotcalib"} is accordingly updated for each CCVA algorithm.
#'
#'
#' @param nocalib.threshold Numeric between 0 and 1. The value used for screening VA causes that cannot be calibrated when \code{donot.calib_type = "learn"}. Default: 0.1.
#'
#'
#' @param stable Logical. \code{TRUE} (default) or \code{FALSE}. Setting \code{TRUE} improves stability in calibration.
#'
#'
#' @param ensemble Logical. \code{TRUE} (default) or \code{FALSE}.
#'
#' Whether to perform ensemble calibration when outputs from multiple algorithms are provided.
#'
#'
#' @param pss Positive numeric. Degree of shrinkage of calibrated cause-specific mortality fraction (CSMF) estimate towards uncalibrated estimates.
#'
#' Always 0 when \code{stable=TRUE}. Defaults to 4 when \code{stable=FALSE}.
#'
#'
#' @param nMCMC Positive integer. Total number of posterior samples to perform inference on.
#'
#' Total number of iterations are \code{nBurn + nMCMC*nThin}.
#'
#' Default 5000.
#'
#' @param nBurn Positive integer. Total burn-in in posterior sampling.
#'
#' Total number of iterations are \code{nBurn + nMCMC*nThin}.
#'
#' Default 5000.
#'
#' @param nThin Positive integer. Number of thinning in posterior sampling.
#'
#' Total number of iterations are \code{nBurn + nMCMC*nThin}.
#'
#' Default 1.
#'
#'
#' @param adapt_delta_stan Positive numeric between 0 and 1. \code{"adapt_delta"} parameter in \code{rstan}.
#'
#' Influences the behavior of the No-U-Turn Sampler (NUTS), the primary MCMC sampling algorithm in Stan.
#'
#' Default 0.9.
#'
#'
#' @param refresh.stan Positive integer. Report progress at every \code{refresh.stan}-th iteration.
#'
#' Default \code{(nBurn + nMCMC*nThin)/10}, that is at every 10% progress.
#'
#'
#' @param seed Numeric. \code{"seed"} parameter in rstan.
#'
#' Default 1.
#'
#'
#' @param verbose Logical. Reports progress or not.
#'
#' \code{TRUE} (default) or \code{FALSE}.
#'
#'
#' @param saveoutput Logical. Save output or not.
#'
#' \code{TRUE} (default) or \code{FALSE}.
#'
#'
#' @param output_filename Character. Output name to save as.
#'
#' Default \code{paste0("calibratedva_", calibmodel.type)}. That is \code{"calibratedva_Mmatprior"} or \code{"calibratedva_Mmatfixed"}.
#'
#'
#' @param plot_it Logical. Whether to return comparison plot for summary.
#'
#' \code{TRUE} (default) or \code{FALSE}.
#'
#' @return A named list:
#'  \describe{
#'      \item{input}{A named list of input data}
#'      \item{p_uncalib}{Uncalibrated cause-specific mortality fractions (CSMF) estimates as observed in the data}
#'      \item{p_calib}{Posterior samples of calibrated CSMF estimates}
#'      \item{pcalib_postsumm}{Posterior summaries (mean and 95% credible interval) of calibrated CSMF estimates}
#'      \item{va_deaths_uncalib}{Uncalibrated cause-specific death counts as observed in the data}
#'      \item{va_deaths_calib_algo}{Algorithm-specific calibrated cause-specific death counts}
#'      \item{va_deaths_calib_ensemble}{Ensemble calibrated cause-specific death counts}
#'      \item{donotcalib}{A logical indicator of causes that are not calibrated for each algorithm}
#'      \item{causes_notcalibrated}{Causes that are not calibrated for each algorithm}
#'  }
#'
#' @examples
#'
#' \donttest{
#'
#' ######### VA input as specific causes #########
#' # output from codEAVA() function in the EAVA package and crossVA() function in openVA package
#'
#' # COMSA-Mozambique: Example (Publicly Available Version)
#' # Individual-Level Specific (High-Resolution) Cause of Death Data
#' data(comsamoz_public_openVAout)
#' head(comsamoz_public_openVAout$data)  # head of the data
#' comsamoz_public_openVAout$data[1,]  # ID and specific cause of death for individual 1
#'
#' # VA-calibration for the "neonate" age group and InSilicoVA algorithm
#' calib_out_specific = vacalibration(va_data =
#'                                             setNames(list(comsamoz_public_openVAout$data),
#'                                                      list(comsamoz_public_openVAout$va_algo)),
#'                                      age_group = comsamoz_public_openVAout$age_group,
#'                                      country = "Mozambique")
#'
#' ### comparing uncalibrated CSMF estimates and posterior summary of calibrated CSMF estimates
#' calib_out_specific$p_uncalib # uncalibrated
#' calib_out_specific$pcalib_postsumm["insilicova",,]
#'
#' ######### VA input as broad causes (output from cause_map()) #########
#'
#' # COMSA-Mozambique: Example (Publicly Available Version)
#' # Individual-Level Broad Cause of Death Data
#' data(comsamoz_public_broad)
#' head(comsamoz_public_broad$data)
#' comsamoz_public_broad$data[1,]  # binary vector indicating cause of death for individual 1
#'
#' # VA-calibration for the "neonate" age group and InSilicoVA algorithm
#' calib_out_broad = vacalibration(va_data = setNames(list(comsamoz_public_broad$data),
#'                                                      list(comsamoz_public_broad$va_algo)),
#'                                   age_group = comsamoz_public_broad$age_group,
#'                                   country = "Mozambique")
#'
#' ### comparing uncalibrated CSMF estimates and posterior summary of calibrated CSMF estimates
#' calib_out_broad$p_uncalib # uncalibrated
#' calib_out_broad$pcalib_postsumm["insilicova",,]
#'
#' ######### VA input as national death counts for different broad causes #########
#' calib_out_asdeathcount = vacalibration(va_data =
#'                                            setNames(list(colSums(comsamoz_public_broad$data)),
#'                                                     list(comsamoz_public_broad$va_algo)),
#'                                          age_group = comsamoz_public_broad$age_group,
#'                                          country = "Mozambique")
#'
#' ### comparing uncalibrated CSMF estimates and posterior summary of calibrated CSMF estimates
#' calib_out_asdeathcount$p_uncalib # uncalibrated
#' calib_out_asdeathcount$pcalib_postsumm["insilicova",,]
#'
#'
#' ######### Example of data based on EAVA and InSilicoVA for neonates in Mozambique #########
#' ## example VA national death count data from EAVA and InSilicoVA
#' va_data_example = list("eava" = c("congenital_malformation" = 40, "pneumonia" = 175,
#'                                   "sepsis_meningitis_inf" = 265, "ipre" = 220,
#'                                   "other" = 30, "prematurity" = 170),
#'                        "insilicova" = c("congenital_malformation" = 5, "pneumonia" = 145,
#'                                         "sepsis_meningitis_inf" = 370, "ipre" = 330,
#'                                         "other" = 60, "prematurity" = 290))
#'
#' ## algorithm-specific and ensemble calibration of EAVA and InSilicoVA
#' calib_out_ensemble = vacalibration(va_data = va_data_example,
#'                                    age_group = "neonate", country = "Mozambique")
#'
#' ### comparing uncalibrated CSMF estimates and posterior summary of calibrated CSMF estimates
#' calib_out_ensemble$p_uncalib # uncalibrated
#' calib_out_ensemble$pcalib_postsumm["eava",,] # EAVA-specific calibration
#' calib_out_ensemble$pcalib_postsumm["insilicova",,] # InSilicoVA-specific calibration
#' calib_out_ensemble$pcalib_postsumm["ensemble",,] # Ensemble calibration
#'
#' }
#'
#' @export
vacalibration <- function(va_data = NULL, age_group = NULL, country = NULL,
                          calibmodel.type = c("Mmatprior", "Mmatfixed")[1],
                          Mmat.asDirich = NULL, Mmat.fixed = NULL,
                          donotcalib = NULL,
                          donot.calib_type = c("learn", "fixed")[1], nocalib.threshold = 0.1,
                          stable = TRUE, ensemble = NULL,
                          pss = NULL,
                          nMCMC = 5000, nBurn = 5000, nThin = 1,
                          adapt_delta_stan = .9, refresh.stan = NULL,
                          seed = 1, verbose = TRUE, saveoutput = FALSE,
                          output_filename = NULL,
                          plot_it = TRUE){


  input.list = as.list(environment())

  if(is.null(va_data)){

    message("Provide unlabeled VA-only data in 'va_data' that you want to calibrate.")
    message("For example, 'va_data'=list('eava'=eava_output, 'insilicova'=insilicova_output, 'interva'=interva_output).")
    message("")
    message("Note:")
    message("     1. eava_output for the eava algorithm are the output from the codEAVA() function in EAVA package.")
    message("     2. insilicova_output and interva_output for the insilicova and inteva algorithms are outputs from the openVA package.")
    message("")
    return(NULL)

  }

  if(is.null(age_group)) stop("Provide age group. Can be 'neonate' (for 0-27 days) or 'child' (for 1-59 months)")

  if(is.null(country)) stop("Provide the country the 'va_data' is from")

  if((calibmodel.type=="Mmatprior")&(is.null(Mmat.asDirich))){

    utils::data("Mmat_champs", package = "vacalibration", envir = environment())
    helper_Mmat = Mmat_champs

    Mmat.asDirich = lapply(1:length(va_data),
                           FUN = function(k){

                             helper_Mmat[[age_group]][[names(va_data)[k]]]$asDirich[[ifelse(country %in% names(helper_Mmat[[age_group]][[names(va_data)[k]]]$asDirich), country, "other")]]

                           })
    names(Mmat.asDirich) = names(va_data)
    # print(Mmat.asDirich)
    # print(helper_Mmat[[age_group]][[names(va_data)[k]]]$asDirich[[ifelse(country %in% helper_Mmat[[age_group]][[names(va_data)[k]]]$asDirich, country, "other")]])

  }

  if((calibmodel.type=="Mmatfixed")&(is.null(Mmat.fixed))){

    utils::data("Mmat_champs", package = "vacalibration", envir = environment())
    helper_Mmat = Mmat_champs

    Mmat.fixed = lapply(1:length(va_data),
                        FUN = function(k){

                          helper_Mmat[[age_group]][[names(va_data)[k]]]$postmean[[ifelse(country %in% names(helper_Mmat[[age_group]][[names(va_data)[k]]]$postmean), country, "other")]]

                        })
    names(Mmat.fixed) = names(va_data)

  }

  if(is.null(donotcalib)) donotcalib = 'other'


  # reducing va_data to list of death counts (sufficient for calibration)
  va_data_tomodel = va_data
  for(k in 1:length(va_data)){

    if(length(dim(va_data[[k]]))==0){

      va_data_tomodel[[k]] = va_data[[k]]

    }else if(length(dim(va_data[[k]]))==2){

      if((ncol(va_data[[k]])==2)&&(isTRUE(all.equal(sort(colnames(va_data[[k]])), sort(c("ID", "cause")))))){

        va_data_tomodel[[k]] = colSums(cause_map(df = va_data[[k]], age_group = age_group))

      }else{

        va_data_tomodel[[k]] = colSums(va_data[[k]])

      }

    }else{

      message("")
      message("'va_data' for each algorithm can be:")
      message("     1. Output from the codEAVA() function in the EAVA package for EAVA algorithm, or Output from the crossVA() function in the openVA package for InSilicoVA and InterVA algorithm (matrix: individual X 2)")
      message("     2. Output of cause_map() function in this package after broad cause mapping (matrix: individual X broad causes)")
      message("     3. A integer vector of death counts for each broad cause (column sum of output of cause_map() function in this package)")
      message("")

      stop(paste0("'va_data' input for algorithm ", names(va_data)[k], " does not meet the expected data structure"))

    }

  }

  ## run calibration ----
  if(verbose){

    # cat("\n")
    message("Calibrating ...")

  }

  modular_vacalib_output = modular.vacalib(va_unlabeled = va_data_tomodel, age_group = age_group,
                                           calibmodel.type = calibmodel.type,
                                           Mmat.asDirich = Mmat.asDirich, Mmat.fixed = Mmat.fixed,
                                           donotcalib = donotcalib,
                                           donot.calib_type = donot.calib_type, nocalib.threshold = nocalib.threshold,
                                           stable = stable, ensemble = ensemble,
                                           pss = pss,
                                           nMCMC = nMCMC, nBurn = nBurn, nThin = nThin,
                                           adapt_delta_stan = adapt_delta_stan, refresh.stan = refresh.stan,
                                           seed = seed, verbose = verbose, saveoutput = saveoutput,
                                           output_filename = output_filename,
                                           plot_it = plot_it)

  if(verbose){

    message("Done.")
    # cat("\n")

  }

  return(modular_vacalib_output)

}


