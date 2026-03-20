#' VA-Calibration
#'
#' This is the main function in the package. It calibrates population-level cause-specific mortality fractions (CSMFs) that are derived using computer-coded verbal autopsy (CCVA) algorithms. For VA-Calibration, the function utilizes the inventory of misclassification matrix estimates \link{CCVA_missmat}.
#' The outputs from \href{https://CRAN.R-project.org/package=EAVA}{EAVA} and \href{https://CRAN.R-project.org/package=openVA}{openVA} for InSilicoVA and InterVA can be input directly (see below). This seamlessly supports VA-Calibration for \href{https://www.jogh.org/documents/issue201601/jogh-06-010601.pdf}{EAVA}, \href{https://www.tandfonline.com/doi/abs/10.1080/01621459.2016.1152191}{InSilicoVA}, and InterVA (\doi{10.3402/gha.v5i0.19281}).
#' For other CCVA algorithms, the input expects either an individual by cause matrix, or cause-specific death count vector (see below).
#' When broad-cause-specific death counts are input and they do not match the broad causes in the stored misclassification estimates, then either \code{studycause_map} or the misclassification matrices (fixed or as row-specific Dirichlet priors) need to be provided.
#' More generally, this allows us to calibrate population-level prevalence derived from single-class predictions of discrete classifiers. For this, users need to provide fixed or uncertainty-quantified misclassification matrices.
#'
#' @param va_data Named list. Algorithm-specific unlabeled VA data.
#'
#' It expects a named list, such as \code{list("algo1" = algo1_output, "algo2" = algo2_output, ...)}.
#'
#' Misclassification matrix estimates in \code{CCVA_missmat} are only available for CCVA algorithms \href{https://www.jogh.org/documents/issue201601/jogh-06-010601.pdf}{EAVA}, \href{https://www.tandfonline.com/doi/abs/10.1080/01621459.2016.1152191}{InSilicoVA}, and InterVA (\doi{10.3402/gha.v5i0.19281}). For them the algorithm names in input data must be \code{"eava"}, \code{"insilicova"}, and \code{"interva"}. Otherwise, users must input misclassification matrices in \code{missmat} (see more details in \code{missmat}).
#'
#' VA data provided for each algorithm (\code{algo1_output}, \code{algo2_output}, ...) can be either
#'    1. outputs of CCVA algorithms (output from \code{codEAVA()} in \code{EAVA} for EAVA, and \code{codeVA()} and \code{prepCalibration()} in \code{openVA} for InSilicoVA and InterVA), or
#'    2. individual broad cause of deaths (output from \code{cause_map}), or
#'    3. a vector of cause-specific death counts.
#'
#' More generally, it can calibrate for any discrete classifier. In that case, the input must be one of these two types:
#'    1. A binary matrix arranged as individuals along rows and class labels as columns. For each individual (row), 1 occurs exactly once and it indicates the estimated class label. Other elements in the row are 0.
#'    2. A vector of label-specific counts. This indicates the estimated number of individuals for each label.
#'
#'
#' @param age_group Character.
#'
#' When \code{missmat} is \code{NULL}, this indicates the age group for which the misclassification matrix estimates in \code{"CCVA_missmat"} should be applied (default).
#'
#' It can be either \code{"neonate"} for neonatal deaths occurring between 0-27 days after birth, and \code{"child"} for deaths among children occurring between 1-59 months.
#'
#'
#' @param country Character.
#'
#' When \code{missmat} is \code{NULL}, this indicates the country for which the misclassification matrix estimates in \code{"CCVA_missmat"} should be applied (default).
#'
#' If input is "Bangladesh", "Ethiopia", "Kenya", "Mali", "Mozambique", "Sierra Leone", or "South Africa", then their corresponding misclassification matrix is applied. For any other country, the estimate for "other" is applied (see \code{"CCVA_missmat"} for more details).
#'
#'
#' @param missmat_type Character. Indicates the type of misclassification matrix estimates provided in \code{missmat}.
#'
#' \code{"prior"} (default) Dirichlet priors for each row of the misclassification matrix.
#'
#' \code{"fixed"} A fixed misclassification matrix.
#'
#' \code{"samples"} Random samples of misclassification matrix.
#'
#' Uncertainty in misclassification matrix estimates is only propagated for \code{"prior"} or \code{"samples"}.
#'
#' @param studycause_map Named character vector. A mapping of observed causes (in \code{va_data}) to broad causes (for which misclassification estimates are available in \code{"CCVA_missmat"}).
#'
#' Required only when \code{missmat} is \code{NULL}, and causes observed in \code{va_data} are not a subset of broad causes in \code{"CCVA_missmat"} (see \code{"CCVA_missmat"} for list of causes).
#'
#' For example, if causes observed in \code{va_data} for neonates are "cause1", "cause2", "cause3", and "cause4", \code{studycause_map} expects input as \code{c("cause1" = "pneumonia", "cause2" = "ipre", "cause3" = "other", "cause4" = "other")}.
#'
#' @param missmat Named list. Similarly structured as \code{va_data}. For example, \code{list("algo1" = missmat_algo1, "algo2" = missmat_algo2, ...)}.
#'
#' For \code{missmat_type = "prior"}, \code{missmat_algo1}, \code{missmat_algo2}, ... are matrices with positive entries and arranged as CHAMPS cause \out{×} VA cause. Each row of the matrix is a vector of Dirichlet scale parameters. This the Dirichlet prior assumed on the corresponding row of the misclassification matrix. See stored estimates \code{CCVA_missmat$neonate$eava$asDirich$Mozambique} for example.
#'
#' For \code{missmat_type = "fixed"}, \code{missmat_algo1}, \code{missmat_algo2}, ... are misclassification matrices arranged as CHAMPS cause \out{×} VA cause. See stored estimates \code{CCVA_missmat$neonate$eava$postmean$Mozambique} for example.
#'
#' For \code{missmat_type = "samples"}, \code{missmat_algo1}, \code{missmat_algo2}, ... are arrays of misclassification matrix samples arranged as samples \out{×} CHAMPS cause \out{×} VA cause. \code{missmat_algo1[i,,]} is the i-th sample of misclassification matrix for \code{algo1}. See the samples stored in the \href{https://github.com/sandy-pramanik/CCVA-Misclassification-Matrices}{CCVA-Misclassification-Matrices} GitHub repository for example.
#'
#' Names and length of \code{missmat} must be identical to \code{va_data}.
#'
#' Users are not required to provide \code{missmat} for using the stored estimates in \code{"CCVA_missmat"}. They can simply input the required \code{age_group}, \code{country}, \code{missmat_type}, and \code{studycause_map} accordingly.
#'
#' \code{missmat} needs to be input when causes observed in \code{va_data} are not a subset of CHAMPS broad causes (in \code{"CCVA_missmat"}) and \code{studycause_map} is not provided.
#'
#' For a general purpose of calibrating categorical classifiers, CHAMPS and VA causes can be interpreted as true and estimated labels and users must input \code{missmat}.
#'
#' @param donotcalib Named list. List of causes for each algorithm that users do not want to calibrate. The set of causes can differ across algorithms.
#'
#' Default: \code{list("eava"="other", "insilicova"="other", "interva"="other")}. When using the stored estimates in \code{CCVA_missmat}, this implies that the cause-specific mortality fractions (CSMF) for CHAMPS broad cause \code{"other"} is not calibrated.
#'
#' When causes observed in \code{va_data} are not a subset of CHAMPS broad causes and \code{studycause_map} is provided, all observed causes in \code{va_data} that match with the causes in \code{donotcalib} are not calibrated.
#'
#' Set \code{list("eava"=NULL, "insilicova"=NULL, "interva"=NULL)} to calibrate all causes.
#'
#' For a general purpose of calibrating categorical classifiers, causes can be interpreted as class labels and specified accordingly.
#'
#' @param donotcalib_type Character. \code{"learn"} (default) or \code{"fixed"}.
#'
#' For \code{donotcalib_type="fixed"}, only the causes specified in \code{"donotcalib"} are not calibrated.
#'
#' For \code{donotcalib_type="learn"}, it learns additional causes from misclassification matrix in \code{"missmat"} that cannot be calibrated.
#'
#' When misclassification rates for a VA cause do not change across CHAMPS causes, the calibration equation becomes underdetermined (see the footnote on pg. 1227 in \href{https://projecteuclid.org/journals/annals-of-applied-statistics/volume-19/issue-2/Modeling-structure-and-country-specific-heterogeneity-in-misclassification-matrices-of/10.1214/24-AOAS2006.full}{Pramanik et al. (2025)}). When \code{donotcalib_type="learn"}, it screens VA causes that do not vary beyond \code{nocalib.threshold}. These causes are added to the \code{donotcalib} list.
#'
#' For a general purpose of calibrating categorical classifiers, causes can be interpreted as class labels and specified accordingly.
#'
#' @param nocalib.threshold Numeric in \code{(0,1)}.
#'
#' The threshold used to screen VA causes when \code{donotcalib_type="learn"}.
#'
#' Default: 0.1.
#'
#' @param path_correction Logical. Setting \code{TRUE} shrinks misclassification matrix towards the identity matrix to improve stability in VA-Calibration.
#'
#' Default is \code{TRUE}.
#'
#' @param ensemble Logical. Whether to perform ensemble calibration when outputs from multiple algorithms are provided.
#'
#' Default is \code{TRUE}.
#'
#' @param pshrink_strength Positive numeric. Degree of shrinkage of calibrated CSMF estimates towards its uncalibrated estimates. This is the parameter \code{eta} in the prior of calibrated CSMF \code{p} (see pg. 1226 in \href{https://projecteuclid.org/journals/annals-of-applied-statistics/volume-19/issue-2/Modeling-structure-and-country-specific-heterogeneity-in-misclassification-matrices-of/10.1214/24-AOAS2006.full}{Pramanik et al. (2025)}).
#'
#' Only used when \code{path_correction=FALSE}. \code{pshrink_strength} is set to 0 when \code{path_correction=TRUE}.
#'
#' Defaults to 4 when \code{path_correction=FALSE}.
#'
#'
#' @param nMCMC Positive integer. Total number of posterior samples to perform inference on.
#'
#' Total number of iterations are \code{nBurn + nMCMC*nThin}. Default 5000.
#'
#'
#' @param nBurn Positive integer. Total burn-in in posterior sampling.
#'
#' Total number of iterations are \code{nBurn + nMCMC*nThin}. Default 5000.
#'
#'
#' @param nThin Positive integer. Number of thinning in posterior sampling.
#'
#' Total number of iterations are \code{nBurn + nMCMC*nThin}. Default 1.
#'
#'
#' @param nChain Positive integer. Number of chains for Stan sampling. Default 1.
#'
#'
#' @param nCore Positive integer. Number of cores to run multiple chains in parallel for Stan sampling. Default 1.
#'
#'
#' @param adapt_delta_stan Numeric in \code{(0,1)}. \code{adapt_delta} parameter in \code{rstan}.
#'
#' Influences the behavior of the No-U-Turn Sampler (NUTS) in Stan.
#'
#' Default 0.9.
#'
#'
#' @param refresh_stan Positive integer. Print every \code{refresh_stan}% progress.
#'
#' Default 20.
#'
#'
#' @param seed Numeric. \code{seed} parameter in rstan. Default 1.
#'
#'
#' @param verbose Logical. Whether to report progress (\code{TRUE}) or not (\code{FALSE}).
#'
#' Default \code{TRUE}.
#'
#'
#' @param saveoutput Logical. Save output (\code{TRUE}) or not (\code{FALSE}).
#'
#' Default \code{TRUE}.
#'
#'
#' @param output_filename Character. Output name to save as.
#'
#' Default \code{vacalibration_out}.
#'
#'
#' @param output_dir Output directory or file path to save at.
#'
#' Default \code{getwd()}, the working directory.
#'
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{calib_MCMCout} — Output from Stan fits.
#'   \item \code{p_uncalib} — Uncalibrated estimates of CSMF. It is a matrix arranged as algorithm \out{×} VA causes (estimated labels).
#'   \item \code{p_calib} — Posterior samples of calibrated CSMF. It is an array arranged as algorithm \out{×} samples \out{×} VA causes (or estimated labels).
#'   \item \code{pcalib_postsumm} — Posterior summaries (mean and 95% credible interval) of calibrated CSMF. It is an array arranged as algorithm \out{×} summary measures \out{×} VA causes (or estimated labels).
#'   \item \code{va_deaths_uncalib} — Uncalibrated cause-specific death counts. It is a matrix arranged as algorithm \out{×} VA causes (or estimated labels).
#'   \item \code{va_deaths_calib_algo} — Calibrated cause-specific death counts from algorithm-specific calibration. It is a matrix arranged as algorithm \out{×} VA causes (or estimated labels).
#'   \item \code{va_deaths_calib_ensemble} — Calibrated cause-specific death counts from ensemble calibration. It is a matrix arranged as algorithm \out{×} VA causes (or estimated labels).
#'   \item \code{Mmat_input} — \code{"missmat"} as provided in the input. It is an array arranged as algorithm \out{×} CHAMPS cause (or true labels) \out{×} VA causes (or estimated labels).
#'   \item \code{Mmat_study} — Modified \code{Mmat_input} if \code{studycause_map} is provided. It is an array arranged in the same way as \code{Mmat_input}.
#'   \item \code{Mmat_tomodel} — Modified \code{Mmat_study} if \code{path_correction} is \code{TRUE}. This is used for calibration. It is an array arranged in the same way as \code{Mmat_input} and \code{Mmat_study}.
#'   \item \code{donotcalib_study} — This indicates causes that are not calibrated for each algorithm, as specified in the input \code{donotcalib}. It is a logical matrix arranged as algorithm \out{×} VA causes (or estimated labels).
#'   \item \code{donotcalib_tomodel} — This indicates causes that are not calibrated in each calibration. This is a modified \code{donotcalib_study} if \code{donotcalib_type} is provided and \code{ensemble=TRUE}. It is a logical matrix arranged as algorithm \out{×} VA causes (or estimated labels).
#'   \item \code{calibrated} — \code{TRUE} or \code{FALSE} indicating whether Stan sampling was performed for calibration.
#'   \item \code{lambda_calibpath} — When \code{path_correction=TRUE}, this indicates the degree of shrinkage of CSMF for each algorithm towards uncalibrated estimates. This is a vector of numerics in \code{[0,1]} showing degrees of shrinkage for each algorithm.
#'   \item \code{K} — Number of algorithms.
#'   \item \code{nCause} — Number of causes.
#'   \item \code{causes} — Name of causes.
#'   \item \code{input} — List of inputs.
#' }
#'
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
#' Fiksel, J., et al. (2022)
#' Generalized Bayes Quantification Learning under Dataset Shift
#' Journal of the American Statistical Association
#' \href{https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1909599}{Link}
#'
#' Datta, A, et al. (2021)
#' Regularized Bayesian transfer learning for population-level etiological distributions.
#' Biostatistics
#' \doi{10.1093/biostatistics/kxaa001}
#'
#' @examples
#'
#' \donttest{
#' ######### COMSA-Mozambique VA-COD data #########
#' data(comsamoz_CCVAoutput)
#'
#' # neonatal deaths
#' comsamoz_CCVAoutput$neonate$eava  # output from running EAVA
#' comsamoz_CCVAoutput$neonate$insilicova  # output from running InSilicoVA
#' comsamoz_CCVAoutput$neonate$interva  # output from running InterVA
#'
#'
#'
#' ######### Algorithm-Specific Calibration #########
#'
#' # EAVA
#' vacalib_out_eava = vacalibration(va_data = comsamoz_CCVAoutput$neonate[1],
#'                                  age_group = "neonate", country = "Mozambique",
#'                                  saveoutput = FALSE)
#'
#' ## CSMF
#' vacalib_out_eava$p_uncalib   # uncalibrated
#' vacalib_out_eava$p_calib   # calibrated
#' vacalib_out_eava$pcalib_postsumm   # summary of calibrated estimates
#'
#' ## death counts
#' vacalib_out_eava$va_deaths_uncalib   # uncalibrated
#' vacalib_out_eava$va_deaths_calib_algo   # calibrated
#'
#'
#' # InSilicoVA
#' vacalib_out_insilicova = vacalibration(va_data = comsamoz_CCVAoutput$neonate[2],
#'                                        age_group = "neonate", country = "Mozambique",
#'                                        saveoutput = FALSE)
#'
#' ## CSMF
#' vacalib_out_insilicova$p_uncalib   # uncalibrated
#' vacalib_out_insilicova$p_calib   # calibrated
#' vacalib_out_insilicova$pcalib_postsumm   # summary of calibrated estimates
#'
#' ## death counts
#' vacalib_out_insilicova$va_deaths_uncalib   # uncalibrated
#' vacalib_out_insilicova$va_deaths_calib_algo   # calibrated
#'
#'
#' # InterVA
#' vacalib_out_interva = vacalibration(va_data = comsamoz_CCVAoutput$neonate[3],
#'                                     age_group = "neonate", country = "Mozambique",
#'                                     saveoutput = FALSE)
#'
#' ## CSMF
#' vacalib_out_interva$p_uncalib   # uncalibrated
#' vacalib_out_interva$p_calib   # calibrated
#' vacalib_out_interva$pcalib_postsumm   # summary of calibrated estimates
#'
#' ## death counts
#' vacalib_out_interva$va_deaths_uncalib   # uncalibrated
#' vacalib_out_interva$va_deaths_calib_algo   # calibrated
#'
#'
#'
#' ######### Ensemble Calibration #########
#' vacalib_out_ensemble = vacalibration(va_data = comsamoz_CCVAoutput$neonate,
#'                                      age_group = "neonate", country = "Mozambique",
#'                                      saveoutput = FALSE)
#'
#' ## CSMF
#' vacalib_out_ensemble$p_uncalib   # uncalibrated
#' vacalib_out_ensemble$p_calib   # calibrated
#' vacalib_out_ensemble$pcalib_postsumm   # summary of calibrated estimates
#'
#' ## death counts
#' vacalib_out_ensemble$va_deaths_uncalib   # uncalibrated
#' vacalib_out_ensemble$va_deaths_calib_algo   # algorithm-specific calibrated death counts
#' vacalib_out_ensemble$va_deaths_calib_ensemble   # ensemble calibrated death counts
#' }
#'
#' @importFrom stats nlminb
#' @importFrom LaplacesDemon ddirichlet
#'
#' @export
vacalibration <- function(va_data = NULL, age_group = NULL, country = NULL,
                          missmat_type = c("prior", "fixed", "samples")[1], studycause_map = NULL,
                          missmat = NULL,
                          donotcalib = NULL,
                          donotcalib_type = c("learn", "fixed")[1], nocalib.threshold = 0.1,
                          path_correction = TRUE, ensemble = NULL,
                          pshrink_strength = NULL,
                          nMCMC = 5000, nBurn = 5000, nThin = 1,
                          nChain = 1, nCore = 1,
                          adapt_delta_stan = .9, refresh_stan = NULL,
                          seed = 1, verbose = TRUE, saveoutput = FALSE,
                          output_filename = NULL, output_dir = NULL){


  input_vacalib = as.list(environment())

  # default
  # how frequently print stan sampling status
  if(is.null(refresh_stan)){refresh_stan = .25}

  if(verbose){

    refresh_stan = max(round((nBurn + nMCMC*nThin)*refresh_stan), 1)

  }else{

    refresh_stan = 0

  }

  # default shrinkage towards uncalibrated csmf estimate
  if(path_correction){

    pshrink_strength = 0

  }else{

    if(is.null(pshrink_strength)) pshrink_strength = 4

  }


  # va_data preparation for input into modular_vacalib() ----
  # rewriting as a list of broad-cause-specific death count vector
  if(is.null(va_data)){

    message("")
    message("Provide VA-only (or unlabeled) data in 'va_data' for calibrating.")
    message("")
    message("Input should be as 'va_data'=list('eava'=eava_output, 'insilicova'=insilicova_output, 'interva'=interva_output).")
    message("")
    message("For each algorithm,")
    message("     1. either provide outputs of CCVA algorithms (that is, output from codEAVA() function in EAVA package, or codVA() function in openVA package), or")
    message("     2. provide a named vector of death counts with causes as names.")
    message("")
    message("See the documentation for more details.")
    stop("")

  }else if(!is.list(va_data)){

    message("")
    message("Provide VA-only (or unlabeled) data in 'va_data' for calibrating.")
    message("")
    message("Input should be as 'va_data'=list('eava'=eava_output, 'insilicova'=insilicova_output, 'interva'=interva_output).")
    message("")
    message("For each algorithm,")
    message("     1. either provide outputs of CCVA algorithms (that is, output from codEAVA() function in EAVA package, or codVA() function in openVA package), or")
    message("     2. provide a named vector of death counts with causes as names.")
    message("")
    message("See the documentation for more details.")
    stop("")

  }else if(is.null(names(va_data))){

    message("'va_data' must be a named list, where the names indicate the CCVA algorithms.")
    message("For example, 'va_data'=list('eava'=eava_output, 'insilicova'=insilicova_output, 'interva'=interva_output).")
    message("For using the stored misclassification matrix estimates 'CCVA_missmat', use 'eava', 'insilicova', and 'interva' for EAVA, InSilicoVA, and InterVA.")
    stop("")

  }else{

    if(verbose){

      message("* Preparing VA data for calibration")

    }

    # message(paste0(names(va_data)))

    va_data_tomodel = va_data
    # print(length(va_data))
    # pb = txtProgressBar(min = 1, max = length(va_data), style = 3)
    for(k in 1:length(va_data)){

      if(length(dim(va_data[[k]]))==2){

        if(ncol(va_data[[k]])==2){

          # outputs from EAVA and OpenVA

          ## hard check: whether has columns names like cause_map() outputs
          if((!is.null(colnames(va_data[[k]])))&&(isTRUE(all.equal(sum(!(colnames(va_data[[k]]) %in% c("ID", "cause", "cause1"))), 0)))){

            # broad-cause specific death count vector
            va_data_tomodel[[k]] = colSums(cause_map(df = va_data[[k]], age_group = age_group))

          }else{

            message(paste0("'va_data' for ", names(va_data)[k], " algorithm matches the format of outputs from EAVA or OpenVA, but it must have 'ID' and 'cause'/'cause1' as column names."))
            stop("")

          }

        }else{

          # outputs from cause_map()

          ## hard checks
          ### whether has column names
          if(is.null(colnames(va_data[[k]]))){

            message(paste0("'va_data' for ", names(va_data)[k], " algorithm matches the output format of 'cause_map()'. It must have broad causes as column names."))
            stop("")

          }

          # ### whether has rownames
          # if(is.null(rownames(va_data[[k]]))){
          #
          #   message("")
          #   message(paste0("Note: Rows of 'va_data' for ", names(va_data)[k], " algorithm doesn't have names."))
          #   message("Assuming they correspond to separate individuals.")
          #   message("")
          #
          # }

          ### whether single-cause
          singlecause.check = apply(X = va_data[[k]], 1,
                                    FUN = function(v){

                                      isFALSE(all.equal(sum(v), 1))|(sum(v!=0)>1)

                                    })

          # if(sum(singlecause.check)>=2){
          if(sum(singlecause.check)>0){

            message(va_data[[k]][singlecause.check,])
            message(paste0("'va_data' for ", names(va_data)[k], " algorithm for these rows does not look like single-cause predictions."))
            stop("")

          }

          # broad-cause specific death count vector
          va_data_tomodel[[k]] = colSums(va_data[[k]])

        }

      }else if(length(dim(va_data[[k]]))==0){

        # broad-cause-specific death counts

        ## hard checks: whether is a count vector
        count.check = sum(va_data[[k]])==sum(floor(va_data[[k]]))
        if(!count.check) {

          message(paste0("'va_data' for ", names(va_data)[k], " algorithm must be a named vector of death counts with causes as names."))
          stop("")

        }

        # broad-cause specific death count vector
        va_data_tomodel[[k]] = va_data[[k]]

      }else{

        message(paste0("Unknown format of 'va_data' for ", names(va_data)[k], " algorithm."))
        message("")
        message("Valid input formats for 'va_data' by algorithm:")
        message("     1. For EAVA: Output of codEAVA() in the EAVA package.")
        message("        For InSilicoVA and InterVA: Output of crossVA() in the openVA package.")
        message("        Matrix structured as individuals X 2")
        message("")
        message("     2. Output of cause_map() function in this package.")
        message("        They are broad cause mapping of outputs from codEAVA(), or crossVA().")
        message("        Matrix structured as individual X broad causes")
        message("")
        message("     3. Vector of broad-cause-specific death counts (column sums of output from cause_map()).")
        stop("")

      }
      # print(va_data_tomodel)

      if(k>1){

        if(isTRUE(all.equal(sort(names(va_data_tomodel[[1]])), sort(names(va_data_tomodel[[k]]))))){

          va_data_tomodel[[k]] = va_data_tomodel[[k]][names(va_data_tomodel[[1]])]

        }else{

          message("Causes specified in 'va_data':")
          message(paste0("For algorithm ", names(va_data_tomodel)[1], ": ", paste0(names(va_data_tomodel[[1]]), collapse = ', ')))
          message(paste0("For algorithm ", names(va_data_tomodel)[k], ": ", paste0(names(va_data_tomodel[[k]]), collapse = ', ')))
          message("They must be identical.")
          stop("")

        }

      }

      # if(verbose){setTxtProgressBar(pb, k)}

    } # end of k loop to prepare va data

    if(verbose){

      # message("Preparing VA data for calibration ... Done.")
      message("")

    }

  }
  # print(va_data_tomodel)


  # misclassification matrix for calibration ----
  # missmat.asDirich_study = missmat.fixed_tomodel = NULL
  if(!(missmat_type %in% c("fixed", "prior", "samples"))){

    message("Below are valid options for 'missmat_type':")
    message("     1. 'fixed': Fixed misclassification matrix.")
    message("                 Has nonnegative entries, and each row sums to 1.")
    message("")
    message("     2. 'prior': Prior on the misclassification matrix.")
    message("                 A matrix specifying the Dirichlet prior row by row.")
    message("                 Each row specifies the scale parameters for that row of the misclassification matrix.")
    message("")
    message("     3. 'samples': Misclassification matrix samples.")
    message("                   A 3-dimensional array.")
    message("                   Dimensions respectively represent:")
    message("                       (1) samples,")
    message("                       (2) row of misclassification matrix, and")
    message("                       (3) column of misclassification matrix.")
    stop("")

  }else if(missmat_type=="fixed"){

    ## fixed misclassification matrix ----

    if(is.null(missmat)){

      ### default: CHAMPS ----

      utils::data("CCVA_missmat", package = "vacalibration", envir = environment())
      # print(CCVA_missmat$neonate$eava$postmean$Mozambique)

      if(is.null(age_group)){

        message("'missmat' is not provided. Provide 'age_group' to use CHAMPS-based misclassification matrix estimates stored in 'CCVA_missmat'.")
        message("'age_group' must be 'neonate' (for 0-27 days) or 'child' (for 1-59 months).")
        message("Otherwise specify your own misclassification matrix in 'missmat'. See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
        stop("")

      }

      if(is.null(country)){

        message("'missmat' is not provided. Provide 'country' to use CHAMPS-based misclassification matrix estimates stored in 'CCVA_missmat'.")
        message("CHAMPS-based estimates are available for 'Bangladesh', 'Ethiopia', 'Kenya', 'Mali', 'Mozambique', 'Sierra Leone', or 'South Africa'. For any other country, CHAMPS estimate for 'other' is used and this is done within the package.")
        message("Otherwise specify your own misclassification matrix in 'missmat'. See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
        stop("")

      }

      if(!all(names(va_data_tomodel) %in% c("eava", "insilicova", "interva"))){

        message("'missmat' is not provided. Algorithm names must be provided as component names in 'va_data' to use CHAMPS-based misclassification matrix estimates stored in 'CCVA_missmat'.")
        message("CHAMPS-based estimates are available for 'eava', 'insilicova', and 'interva'.")
        message("Otherwise specify your own misclassification matrix in 'missmat'. See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
        stop("")

      }

      country_tomodel = ifelse(country %in% names(CCVA_missmat[[1]][[1]]$postmean), country, "other")
      # print(country_tomodel)
      if(verbose){

        message(paste0("* Using the misclassification matrix of ", country_tomodel, " for calibration"))
        message("")

      }

      Mmat_input = lapply(1:length(va_data_tomodel),
                          FUN = function(k){

                            CCVA_missmat[[age_group]][[names(va_data_tomodel)[k]]]$postmean[[country_tomodel]]

                          })
      names(Mmat_input) = names(va_data_tomodel)
      # print(Mmat.asDirich)
      # print(CCVA_missmat[[age_group]][[names(va_data)[k]]]$asDirich[[ifelse(country %in% CCVA_missmat[[age_group]][[names(va_data)[k]]]$asDirich, country, "other")]])

    }else{

      ### user provided ----
      if(verbose){

        message(paste0("* Using the user-specified misclassification matrix for calibration"))
        message("")

      }

      # hard checks
      if(!is.list(missmat)){

        message("'missmat' must be a named list like 'va_data'. See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
        stop("")

      }else if(is.null(names(missmat))){

        message("'missmat' must be a named list like 'va_data'. The component names indicate CCVA algorithms, and they must be identical to the component names provided in 'va_data'.")
        message("See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
        stop("")

      }else if(isFALSE(all.equal(sort(names(missmat)), sort(names(va_data_tomodel))))){

        message("'missmat' must be a named list like 'va_data' with component names indicating CCVA algorithms. The algorithm names provided in 'missmat' and 'va_data' must be identical.")
        stop("")

      }else{

        Mmat = missmat[names(va_data_tomodel)]  # ordering missmat and va_data_tomodel according to the algorithms in the same way

        for(k in 1:length(Mmat)){

          nonsimplex.entries = apply(Mmat[[k]], 1,
                                     FUN = function(v){

                                       any(v<0)||isFALSE(all.equal(sum(v), 1))

                                     })

          if(sum(nonsimplex.entries)>0){

            message(paste0("Row index in misclassification matrix: ", paste0(which(nonsimplex.entries), collapse = ", "), "."))
            message(paste0("The above rows in the fixed misclassification matrix provided for ", names(va_data_tomodel)[k], " algorithm are not in the simplex."))
            message("Each row of the matrix must be non-negative and sum to 1.")
            stop("")

          }

          if(!is.matrix(Mmat[[k]])) stop(paste0("Fixed misclassification matrix for ", names(Mmat)[k]," algorithm must be a matrix."))
          if(is.null(rownames(Mmat[[k]]))) stop(paste0("Fixed misclassification matrix for ", names(Mmat)[k], " algorithm must have causes as row names."))
          if(is.null(colnames(Mmat[[k]]))) stop(paste0("Fixed misclassification matrix for ", names(Mmat)[k], " algorithm must have causes as column names."))
          if(isFALSE(all.equal(rownames(Mmat[[k]]), colnames(Mmat[[k]]), names(va_data_tomodel[[k]])))) stop("Causes specified in rows and columns of 'missmat' do not match.")

        }

        Mmat_input = Mmat

      }

    }

    if(all(names(va_data_tomodel[[1]]) %in% colnames(Mmat_input[[1]]))){

      ## study causes subset of CHAMPS ----
      Mmat_temp = lapply(1:length(va_data_tomodel),
                         FUN = function(k){

                           tempmat = Mmat_input[[k]][names(va_data_tomodel[[k]]),names(va_data_tomodel[[k]])]
                           tempmat = tempmat/rowSums(tempmat)
                           rownames(tempmat) = colnames(tempmat) = names(va_data_tomodel[[k]])
                           return(tempmat)

                         })
      names(Mmat_temp) = names(va_data_tomodel)
      Mmat_input = Mmat_temp

    }else{

      ## study causes outside CHAMPS ----

      if(is.null(studycause_map)){

        message(paste0("Study causes: ", names(va_data_tomodel[[1]])))
        message(paste0("CHAMPS causes: ", colnames(CCVA_missmat[[1]][[1]]$postmean[[1]])))
        message("Study causes are not a subset of CHAMPS causes.")
        message("")
        message("Either map study causes to CHAMPS causes. For this specify 'studycause_map' = c('study_cause1' = 'pneumonia', 'study_cause2' = 'ipre', 'study_cause3' = 'other', ...).")
        message("")
        message("Or, specify your own misclassification matrix in 'missmat'. See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
        stop("")

      }

      if(!is.character(studycause_map)){

        message("'studycause_map' must be a named character vector and should provided as 'studycause_map' = c('study_cause1' = 'pneumonia', 'study_cause2' = 'ipre', 'study_cause3' = 'other', ...)")
        stop("")

      }

      if(isFALSE(all(names(va_data_tomodel[[1]]) %in% names(studycause_map)))){

        message(paste0("Match not provided for ", paste0(names(va_data_tomodel[[1]])[which(!(names(va_data_tomodel[[1]]) %in% names(studycause_map)))], collapse = ', ')))
        stop("")

      }

      studycause_map = studycause_map[names(va_data_tomodel[[1]])]

    }

  }else if(missmat_type=="prior"){

    ## prior on misclassification matrix ----

    if(is.null(missmat)){

      ### default: CHAMPS ----

      utils::data("CCVA_missmat", package = "vacalibration", envir = environment())
      # print(CCVA_missmat$neonate$eava$postmean$Mozambique)

      if(is.null(age_group)){

        message("'missmat' is not provided. Provide 'age_group' to use CHAMPS-based misclassification matrix estimates stored in 'CCVA_missmat'.")
        message("'age_group' must be 'neonate' (for 0-27 days) or 'child' (for 1-59 months).")
        message("Otherwise specify your own misclassification matrix in 'missmat'. See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
        stop("")

      }

      if(is.null(country)){

        message("'missmat' is not provided. Provide 'country' to use CHAMPS-based misclassification matrix estimates stored in 'CCVA_missmat'.")
        message("CHAMPS-based estimates are available for 'Bangladesh', 'Ethiopia', 'Kenya', 'Mali', 'Mozambique', 'Sierra Leone', or 'South Africa'. For any other country, CHAMPS estimate for 'other' is used and this is done within the package.")
        message("Otherwise specify your own misclassification matrix in 'missmat'. See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
        stop("")

      }

      if(!all(names(va_data_tomodel) %in% c("eava", "insilicova", "interva"))){

        message("'missmat' is not provided. Algorithm names must be provided as component names in 'va_data' to use CHAMPS-based misclassification matrix estimates stored in 'CCVA_missmat'.")
        message("CHAMPS-based estimates are available for 'eava', 'insilicova', and 'interva'.")
        message("Otherwise specify your own misclassification matrix in 'missmat'. See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
        stop("")

      }

      country_tomodel = ifelse(country %in% names(CCVA_missmat[[1]][[1]]$asDirich), country, "other")
      # print(country_tomodel)
      if(verbose){

        message(paste0("* Using the misclassification matrix of ", country_tomodel, " for calibration"))
        message("")

      }

      # print(1:length(va_data_tomodel))
      # print(names(va_data_tomodel))
      # print(CCVA_missmat[[age_group]][[names(va_data_tomodel)[1]]]$asDirich[[country_tomodel]])
      Mmat_input = lapply(1:length(va_data_tomodel),
                          FUN = function(k){

                            # print(k)
                            # print(CCVA_missmat[[age_group]][[names(va_data_tomodel)[k]]]$asDirich)
                            # print(ifelse(country %in% names(CCVA_missmat[[age_group]][[names(va_data_tomodel)[k]]]$asDirich), country, "other"))
                            # return(CCVA_missmat[[age_group]][[names(va_data_tomodel)[k]]]$asDirich[[country_tomodel]])

                            CCVA_missmat[[age_group]][[names(va_data_tomodel)[k]]]$asDirich[[country_tomodel]]

                          })
      # print(Mmat.asDirich)
      names(Mmat_input) = names(va_data_tomodel)
      # print(Mmat.asDirich)
      # print(CCVA_missmat[[age_group]][[names(va_data)[k]]]$asDirich[[ifelse(country %in% CCVA_missmat[[age_group]][[names(va_data)[k]]]$asDirich, country, "other")]])

    }else{

      ### user provided ----
      if(verbose){

        message(paste0("* Using the user-specified misclassification matrix for calibration"))
        message("")

      }

      # hard checks
      if(!is.list(missmat)){

        message("'missmat' must be a named list like 'va_data'. See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
        stop("")

      }else if(is.null(names(missmat))){

        message("'missmat' must be a named list like 'va_data'. The component names indicate CCVA algorithms, and they must be identical to the component names provided in 'va_data'.")
        message("See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
        stop("")

      }else if(isFALSE(all.equal(sort(names(missmat)), sort(names(va_data_tomodel))))){

        message("'missmat' must be a named list like 'va_data' with component names indicating CCVA algorithms. The algorithm names provided in 'missmat' and 'va_data' must be identical.")
        stop("")

      }else{

        Mmat = missmat[names(va_data_tomodel)]  # ordering missmat and va_data_tomodel according to the algorithms in the same way

        for(k in 1:length(Mmat)){

          nonpositive.entries = apply(Mmat[[k]], 1, FUN = function(v){sum(v<=0)>0})
          if(sum(nonpositive.entries)>0){

            message(paste0("Row index in misclassification matrix: ", paste0(which(nonpositive.entries), collapse = ", "), "."))
            message(paste0("The above rows in the prior misclassification matrix provided for ", names(va_data_tomodel)[k], " algorithm are not Dirichlet scale parameters."))
            message("Each row of the matrix must be strictly positive.")
            stop("")

          }

          if(!is.matrix(Mmat[[k]])) stop(paste0("Prior for ", names(Mmat)[k]," algorithm must be a matrix."))
          if(is.null(rownames(Mmat[[k]]))) stop(paste0("Prior for ", names(Mmat)[k], " algorithm must have causes as row names."))
          if(is.null(colnames(Mmat[[k]]))) stop(paste0("Prior for ", names(Mmat)[k], " algorithm must have causes as column names."))
          if(isFALSE(all.equal(rownames(Mmat[[k]]), colnames(Mmat[[k]]), names(va_data_tomodel[[k]])))) stop("Causes specified in rows and columns of 'missmat' do not match.")

        }

        Mmat_input = Mmat

      }

    }

    if(all(names(va_data_tomodel[[1]]) %in% colnames(Mmat_input[[1]]))){

      ## study causes subset of CHAMPS ----
      Mmat_temp = lapply(1:length(va_data_tomodel),
                         FUN = function(k){

                           tempmat = Mmat_input[[k]][names(va_data_tomodel[[k]]),names(va_data_tomodel[[k]])]
                           rownames(tempmat) = colnames(tempmat) = names(va_data_tomodel[[k]])
                           return(tempmat)

                         })
      names(Mmat_temp) = names(va_data_tomodel)
      Mmat_input = Mmat_temp

    }else{

      ## study causes outside CHAMPS ----

      if(is.null(studycause_map)){

        message(paste0("Study causes: ", names(va_data_tomodel[[1]])))
        message(paste0("CHAMPS causes: ", colnames(CCVA_missmat[[1]][[1]]$postmean[[1]])))
        message("Study causes are not a subset of CHAMPS causes.")
        message("")
        message("Either map study causes to CHAMPS causes. For this specify 'studycause_map' = c('study_cause1' = 'pneumonia', 'study_cause2' = 'ipre', 'study_cause3' = 'other', ...).")
        message("")
        message("Or, specify your own misclassification matrix in 'missmat'. See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
        stop("")

      }

      if(!is.character(studycause_map)){

        message("'studycause_map' must be a named character vector and should provided as 'studycause_map' = c('study_cause1' = 'pneumonia', 'study_cause2' = 'ipre', 'study_cause3' = 'other', ...)")
        stop("")

      }

      if(isFALSE(all(names(va_data_tomodel[[1]]) %in% names(studycause_map)))){

        message(paste0("Match not provided for ", paste0(names(va_data_tomodel[[1]])[which(!(names(va_data_tomodel[[1]]) %in% names(studycause_map)))], collapse = ', ')))
        stop("")

      }

      studycause_map = studycause_map[names(va_data_tomodel[[1]])]

    }

  }else if(missmat_type=="samples"){

    ## misclassification matrix samples ----

    if(is.null(missmat)){

      stop("Provide samples of misclassification matrix.")

    }else{

      # hard checks
      if(!is.list(missmat)){

        message("'missmat' must be a named list like 'va_data'. See the documentation for more details on 'missmat' for missmat_type=='samples'.")
        stop("")

      }else if(is.null(names(missmat))){

        message("'missmat' must be a named list like 'va_data'. The component names indicate CCVA algorithms, and they must be identical to the component names provided in 'va_data'.")
        message("See the documentation for more details on 'missmat' for missmat_type=='samples'.")
        stop("")

      }else if(isFALSE(all.equal(sort(names(missmat)), sort(names(va_data_tomodel))))){

        message("'missmat' must be a named list like 'va_data' with component names indicating CCVA algorithms. The algorithm names provided in 'missmat' and 'va_data' must be identical.")
        stop("")

      }else{

        if(verbose){

          message(paste0("* Using the user-specified misclassification matrix for calibration"))
          message("")

        }

        Mmat = missmat[names(va_data_tomodel)]

        Mmat_input = va_data_tomodel
        for(k in 1:length(Mmat)){

          if(!is.array(Mmat[[k]])){

            message(paste0("Samples for ", names(Mmat)[k]," algorithm must be a 3-dimensional array."))
            message("     Dimensions respectively represent:")
            message("         (1) samples,")
            message("         (2) row of misclassification matrix, and")
            message("         (3) column of misclassification matrix.")
            stop("")

          }

          if(is.null(dimnames(Mmat[[k]])[[2]])) stop(paste0("Each misclassification matrix sample for ", names(Mmat)[k], " algorithm must have causes as row names."))

          if(is.null(dimnames(Mmat[[k]])[[3]])) stop(paste0("Each misclassification matrix sample for ", names(Mmat)[k], " algorithm must have causes as column names."))

          if(isFALSE(all.equal(dimnames(Mmat[[k]])[[2]], dimnames(Mmat[[k]])[[3]]))){

            message(paste0("Row names: ", paste0(dimnames(Mmat[[k]])[[2]], collapse = ", ")))
            message(paste0("Column names: ", paste0(dimnames(Mmat[[k]])[[3]], collapse = ", ")))
            message(paste0("Causes specified in rows and columns of misclassification matrix samples for ", names(Mmat)[k], " algorithm must be identical."))
            stop("")

          }

          # Mmat[[k]] = Mmat[[k]][,names(va_data_tomodel[[k]]),names(va_data_tomodel[[k]])]

          nonsimplex.entries = apply(Mmat[[k]], 1:2,
                                     FUN = function(v){

                                       any(v<0)||isFALSE(all.equal(sum(v), 1))

                                     })

          if(sum(nonsimplex.entries)>0){

            arr.ind = which(nonsimplex.entries, arr.ind = T)

            message(paste0("Sample index: ", paste0(arr.ind[,"row"], collapse = ", "), "."))
            message(paste0("Rows in misclassification matrix: ", paste0(arr.ind[,"col"], collapse = ", "), "."))
            message("Above misclassification matrix samples are not in the simplex.")
            message("For each sample, each row of misclassification matrix must be non-negative and sum to 1.")
            stop("")

          }

          # if(isFALSE(all.equal(dimnames(Mmat[[k]])[[2]], names(va_data_tomodel[[k]])))){
          #
          #   message(paste0("Causes in misclassification samples: ", paste0(dimnames(Mmat[[k]])[[2]], collapse = ", ")))
          #   message(paste0("Causes in 'va_data': ", paste0(names(va_data_tomodel[[k]]), collapse = ", ")))
          #   message(paste0("Causes specified in misclassification matrix samples and 'va_data' for ", names(Mmat)[k], " algorithm do not match."))
          #   stop("")
          #
          # }

          # posterior mean
          Mmat.mean = apply(Mmat[[k]], 2:3, mean)

          Mmat_input[[k]] = do.call("rbind",
                                    lapply(1:dim(Mmat[[k]])[2],
                                           FUN = function(i){

                                             mle.out = nlminb(start = 1,
                                                              lower = 0,
                                                              upper = Inf,
                                                              objective = function(lambda_mmat){

                                                                -sum(LaplacesDemon::ddirichlet(x = (Mmat[[k]])[,i,],
                                                                                               alpha = dim(Mmat[[k]])[2]*lambda_mmat*Mmat.mean[i,],
                                                                                               log = T))

                                                              })

                                             dim(Mmat[[k]])[2]*mle.out$par*Mmat.mean[i,]

                                           }))
          rownames(Mmat_input[[k]]) = colnames(Mmat_input[[k]]) = colnames(Mmat.mean)

        }

      }

      if(all(names(va_data_tomodel[[1]]) %in% colnames(Mmat_input[[1]]))){

        ## study causes subset of CHAMPS ----
        Mmat_temp = lapply(1:length(va_data_tomodel),
                           FUN = function(k){

                             tempmat = Mmat_input[[k]][names(va_data_tomodel[[k]]),names(va_data_tomodel[[k]])]
                             rownames(tempmat) = colnames(tempmat) = names(va_data_tomodel[[k]])
                             return(tempmat)

                           })
        names(Mmat_temp) = names(va_data_tomodel)
        Mmat_input = Mmat_temp

      }else{

        ## study causes outside CHAMPS ----

        if(is.null(studycause_map)){

          message(paste0("Study causes: ", names(va_data_tomodel[[1]])))
          message(paste0("CHAMPS causes: ", colnames(CCVA_missmat[[1]][[1]]$postmean[[1]])))
          message("Study causes are not a subset of CHAMPS causes.")
          message("")
          message("Either map study causes to CHAMPS causes. For this specify 'studycause_map' = c('study_cause1' = 'pneumonia', 'study_cause2' = 'ipre', 'study_cause3' = 'other', ...).")
          message("")
          message("Or, specify your own misclassification matrix in 'missmat'. See the documentation for more details on 'missmat' for missmat_type=='fixed'.")
          stop("")

        }

        if(!is.character(studycause_map)){

          message("'studycause_map' must be a named character vector and should provided as 'studycause_map' = c('study_cause1' = 'pneumonia', 'study_cause2' = 'ipre', 'study_cause3' = 'other', ...)")
          stop("")

        }

        if(isFALSE(all(names(va_data_tomodel[[1]]) %in% names(studycause_map)))){

          message(paste0("Match not provided for ", paste0(names(va_data_tomodel[[1]])[which(!(names(va_data_tomodel[[1]]) %in% names(studycause_map)))], collapse = ', ')))
          stop("")

        }

        studycause_map = studycause_map[names(va_data_tomodel[[1]])]

      }

    }

  }


  # causes not calibrated ----
  if(is.null(donotcalib)){

    donotcalib = lapply(1:length(va_data_tomodel),
                        FUN = function(k){

                          'other'

                        })
    names(donotcalib) = names(va_data_tomodel)
    donotcalib_study = donotcalib

  }else{

    if(!is.list(donotcalib)){

      stop("'donotcalib' must be a named list like 'va_data' where the component names indicate CCVA algorithms. It indicates algorithm-specific causes that will not be calibrated.")

    }else if(is.null(names(donotcalib))||isFALSE(all.equal(sort(names(donotcalib)), sort(names(va_data_tomodel))))){

      stop("'donotcalib' is a named list like 'va_data' where component names indicate CCVA algorithms. The names must match with the names provided in 'va_data'.")

    }else{

      donotcalib = donotcalib[names(va_data_tomodel)]

      donotcalib_study = lapply(1:length(va_data_tomodel),
                                FUN = function(k){

                                  if(is.null(donotcalib[[k]])){

                                    return(donotcalib[[k]])

                                  }else{

                                    if((!all(donotcalib[[k]] %in% names(va_data_tomodel[[k]])))&(!all(donotcalib[[k]] %in% unique(studycause_map)))){

                                      message(paste0("ALL causes provided in 'donotcalib' for ", names(va_data_tomodel)[k], " algorithm must either be observed in 'va_data', or they must be a subset of CHAMPS broad causes if using stored misclassification estimates in 'CCVA_missmat'."))
                                      message("")
                                      message(paste0("Causes not observed in 'va_data': ", paste0(donotcalib[[k]][!(donotcalib[[k]] %in% names(va_data_tomodel[[k]]))], collapse = ', ')))
                                      message(paste0("Causes not CHAMPS broad causes: ", paste0(donotcalib[[k]][!(donotcalib[[k]] %in% unique(studycause_map))], collapse = ', ')))
                                      stop("")

                                    }else{

                                      if(all(donotcalib[[k]] %in% names(va_data_tomodel[[k]]))){

                                        return(donotcalib[[k]])

                                      }else if(all(donotcalib[[k]] %in% unique(studycause_map))){

                                        return(names(studycause_map)[studycause_map %in% donotcalib[[k]]])

                                      }

                                    }

                                  }

                                })
      names(donotcalib_study) = names(va_data_tomodel)

    }

  }

  # update input list
  # print(input_vacalib)
  input_vacalib.now = as.list(environment())
  # print(input_vacalib.now)
  input_vacalib = input_vacalib.now[names(input_vacalib)]
  # print(input_vacalib)


  # run calibration ----
  # if(verbose){
  #
  #   # cat("\n")
  #   message("Calibrating ...")
  #
  # }

  if(missmat_type %in% c('prior', 'samples')){

    modular_vacalib_output = modular_vacalib_prior(va_unlabeled = va_data_tomodel,
                                                   Mmat_calib = Mmat_input, studycause_map = studycause_map,
                                                   donotcalib = donotcalib_study,
                                                   donotcalib_type = donotcalib_type, nocalib.threshold = nocalib.threshold,
                                                   path_correction = path_correction, ensemble = ensemble,
                                                   # shrink_towards = shrink_towards,
                                                   pshrink_strength = pshrink_strength,
                                                   nMCMC = nMCMC, nBurn = nBurn, nThin = nThin,
                                                   nChain = nChain, nCore = nCore,
                                                   adapt_delta_stan = adapt_delta_stan, refresh_stan = refresh_stan,
                                                   seed = seed, verbose = verbose,
                                                   input_vacalib = input_vacalib)

  }else if(missmat_type %in% c('fixed')){

    modular_vacalib_output = modular_vacalib_fixed(va_unlabeled = va_data_tomodel,
                                                   Mmat_calib = Mmat_input, studycause_map = studycause_map,
                                                   donotcalib = donotcalib_study,
                                                   donotcalib_type = donotcalib_type, nocalib.threshold = nocalib.threshold,
                                                   path_correction = path_correction, ensemble = ensemble,
                                                   # shrink_towards = shrink_towards,
                                                   pshrink_strength = pshrink_strength,
                                                   nMCMC = nMCMC, nBurn = nBurn, nThin = nThin,
                                                   nChain = nChain, nCore = nCore,
                                                   adapt_delta_stan = adapt_delta_stan, refresh_stan = refresh_stan,
                                                   seed = seed, verbose = verbose,
                                                   input_vacalib = input_vacalib)

  }

  if(verbose){

    message("* VA-Calibration complete")
    message("")


    # print output ----
    cat("\n")
    message("========================")
    message("Summary:")
    message("========================")
    if(sum(modular_vacalib_output$donotcalib_tomodel)>0){

      message("* Causes not calibrated:")

      for(k in 1:nrow(modular_vacalib_output$donotcalib_tomodel)){

        if(sum(modular_vacalib_output$donotcalib_tomodel[k,])>0){message(paste0("     ", rownames(modular_vacalib_output$donotcalib_tomodel)[k], ": ", paste0((colnames(modular_vacalib_output$donotcalib_tomodel))[modular_vacalib_output$donotcalib_tomodel[k,]], collapse = ', ')))}

      }

    }

    cat("\n")
    message("* Cause-specific mortality fractions (in %):")
    print(knitr::kable(modular_vacalib_output$puncalib_percentage,
                       # format = "markdown",
                       caption = "Uncalibrated estimate"))
    print(knitr::kable(modular_vacalib_output$pcalib_postmean_percentage,
                       # format = "markdown",
                       caption = "Calibrated point estimate"))
    # message(print(knitr::kable(round(100*pcalib_postmean), format = "markdown",
    #                            caption = "Calibrated point estimate of cause-specific mortality fractions (in %)")))
    # message("")

    if(!modular_vacalib_output$input$ensemble){

      cat("\n")
      message("* Number of deaths:")
      print(knitr::kable(data.frame(modular_vacalib_output$va_deaths_uncalib,
                                    "Total" = rowSums(modular_vacalib_output$va_deaths_uncalib)),
                         # format = "markdown",
                         caption = "Observed"))
      print(knitr::kable(data.frame(modular_vacalib_output$va_deaths_calib_algo,
                                    "Total" = rowSums(modular_vacalib_output$va_deaths_calib_algo)),
                         # format = "markdown",
                         caption = "Calibrated point estimate"))
      # message(print(knitr::kable(modular_vacalib_output$va_deaths_calib_algo, format = "markdown",
      #                            caption = "Calibrated point estimate of death counts")))
      message("")

    }else{

      cat("\n")
      message("* Number of deaths:")
      print(knitr::kable(data.frame(modular_vacalib_output$va_deaths_uncalib,
                                    "Total" = rowSums(modular_vacalib_output$va_deaths_uncalib)),
                         # format = "markdown",
                         caption = "Observed"))
      print(knitr::kable(data.frame(modular_vacalib_output$va_deaths_calib_algo,
                                    "Total" = rowSums(modular_vacalib_output$va_deaths_calib_algo)),
                         # format = "markdown",
                         caption = "Point estimate from algorithm-specific calibration"))
      print(knitr::kable(data.frame(modular_vacalib_output$va_deaths_calib_ensemble,
                                    "Total" = rowSums(modular_vacalib_output$va_deaths_calib_ensemble)),
                         # format = "markdown",
                         caption = "Point estimate from ensemble calibration"))

      # message(print(knitr::kable(modular_vacalib_output$va_deaths_calib_algo, format = "markdown",
      #                            caption = "Point estimate of death counts (Algorithm-specific calibration)")))
      # message(print(knitr::kable(modular_vacalib_output$va_deaths_calib_ensemble, format = "markdown",
      #                            caption = "Point estimate of death counts (Ensemble calibration)")))
      message("")

    }

  }


  # return output ----
  class(modular_vacalib_output) = "vacalibration"

  if(saveoutput){

    if(is.null(output_filename)){output_filename = "vacalibration_out"}

    if(is.null(output_dir)){

      saveRDS(modular_vacalib_output, file.path(tempdir(), output_filename))

    }else{saveRDS(modular_vacalib_output, file.path(output_dir, output_filename))}

  }else{

    return(modular_vacalib_output)

  }

}


