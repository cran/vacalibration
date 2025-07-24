#' Modular VA-Calibration
#'
#' @param va_unlabeled A named list. Algorithm-specific unlabeled VA-only data.
#'
#' For example, \code{list("algo1" = algo1_output, "algo2" = algo2_output, ...)}
#'
#' Algorithm names (\code{"algo1"}, \code{"algo2"}, ...) can be "eava", "insilicova", or "interva".
#'
#' Data (\code{algo1_output}, \code{algo2_output}, ...) can be broad causes (output from the \code{cause_map()} function in this package), or broad-cause-specific death counts (integer vector).
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
#' \code{Mmat.asDirich_algo1} is a matrix of dimension CHAMPS ("gold standard") cause by VA cause.
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
#'
#' @return A named list. Use \code{vacalibration()} for general purpose.
#'
#' @import patchwork
#' @import rstan
#' @import loo
#' @import ggplot2
#'
#' @importFrom stats quantile
#' @importFrom utils head
modular.vacalib <- function(va_unlabeled = NULL, age_group = NULL,
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

  # saving all input
  input.list = as.list(environment())

  # how frequently print stan sampling status
  if(is.null(refresh.stan)){

    if(verbose){

      refresh.stan = max((nBurn + nMCMC*nThin)/10, 1)

    }else{

      refresh.stan = 0

    }

  }

  # default shrinkage towards uncalibrated csmf estimate
  # p_calib ~ Dirichlet( 1 + nCause*pss*p_uncalib)
  if(stable){

    pss = 0

  }else{

    if(is.null(pss)) pss = 4

  }

  ## preprocessing va unlabeled data ----
  if(!is.list(va_unlabeled)){
    ### single algorithm ----

    K = 1 # number of algorithms

    #### va unlabeled data ----
    if(!any(is.vector(va_unlabeled), is.matrix(va_unlabeled))){

      stop("'va_unlabeled' must be either a vector (death counts for each cause), or a matrix (individuals along rows, causes along columns).")

    }else if(is.matrix(va_unlabeled)){

      if(is.null(colnames(va_unlabeled))) stop("Columns of 'va_unlabeled' must have causes as names.")

      causes = colnames(va_unlabeled)
      nCause = length(causes)

      if(is.null(rownames(va_unlabeled))) message("Rows of 'va_unlabeled' doesn't have names. Assuming they correspond to individuals.")

      singlecause.check = (rowSums(va_unlabeled)==1)&(apply(X = va_unlabeled, 1, FUN = function(v){sum(v!=0)})==1)
      if((!all(singlecause.check))){

        # cat("\n")
        if(verbose) message(va_unlabeled[singlecause.check,])
        # cat("\n")
        stop("'va_unlabeled' doesn't look like single cause predictions for the above cases.")

      }

      if(verbose) message("No algorithm name provided. Naming it as 'algorithm'.")
      va_deaths = array(dim = c(K, nCause),
                        dimnames = list('algorithm', causes))
      va_deaths[1,] = colSums(va_unlabeled)

    }else if(is.vector(va_unlabeled)){

      if(is.null(names(va_unlabeled))) stop("Components of 'va_unlabeled' must have causes as names.")

      causes = names(va_unlabeled)
      nCause = length(causes)

      count.check = sum(va_unlabeled)==sum(floor(va_unlabeled))
      if(!count.check) stop("'va_unlabeled' must be a vector of death counts for each cause.")

      va_deaths = array(dim = c(K, nCause),
                        dimnames = list('algorithm', causes))
      va_deaths[1,] = va_unlabeled

    }

    # to_calib_flat = to_calib
    # nto_calib = c(0, length(to_calib))


    #### misclassification matrices ----
    if(!(calibmodel.type %in% c("Mmatfixed", "Mmatprior"))) {

      stop("'calibmodel.type' must be 'Mmatfixed' or 'Mmatprior'")

    }else if(calibmodel.type=="Mmatfixed") {

      ##### fixed Mmat ----

      # checks
      if(is.null(Mmat.fixed)) stop("Must specify 'Mmat.fixed' for calibmodel.type 'Mmatfixed'.")

      Mmat.fixed = list('algorithm' = Mmat.fixed)

    }else if(calibmodel.type=="Mmatprior") {

      ##### Mmat as prior ----

      # checks
      if(is.null(Mmat.asDirich)) stop("Must specify 'Mmat.asDirich' for calibmodel.type 'Mmatprior'.")

      Mmat.asDirich = list('algorithm' = Mmat.asDirich)

    }

    #### causes not to calibrate ----
    donotcalib = list('algorithm' = donotcalib)

  }else{
    ### multiple algorithms ----

    if(is.null(names(va_unlabeled))) stop("Components of 'va_unlabeled' must have algorithms as names.")

    K = length(va_unlabeled) # number of algorithms

    for(k in 1:K){

      #### va unlabeled data ----
      if(!any(is.vector(va_unlabeled[[k]]), is.matrix(va_unlabeled[[k]]))){

        stop("'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," must be either a vector (death counts for each cause), or a matrix (individuals along rows, causes along columns).")

      }else if(is.matrix(va_unlabeled[[k]])){

        if(is.null(colnames(va_unlabeled[[k]]))) stop(paste0("Columns of 'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," must have causes as names."))

        if(is.null(rownames(va_unlabeled[[k]]))) message(paste0("Rows of 'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," doesn't have names. Assuming they correspond to individuals."))

        singlecause.check = (rowSums(va_unlabeled[[k]])==1)&(apply(X = va_unlabeled[[k]], 1, FUN = function(v){sum(v!=0)})==1)
        if(!all(singlecause.check)){

          # cat("\n")
          if(verbose) message(va_unlabeled[[k]][singlecause.check,])
          # cat("\n")
          stop(paste0("'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," doesn't look like single cause predictions for the above cases."))

        }

        causes_temp = colnames(va_unlabeled[[k]])

        if(k==1){

          causes = causes_temp
          nCause = length(causes)

          va_deaths = array(dim = c(K, nCause),
                            dimnames = list(names(va_unlabeled), causes))

        }else if(k>1){

          if(!identical(causes_temp, causes)) stop("Causes input through 'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," does not match with algorithm ", names(va_unlabeled)[k-1],".")

        }

        va_deaths[k,] = colSums(va_unlabeled[[k]])

      }else if(is.vector(va_unlabeled[[k]])){

        if(is.null(names(va_unlabeled[[k]]))) stop("Components of 'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," must have causes as names.")

        causes_temp = names(va_unlabeled[[k]])

        count.check = sum(va_unlabeled[[k]])==sum(floor(va_unlabeled[[k]]))
        if(!count.check) stop(paste0("'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," doesn't look like a count vector. It must be a vector of death counts for each cause."))

        if(k==1){

          causes = causes_temp
          nCause = length(causes)

          va_deaths = array(dim = c(K, nCause),
                            dimnames = list(names(va_unlabeled), causes))

        }else if(k>1){

          if(!identical(causes_temp, causes)) stop("Causes input through 'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," does not match with algorithm ", names(va_unlabeled)[k-1],".")

        }

        va_deaths[k,] = va_unlabeled[[k]]

      }

    }

    # to_calib_flat = unlist(to_calib)
    # nto_calib = c(0, cumsum(sapply(vectors, length)))

  }
  # print(va_deaths)

  ## observed uncalibrated csmf ----
  (puncalib = va_deaths/rowSums(va_deaths))
  # if(any(puncalib==0)){
  #
  #   puncalib[which.max(puncalib)] = puncalib[which.max(puncalib)] - 0.001
  #   puncalib[puncalib==0] = 0.001/sum(puncalib==0)
  #   puncalib[which.max(puncalib)] = 1 - sum(puncalib[-which.max(puncalib)])
  #
  # }


  ## calibration modeling ----
  if(!(calibmodel.type %in% c("Mmatfixed", "Mmatprior"))) {

    stop("'calibmodel.type' must be 'Mmatfixed' or 'Mmatprior'")

  }else if(calibmodel.type=="Mmatfixed") {

    ### fixed Mmat ----
    seqcalib_mmat_model = readRDS(system.file("stan/seqcalib_mmat.rds", package = "vacalibration"))

    #### hard checks ----
    if(!is.list(Mmat.fixed)) stop("'Mmat.fixed' must be a list of matrices for each algorithm. Rows and columns of each matrix must correspond to CHAMPS and VA causes.")
    if(is.null(names(Mmat.fixed))) stop("Components of 'Mmat.fixed' must have algorithms as names.")
    if(!identical(names(Mmat.fixed), dimnames(va_deaths)[[1]])) stop("Algorithm names specified in 'Mmat.fixed' and 'va_unlabeled' do not match.")

    #### causes not calibrated ----
    if(is.null(donotcalib)){

      donotcalib = lapply(1:K, FUN = function(k){NULL})
      names(donotcalib) = dimnames(va_deaths)[[1]]

    }else{

      if(!is.list(donotcalib)){

        if(is.character(donotcalib)){

          donotcalib = lapply(1:K, FUN = function(k){donotcalib})
          names(donotcalib) = dimnames(va_deaths)[[1]]

        }

        # stop("'donotcalib' must be a list of character vectors for each algorithm. Each vector indicates the causes not to calibrate.")

      }else{

        if(is.null(names(donotcalib))){

          stop("Components of 'donotcalib' must have algorithms as names.")

        }else if(!identical(names(donotcalib), dimnames(va_deaths)[[1]])){

          stop("Either names are not provided, or they do not match with that provided in 'va_unlabeled'")

        }

      }

    }

    if(!identical(names(donotcalib), dimnames(va_deaths)[[1]])) stop("Algorithm names specified in 'donotcalib' and 'va_unlabeled' do not match.")



    #### calibration for each algorithm ----

    # storage
    Mmat_input_asarray = array(dim = c(K, nCause, nCause),
                               dimnames = list(names(Mmat.fixed), causes, causes))
    Mmat_asarray = array(dim = c(K, nCause, nCause),
                         dimnames = list(names(Mmat.fixed), causes, causes))
    donotcalib_asmat = matrix(nrow = K, ncol = nCause)
    causes_notcalibrated = calibout = vector(mode = "list", length = K)
    rownames(donotcalib_asmat) = names(causes_notcalibrated) =
      names(calibout) = names(Mmat.fixed)
    colnames(donotcalib_asmat) = causes
    va_deaths_calib = va_deaths
    pcalib_asarray = array(dim = c(K, nMCMC, nCause),
                           dimnames = list(names(Mmat.fixed), NULL, causes))
    pcalib_postsumm_asarray = array(dim = c(K, 3, nCause),
                                    dimnames = list(names(Mmat.fixed), c('postmean', 'lowcredI', 'upcredI'), causes))

    for(k in 1:K){

      # hard checks
      if(!is.matrix(Mmat.fixed[[k]])) stop(paste0("'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k]," must be a matrix with CHAMPS causes along the rows and VA causes along the columns."))
      if(is.null(colnames(Mmat.fixed[[k]]))) stop(paste0("'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k], " must have VA causes as column names."))
      if(!identical(colnames(Mmat.fixed[[k]]), causes)) stop(paste0("VA causes specified in 'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k], " do not match with that specified in 'va_unlabeled'."))
      if(is.null(rownames(Mmat.fixed[[k]]))) stop(paste0("'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k], " must have CHAMPS causes as row names."))
      if(!identical(rownames(Mmat.fixed[[k]]), causes)) stop(paste0("CHAMPS causes specified in 'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k], " do not match with that specified in 'va_unlabeled'."))
      if(!all.equal(target = unname(rowSums(Mmat.fixed[[k]])), current = rep(1, nrow(Mmat.fixed[[k]])))) stop(paste0("Each row of 'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k], " must sum to 1."))
      # if(any(rowSums(Mmat.fixed)!=1)) stop(paste0("Each row of 'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k], " must sum to 1."))

      # causes not to calibrate
      if((!is.null(donotcalib[[k]])) & (sum(causes %in% donotcalib[[k]])==0)){

        stop("'donotcalib' for algorithm ", names(donotcalib)[k],"  doesn't match with causes named in 'va_unlabeled' and 'Mmat.fixed'")

      }
      # if(('other' %in% causes) & (!('other' %in% donotcalib[[k]]))) donotcalib[[k]] = c(donotcalib[[k]], 'other')

      (donotcalib_asmat[k,] = causes %in% donotcalib[[k]])

      # learn more causes not to calibrate
      if(donot.calib_type=="learn"){

        (donot.calib_Mmat_k = apply(Mmat.fixed[[k]], 2,
                                    FUN = function(v){

                                      # mean(abs(v - mean(v)))
                                      diff(range(v))<=nocalib.threshold

                                    }))
        (donotcalib_asmat[k,] = donotcalib_asmat[k,]|donot.calib_Mmat_k)

      }
      causes_notcalibrated[[k]] = causes[donotcalib_asmat[k,]]

      ##### preparing Mmat for calibration ----
      if(sum(!donotcalib_asmat[k,])<=1){

        calibrated_k = FALSE
        lambda_k = 1
        Mmat_input_asarray[k,,] = Mmat_asarray[k,,] = diag(nCause)

      }else{

        calibrated_k = TRUE

        # shrink Mmat for stability
        if(!stable){

          lambda_k = 0

          # sub Mmat for calibrates causes
          Mmat_fixed_temp_k = diag(nCause)

          Mmat_fixed_sub_k = Mmat.fixed[[k]][!donotcalib_asmat[k,],!donotcalib_asmat[k,]]
          Mmat_fixed_sub_k = Mmat_fixed_sub_k/rowSums(Mmat_fixed_sub_k)

          # print(k)
          # print(puncalib)
          # print(donotcalib_asmat)
          puncalib_sub_k = puncalib[k,!donotcalib_asmat[k,]]/sum(puncalib[k,!donotcalib_asmat[k,]])

          Mmat_fixed_temp_k[!donotcalib_asmat[k,],!donotcalib_asmat[k,]] = Mmat_fixed_sub_k

          Mmat_input_asarray[k,,] = Mmat_asarray[k,,] = Mmat_fixed_temp_k

        }else{

          Mmat_fixed_temp_k = diag(nCause)

          # sub Mmat for calibrates causes
          Mmat_fixed_sub_k = Mmat.fixed[[k]][!donotcalib_asmat[k,],!donotcalib_asmat[k,]]
          Mmat_fixed_sub_k = Mmat_fixed_sub_k/rowSums(Mmat_fixed_sub_k)

          # print(k)
          # print(puncalib)
          # print(donotcalib_asmat)
          puncalib_sub_k = puncalib[k,!donotcalib_asmat[k,]]/sum(puncalib[k,!donotcalib_asmat[k,]])

          Mmat_fixed_temp_k[!donotcalib_asmat[k,],!donotcalib_asmat[k,]] = Mmat_fixed_sub_k

          Mmat_input_asarray[k,,] = Mmat_fixed_temp_k

          count.lambda_k = 0
          got_it_k = FALSE
          while(!got_it_k){

            # shrink coeff
            count.lambda_k = count.lambda_k + 1
            if(count.lambda_k==1){

              lambda_k = 1

            }else{

              lambda_k = lambda_k - .01

            }

            # shrinking towards identity
            (Mmat_lambda_k = lambda_k*diag(nrow(Mmat_fixed_sub_k)) + (1-lambda_k)*Mmat_fixed_sub_k)

            # solving calibration eq
            pcalib_lambda_k = as.numeric(solve(Mmat_lambda_k %*% t(Mmat_lambda_k)) %*%
                                           Mmat_lambda_k %*% as.matrix(puncalib_sub_k))

            # pcalib_lambda = puncalib
            # pcalib_lambda[!donotcalib_asmat[k,]] = (1-sum(puncalib[donotcalib_asmat[k,]]))*pcalib_lambda

            # checking if within simplex
            within_simplex_k = (sum(pcalib_lambda_k<0)==0)&all.equal(sum(pcalib_lambda_k), 1)
            got_it_k = (!within_simplex_k)|isTRUE(all.equal(lambda_k, 0))
            # got_it = !within_simplex

          }

          if((lambda_k<1) & (!within_simplex_k)){

            lambda_k = lambda_k + .01

          }

          # shrinked Mmat
          (Mmat_lambda_opt_k = lambda_k*diag(nrow(Mmat_fixed_sub_k)) + (1-lambda_k)*Mmat_fixed_sub_k)
          rownames(Mmat_lambda_opt_k) = colnames(Mmat_lambda_opt_k) = causes[!donotcalib_asmat[k,]]

          Mmat_fixed_temp_k[!donotcalib_asmat[k,],!donotcalib_asmat[k,]] = Mmat_lambda_opt_k

          Mmat_asarray[k,,] = Mmat_fixed_temp_k

        }

      }


      #### calibration output ----
      if(!calibrated_k){

        ###### not calibrated ----
        MCMCout_k = list('p_calib' = matrix(data = puncalib[k,],
                                            nrow = nMCMC, ncol = nCause,
                                            byrow = TRUE),
                         'loglik' = NULL,
                         "calibrated" = calibrated_k, "lambda" = lambda_k,
                         'loo.out' = NULL, 'waic.out' = NULL,
                         'ic.df' = NULL, 'mcmc.diagnostic' = NULL,
                         'p_calib_postsumm' = rbind(puncalib[k,],
                                                    puncalib[k,],
                                                    puncalib[k,]))

        # posterior summary of calibrated estimate
        colnames(MCMCout_k$p_calib) = causes
        rownames(MCMCout_k$p_calib_postsumm) = c('postmean', 'lowcredI', 'upcredI')

      }else{

        ###### calibrated ----
        stanfit_k = rstan::sampling(seqcalib_mmat_model,
                                    pars = c('p_calib', 'loglik'),
                                    include = TRUE,
                                    data = list('nCause' = sum(!donotcalib_asmat[k,]),
                                                'nAlgo' = 1,
                                                'aj' = va_deaths[k,!donotcalib_asmat[k,], drop = FALSE],
                                                'p_uncalib' = puncalib_sub_k,
                                                'Mmat_bysimplex' = Mmat_asarray[k,!donotcalib_asmat[k,],!donotcalib_asmat[k,], drop = FALSE],
                                                'pss' = pss
                                    ),
                                    chains = 1,
                                    iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                                    control = list('adapt_delta' = adapt_delta_stan),
                                    seed = seed,
                                    refresh = refresh.stan#,
                                    # init = list(list('Mmat' = Mmat.init,
                                    #                  'p_calib' = puncalib))
        )

        MCMCout_k = rstan::extract(stanfit_k) # MCMC output
        MCMCout_k$calibrated = calibrated_k
        MCMCout_k$lambda = lambda_k

        # mcmc diagnostic
        # max Rhat
        max_Rhat_k = max(apply(X = MCMCout_k$p_calib, 2,
                               FUN = function(v){

                                 rstan::Rhat(v)

                               }))

        # min bulk ESS
        min_ess_bulk_k = min(apply(X = MCMCout_k$p_calib, 2,
                                   FUN = function(v){

                                     rstan::ess_bulk(v)

                                   }))/nMCMC

        # loo ic
        MCMCout_k$loo.out = loo::loo(MCMCout_k$loglik,
                                     r_eff = loo::relative_eff(exp(MCMCout_k$loglik),
                                                               chain_id = rep(1, nrow(MCMCout_k$loglik))),
                                     cores = 1)

        # waic
        MCMCout_k$waic.out = loo::waic(MCMCout_k$loglik)

        ic.df_k = rbind(MCMCout_k$waic.out$estimates,
                        MCMCout_k$loo.out$estimates)

        ic.df.melt_k = as.numeric(t(ic.df_k))
        names(ic.df.melt_k) = paste0(rep(rownames(ic.df_k), each = ncol(ic.df_k)),'_',rep(colnames(ic.df_k), nrow(ic.df_k)))

        MCMCout_k$ic.df = ic.df.melt_k

        # mcmc diagnostic summary
        MCMCout_k$mcmc.diagnostic = c('max_Rhat' = max_Rhat_k, 'min_ess_bulk' = min_ess_bulk_k,
                                      'num_divergent' = rstan::get_num_divergent(stanfit_k),
                                      'num_max_treedepth' = rstan::get_num_max_treedepth(stanfit_k))

        if(verbose){

          # cat("\n")
          message(MCMCout_k$mcmc.diagnostic)
          # cat("\n")

        }

        # calibrated posterior for all causes
        # print(head(MCMCout$p_calib))
        if(sum(donotcalib_asmat[k,])>0){

          # csmf
          p_calib_out_k = matrix(data = puncalib[k,],
                                 nrow = nMCMC, ncol = nCause,
                                 byrow = TRUE)
          # print(dim(p_calib_out[,donotcalib_touse, drop = F]))
          # p_calib_out[,donotcalib_touse] = matrix(data = puncalib[donotcalib_touse],
          #                                         nrow = nMCMC, ncol = sum(donotcalib_touse),
          #                                         byrow = T)
          p_calib_out_k[,!donotcalib_asmat[k,]] = (1 - sum(puncalib[k,donotcalib_asmat[k,]]))*MCMCout_k$p_calib
          MCMCout_k$p_calib = p_calib_out_k

        }
        colnames(MCMCout_k$p_calib) = causes
        # MCMCout$p_calib.postmean = colMeans(MCMCout$p_calib)
        # print(head(MCMCout$p_calib))

        # posterior summary of calibrated estimate
        MCMCout_k$p_calib_postsumm = apply(MCMCout_k$p_calib, 2,
                                           FUN = function(v){

                                             c(mean(v),
                                               quantile(x = v, probs = c(.025, .975)))

                                           })
        rownames(MCMCout_k$p_calib_postsumm) = c('postmean', 'lowcredI', 'upcredI')

      }

      # store MCMC output
      calibout[[k]] = MCMCout_k

      #### calibrated number of deaths ----
      va_deaths_calib_k = va_deaths[k,]
      if(calibrated_k){

        va_deaths_calib_k[!donotcalib_asmat[k,]] = round(colMeans(sum(va_deaths[k,!donotcalib_asmat[k,]])*
                                                                    (MCMCout_k$p_calib[,!donotcalib_asmat[k,]]/(1 - sum(puncalib[k,donotcalib_asmat[k,]])))))

        # exactly match total number of deaths
        if(!isTRUE(all.equal(sum(va_deaths_calib_k), sum(va_deaths[k,])))){

          (adjust.amnt_k = sum(va_deaths_calib_k) - sum(va_deaths[k,]))

          adjust.vec_k = rep(0, nCause)
          (adjust.vec_k[!donotcalib_asmat[k,]] = rep(floor(abs(adjust.amnt_k)/sum(!donotcalib_asmat[k,])),
                                                     sum(!donotcalib_asmat[k,])))

          id.extra.adjust_k = head(intersect(order(va_deaths_calib_k, decreasing = TRUE),
                                             which(!donotcalib_asmat[k,])),
                                   n = abs(adjust.amnt_k) %% sum(!donotcalib_asmat[k,]))
          adjust.vec_k[id.extra.adjust_k] = adjust.vec_k[id.extra.adjust_k] + 1
          adjust.vec_k

          if(adjust.amnt_k>0){

            va_deaths_calib_k = va_deaths_calib_k - adjust.vec_k

          }else{

            va_deaths_calib_k = va_deaths_calib_k + adjust.vec_k

          }

        }

      }

      va_deaths_calib[k,] = va_deaths_calib_k
      pcalib_asarray[k,,] = MCMCout_k$p_calib
      pcalib_postsumm_asarray[k,,] = MCMCout_k$p_calib_postsumm

    }# ending calibration for each algorithm


    ### ensemble calibration  ----
    if(K==1){

      if(is.null(ensemble)){

        ensemble = FALSE

      }else{

        if(ensemble){

          ensemble = FALSE
          if(verbose) message("Nothing to ensemble. Only one algorithm provided. Setting ensemble=F.")

        }

      }

    }else if(K>1){

      # defaulting to ensemble if multiple algorithms
      if(is.null(ensemble)) ensemble = TRUE

    }

    # checking if needs to run ensemble
    if(!ensemble){

      # output list
      input.list.now = as.list(environment())
      input.list = c(input.list.now[names(input.list)],
                     "K" = K, "nCause" = nCause, "causes" = causes)

      output = list("calib_MCMCout" = calibout,
                    "p_uncalib" = puncalib,
                    "p_calib" = pcalib_asarray,
                    "pcalib_postsumm" = pcalib_postsumm_asarray,
                    "va_deaths_uncalib" = va_deaths,
                    "va_deaths_calib_algo" = va_deaths_calib,
                    'Mmat.fixed_input' = Mmat_input_asarray,
                    'Mmat.fixed' = Mmat_asarray,
                    'donotcalib' = donotcalib_asmat,
                    'causes_notcalibrated' = causes_notcalibrated,
                    'input' = input.list)

    }else{

      # ensemble uncalibrated csmf
      (puncalib_ens = colSums(va_deaths)/sum(va_deaths))

      donotcalib_ens = donotcalib_asmat[1,]
      for(k in 2:K){

        donotcalib_ens = donotcalib_ens&donotcalib_asmat[k,]

      }

      calibrated_ens = sum(!donotcalib_ens)>1
      if(!calibrated_ens){

        MCMCout_ens = list('p_calib' = matrix(data = puncalib_ens,
                                              nrow = nMCMC, ncol = nCause,
                                              byrow = TRUE),
                           'loglik' = NULL,
                           "calibrated" = calibrated_ens,
                           'loo.out' = NULL, 'waic.out' = NULL,
                           'ic.df' = NULL, 'mcmc.diagnostic' = NULL,
                           'p_calib_postsumm' = rbind(puncalib_ens,
                                                      puncalib_ens,
                                                      puncalib_ens))

        # posterior summary of calibrated estimate
        colnames(MCMCout_ens$p_calib) = causes
        rownames(MCMCout_ens$p_calib_postsumm) = c('postmean', 'lowcredI', 'upcredI')

      }else{

        #### stan fit ----
        stanfit_ens = rstan::sampling(seqcalib_mmat_model,
                                      pars = c('p_calib', 'loglik'),
                                      include = TRUE,
                                      data = list('nCause' = sum(!donotcalib_ens),
                                                  'nAlgo' = K,
                                                  'aj' = va_deaths[,!donotcalib_ens, drop = FALSE],
                                                  'p_uncalib' = puncalib_ens[!donotcalib_ens]/sum(puncalib_ens[!donotcalib_ens]),
                                                  'Mmat_bysimplex' = Mmat_asarray[,!donotcalib_ens,!donotcalib_ens, drop = FALSE],
                                                  'pss' = pss
                                      ),
                                      chains = 1,
                                      iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                                      control = list('adapt_delta' = adapt_delta_stan),
                                      seed = seed,
                                      refresh = refresh.stan#,
                                      # init = list(list('Mmat' = Mmat.init,
                                      #                  'p_calib' = puncalib))
        )

        MCMCout_ens = rstan::extract(stanfit_ens) # MCMC output
        MCMCout_ens$calibrated = calibrated_ens

        #### mcmc diagnostic ----
        # max Rhat
        max_Rhat_ens = max(apply(X = MCMCout_ens$p_calib, 2,
                                 FUN = function(v){

                                   rstan::Rhat(v)

                                 }))

        # min bulk ESS
        min_ess_bulk_ens = min(apply(X = MCMCout_ens$p_calib, 2,
                                     FUN = function(v){

                                       rstan::ess_bulk(v)

                                     }))/nMCMC

        # loo ic
        MCMCout_ens$loo.out = loo::loo(MCMCout_ens$loglik,
                                       r_eff = loo::relative_eff(exp(MCMCout_ens$loglik),
                                                                 chain_id = rep(1, nrow(MCMCout_ens$loglik))),
                                       cores = 1)

        # waic
        MCMCout_ens$waic.out = loo::waic(MCMCout_ens$loglik)

        ic.df_ens = rbind(MCMCout_ens$waic.out$estimates,
                          MCMCout_ens$loo.out$estimates)

        ic.df.melt_ens = as.numeric(t(ic.df_ens))
        names(ic.df.melt_ens) = paste0(rep(rownames(ic.df_ens), each = ncol(ic.df_ens)),'_',rep(colnames(ic.df_ens), nrow(ic.df_ens)))

        MCMCout_ens$ic.df = ic.df.melt_ens

        # mcmc diagnostic summary
        MCMCout_ens$mcmc.diagnostic = c('max_Rhat' = max_Rhat_ens, 'min_ess_bulk' = min_ess_bulk_ens,
                                        'num_divergent' = rstan::get_num_divergent(stanfit_ens),
                                        'num_max_treedepth' = rstan::get_num_max_treedepth(stanfit_ens))

        if(verbose){

          # cat("\n")
          message(MCMCout_ens$mcmc.diagnostic)
          # cat("\n")

        }

        #### calibrated posterior for all causes ----
        # print(head(MCMCout$p_calib))
        if(sum(donotcalib_ens)>0){

          # csmf
          p_calib_out_ens = matrix(data = puncalib_ens,
                                   nrow = nMCMC, ncol = nCause,
                                   byrow = TRUE)
          # print(dim(p_calib_out[,donotcalib_touse, drop = F]))
          # p_calib_out[,donotcalib_touse] = matrix(data = puncalib[donotcalib_touse],
          #                                         nrow = nMCMC, ncol = sum(donotcalib_touse),
          #                                         byrow = T)
          p_calib_out_ens[,!donotcalib_ens] = (1 - sum(puncalib_ens[donotcalib_ens]))*MCMCout_ens$p_calib
          MCMCout_ens$p_calib = p_calib_out_ens

        }
        colnames(MCMCout_ens$p_calib) = causes
        # MCMCout$p_calib.postmean = colMeans(MCMCout$p_calib)
        # print(head(MCMCout$p_calib))

        # posterior summary of calibrated estimate
        MCMCout_ens$p_calib_postsumm = apply(MCMCout_ens$p_calib, 2,
                                             FUN = function(v){

                                               c(mean(v),
                                                 quantile(x = v, probs = c(.025, .975)))

                                             })
        rownames(MCMCout_ens$p_calib_postsumm) = c('postmean', 'lowcredI', 'upcredI')

      }
      # print(1)

      #### calibrated number of deaths ----
      va_deaths_ens = va_deaths
      for(k in 1:K){

        va_deaths_ens_k = va_deaths[k,]

        if(calibrated_ens){

          # print(donotcalib_ens)
          # print(va_deaths[k,!donotcalib_ens])
          va_deaths_ens_k[!donotcalib_ens] = round(colMeans(sum(va_deaths[k,!donotcalib_ens])*
                                                              (MCMCout_ens$p_calib[,!donotcalib_ens]/(1 - sum(puncalib_ens[donotcalib_ens])))))

          # exactly match total number of deaths
          if(!isTRUE(all.equal(sum(va_deaths_ens_k), sum(va_deaths[k,])))){

            (adjust.amnt_ens_k = sum(va_deaths_ens_k) - sum(va_deaths[k,]))

            adjust.vec_ens_k = rep(0, nCause)
            (adjust.vec_ens_k[!donotcalib_ens] = rep(floor(abs(adjust.amnt_ens_k)/sum(!donotcalib_ens)),
                                                     sum(!donotcalib_ens)))

            id.extra.adjust_ens_k = head(intersect(order(va_deaths_ens_k, decreasing = TRUE),
                                                   which(!donotcalib_ens)),
                                         n = abs(adjust.amnt_ens_k) %% sum(!donotcalib_ens))
            adjust.vec_ens_k[id.extra.adjust_ens_k] = adjust.vec_ens_k[id.extra.adjust_ens_k] + 1
            adjust.vec_ens_k

            if(adjust.amnt_ens_k>0){

              va_deaths_ens_k = va_deaths_ens_k - adjust.vec_ens_k

            }else{

              va_deaths_ens_k = va_deaths_ens_k + adjust.vec_ens_k

            }

          }

        }

        va_deaths_ens[k,] = va_deaths_ens_k
        # print(k)

      }


      # ensemble calibration output
      calibout = c(calibout, list('ensemble' = MCMCout_ens))

      puncalib = rbind(puncalib, puncalib_ens)
      rownames(puncalib) = names(calibout)
      colnames(puncalib) = causes

      pcalib_asarray_wens = array(dim = c(K+1, nMCMC, nCause),
                                  dimnames = list(rownames(puncalib), NULL, causes))
      pcalib_asarray_wens[1:K,,] = pcalib_asarray
      pcalib_asarray_wens[K+1,,] = MCMCout_ens$p_calib

      pcalib_postsumm_asarray_wens = array(dim = c(K+1, 3, nCause),
                                           dimnames = list(rownames(puncalib),
                                                           dimnames(pcalib_postsumm_asarray)[[2]], causes))
      pcalib_postsumm_asarray_wens[1:K,,] = pcalib_postsumm_asarray
      pcalib_postsumm_asarray_wens[K+1,,] = MCMCout_ens$p_calib_postsumm

      donotcalib_wens = rbind(donotcalib_asmat, donotcalib_ens)
      rownames(donotcalib_wens) = names(calibout)
      colnames(donotcalib_wens) = causes


      # output list
      input.list.now = as.list(environment())
      input.list = input.list.now[names(input.list)]

      output = list("calib_MCMCout" = calibout,
                    "p_uncalib" = puncalib,
                    "p_calib" = pcalib_asarray_wens,
                    "pcalib_postsumm" = pcalib_postsumm_asarray_wens,
                    "va_deaths_uncalib" = va_deaths,
                    "va_deaths_calib_algo" = va_deaths_calib,
                    "va_deaths_calib_ensemble" = va_deaths_ens,
                    'Mmat.fixed_input' = Mmat_input_asarray,
                    'Mmat.fixed' = Mmat_asarray,
                    'donotcalib' = donotcalib_wens,
                    'causes_notcalibrated' = c(causes_notcalibrated,
                                               list("ensemble" = causes[donotcalib_ens])),
                    'input' = input.list)

      # print(dim(Mmat_input_asarray))
      # print(dim(Mmat_asarray))
      # print(dim(output$Mmat.fixed_input))
      # print(dim(output$Mmat.fixed))

      # end ensemble

    }

    # ending Mmat fixed

  }else if(calibmodel.type=="Mmatprior") {

    ### prior Mmat ----
    seqcalib_model = readRDS(system.file("stan/seqcalib.rds", package = "vacalibration"))

    #### hard checks ----
    if(!is.list(Mmat.asDirich)) stop("'Mmat.asDirich' must be a list of matrices for each algorithm. Rows and columns of each matrix must correspond to CHAMPS and VA causes.")
    if(is.null(names(Mmat.asDirich))) stop("Components of 'Mmat.asDirich' must have algorithms as names.")
    if(!identical(names(Mmat.asDirich), dimnames(va_deaths)[[1]])) stop("Algorithm names specified in 'Mmat.asDirich' and 'va_unlabeled' do not match.")

    #### causes not calibrated ----
    if(is.null(donotcalib)){

      donotcalib = lapply(1:K, FUN = function(k){NULL})
      names(donotcalib) = dimnames(va_deaths)[[1]]

    }else{

      if(!is.list(donotcalib)){

        if(is.character(donotcalib)){

          donotcalib = lapply(1:K, FUN = function(k){donotcalib})
          names(donotcalib) = dimnames(va_deaths)[[1]]

        }

        # stop("'donotcalib' must be a list of character vectors for each algorithm. Each vector indicates the causes not to calibrate.")

      }else{

        if(is.null(names(donotcalib))){

          stop("Components of 'donotcalib' must have algorithms as names.")

        }else if(!identical(names(donotcalib), dimnames(va_deaths)[[1]])){

          stop("Either names are not provided, or they do not match with that provided in 'va_unlabeled'")

        }

      }

    }

    if(!identical(names(donotcalib), dimnames(va_deaths)[[1]])) stop("Algorithm names specified in 'donotcalib' and 'va_unlabeled' do not match.")



    #### calibration for each algorithm ----

    # storage
    Mmat_input_asarray = array(dim = c(K, nCause, nCause),
                               dimnames = list(names(Mmat.asDirich), causes, causes))
    Mmat_asarray = array(dim = c(K, nCause, nCause),
                         dimnames = list(names(Mmat.asDirich), causes, causes))
    donotcalib_asmat = matrix(nrow = K, ncol = nCause)
    causes_notcalibrated = calibout = vector(mode = "list", length = K)
    rownames(donotcalib_asmat) = names(causes_notcalibrated) =
      names(calibout) = names(Mmat.asDirich)
    colnames(donotcalib_asmat) = causes
    va_deaths_calib = va_deaths
    pcalib_asarray = array(dim = c(K, nMCMC, nCause),
                           dimnames = list(names(Mmat.asDirich), NULL, causes))
    pcalib_postsumm_asarray = array(dim = c(K, 3, nCause),
                                    dimnames = list(names(Mmat.asDirich), c('postmean', 'lowcredI', 'upcredI'), causes))

    for(k in 1:K){

      # hard checks
      if(!is.matrix(Mmat.asDirich[[k]])) stop(paste0("'Mmat.asDirich' for algorithm ", names(Mmat.asDirich)[k]," must be a matrix with CHAMPS causes along the rows and VA causes along the columns."))
      if(is.null(colnames(Mmat.asDirich[[k]]))) stop(paste0("'Mmat.asDirich' for algorithm ", names(Mmat.asDirich)[k], " must have VA causes as column names."))
      if(!identical(colnames(Mmat.asDirich[[k]]), causes)) stop(paste0("VA causes specified in 'Mmat.asDirich' for algorithm ", names(Mmat.asDirich)[k], " do not match with that specified in 'va_unlabeled'."))
      if(is.null(rownames(Mmat.asDirich[[k]]))) stop(paste0("'Mmat.asDirich' for algorithm ", names(Mmat.asDirich)[k], " must have CHAMPS causes as row names."))
      if(!identical(rownames(Mmat.asDirich[[k]]), causes)) stop(paste0("CHAMPS causes specified in 'Mmat.asDirich' for algorithm ", names(Mmat.asDirich)[k], " do not match with that specified in 'va_unlabeled'."))
      # if(any(rowSums(Mmat.asDirich)!=1)) stop(paste0("Each row of 'Mmat.asDirich' for algorithm ", names(Mmat.asDirich)[k], " must sum to 1."))

      (Mmat.meanDirich_k = Mmat.asDirich[[k]]/rowSums(Mmat.asDirich[[k]])) # prior mean

      # causes not to calibrate
      if((!is.null(donotcalib[[k]])) & (sum(causes %in% donotcalib[[k]])==0)){

        stop("'donotcalib' for algorithm ", names(donotcalib)[k],"  doesn't match with causes named in 'va_unlabeled' and 'Mmat.asDirich'")

      }
      # if(('other' %in% causes) & (!('other' %in% donotcalib[[k]]))) donotcalib[[k]] = c(donotcalib[[k]], 'other')

      (donotcalib_asmat[k,] = causes %in% donotcalib[[k]])

      # learn more causes not to calibrate
      if(donot.calib_type=="learn"){

        (donot.calib_Mmat_k = apply(Mmat.meanDirich_k, 2,
                                    FUN = function(v){

                                      # mean(abs(v - mean(v)))
                                      diff(range(v))<=nocalib.threshold

                                    }))
        (donotcalib_asmat[k,] = donotcalib_asmat[k,]|donot.calib_Mmat_k)

      }
      causes_notcalibrated[[k]] = causes[donotcalib_asmat[k,]]

      ##### preparing Mmat for calibration ----
      if(sum(!donotcalib_asmat[k,])<=1){

        calibrated_k = FALSE
        lambda_k = 1
        Mmat_input_asarray[k,,] = Mmat_asarray[k,,] = diag(rowSums(Mmat.asDirich[[k]]))

      }else{

        calibrated_k = TRUE
        Mmat_input_asarray[k,,] = Mmat_asarray[k,,] = Mmat.asDirich[[k]]

        idtocalib_k = which(!donotcalib_asmat[k,])
        nTocalib_k = length(idtocalib_k)

        # shrink Mmat for stability
        if(!stable){

          lambda_k = 0

        }else{

          # sub Mmat for calibrates causes
          Mmat_fixed_sub_k = Mmat.asDirich[[k]][!donotcalib_asmat[k,],!donotcalib_asmat[k,]]
          Mmat_fixed_sub_k = Mmat_fixed_sub_k/rowSums(Mmat_fixed_sub_k)

          # print(k)
          # print(puncalib)
          # print(donotcalib_asmat)
          puncalib_sub_k = puncalib[k,!donotcalib_asmat[k,]]/sum(puncalib[k,!donotcalib_asmat[k,]])

          count.lambda_k = 0
          got_it_k = FALSE
          while(!got_it_k){

            # shrink coeff
            count.lambda_k = count.lambda_k + 1
            if(count.lambda_k==1){

              lambda_k = .99

            }else{

              lambda_k = lambda_k - .01

            }

            # shrinking towards identity
            (Mmat_lambda_k = lambda_k*diag(nrow(Mmat_fixed_sub_k)) + (1-lambda_k)*Mmat_fixed_sub_k)

            # solving calibration eq
            pcalib_lambda_k = as.numeric(solve(Mmat_lambda_k %*% t(Mmat_lambda_k)) %*%
                                           Mmat_lambda_k %*% as.matrix(puncalib_sub_k))

            # pcalib_lambda = puncalib
            # pcalib_lambda[!donotcalib_asmat[k,]] = (1-sum(puncalib[donotcalib_asmat[k,]]))*pcalib_lambda

            # checking if within simplex
            within_simplex_k = (sum(pcalib_lambda_k<0)==0)&all.equal(sum(pcalib_lambda_k), 1)
            got_it_k = (!within_simplex_k)|isTRUE(all.equal(lambda_k, 0))
            # got_it = !within_simplex

          }

          if((lambda_k<.99) & (!within_simplex_k)){

            lambda_k = lambda_k + .01

          }

          # shrinked Mmat
          (Mmat_lambda_opt_k = lambda_k*diag(nrow(Mmat_fixed_sub_k)) + (1-lambda_k)*Mmat_fixed_sub_k)
          rownames(Mmat_lambda_opt_k) = colnames(Mmat_lambda_opt_k) = causes[!donotcalib_asmat[k,]]

          Mmat_asarray[k,!donotcalib_asmat[k,],!donotcalib_asmat[k,]] =
            rowSums(Mmat.asDirich[[k]][!donotcalib_asmat[k,],!donotcalib_asmat[k,]])*Mmat_lambda_opt_k

        }

      }



      #### calibration output ----
      if(!calibrated_k){

        ###### not calibrated ----
        MCMCout_k = list('p_calib' = matrix(data = puncalib[k,],
                                            nrow = nMCMC, ncol = nCause,
                                            byrow = TRUE),
                         'loglik' = NULL,
                         "calibrated" = calibrated_k, "lambda" = lambda_k,
                         'loo.out' = NULL, 'waic.out' = NULL,
                         'ic.df' = NULL, 'mcmc.diagnostic' = NULL,
                         'p_calib_postsumm' = rbind(puncalib[k,],
                                                    puncalib[k,],
                                                    puncalib[k,]))

        # posterior summary of calibrated estimate
        colnames(MCMCout_k$p_calib) = causes
        rownames(MCMCout_k$p_calib_postsumm) = c('postmean', 'lowcredI', 'upcredI')

      }else{

        # print(nTocalib_k)
        # print(idtocalib_k)

        ###### calibrated ----
        stanfit_k = rstan::sampling(seqcalib_model,
                                    pars = c('Mmat', 'p_calib', 'loglik'),
                                    include = TRUE,
                                    data = list('nCause' = nCause,
                                                'nAlgo' = 1,
                                                'aj' = va_deaths[k,,drop = FALSE],
                                                'nTocalib' = as.array(nTocalib_k),
                                                'idtocalib' = idtocalib_k,
                                                'nCumtocalib' = c(0, cumsum(nTocalib_k)),
                                                'p_uncalib' = puncalib[k,],
                                                'Mmatprior_asDirich' = Mmat_asarray[k,,,drop = FALSE],
                                                'pss' = pss
                                    ),
                                    chains = 1,
                                    iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                                    control = list('adapt_delta' = adapt_delta_stan),
                                    seed = seed,
                                    refresh = refresh.stan#,
                                    # init = list(list('Mmat' = Mmat.init,
                                    #                  'p_calib' = puncalib))
        )

        MCMCout_k = rstan::extract(stanfit_k) # MCMC output
        MCMCout_k$calibrated = calibrated_k
        MCMCout_k$lambda = lambda_k

        # mcmc diagnostic
        # max Rhat
        max_Rhat_k = max(apply(X = MCMCout_k$p_calib, 2,
                               FUN = function(v){

                                 rstan::Rhat(v)

                               }))

        # min bulk ESS
        min_ess_bulk_k = min(apply(X = MCMCout_k$p_calib, 2,
                                   FUN = function(v){

                                     rstan::ess_bulk(v)

                                   }))/nMCMC

        # loo ic
        MCMCout_k$loo.out = loo::loo(MCMCout_k$loglik,
                                     r_eff = loo::relative_eff(exp(MCMCout_k$loglik),
                                                               chain_id = rep(1, nrow(MCMCout_k$loglik))),
                                     cores = 1)

        # waic
        MCMCout_k$waic.out = loo::waic(MCMCout_k$loglik)

        ic.df_k = rbind(MCMCout_k$waic.out$estimates,
                        MCMCout_k$loo.out$estimates)

        ic.df.melt_k = as.numeric(t(ic.df_k))
        names(ic.df.melt_k) = paste0(rep(rownames(ic.df_k), each = ncol(ic.df_k)),'_',rep(colnames(ic.df_k), nrow(ic.df_k)))

        MCMCout_k$ic.df = ic.df.melt_k

        # mcmc diagnostic summary
        MCMCout_k$mcmc.diagnostic = c('max_Rhat' = max_Rhat_k, 'min_ess_bulk' = min_ess_bulk_k,
                                      'num_divergent' = rstan::get_num_divergent(stanfit_k),
                                      'num_max_treedepth' = rstan::get_num_max_treedepth(stanfit_k))

        if(verbose){

          # cat("\n")
          message(MCMCout_k$mcmc.diagnostic)
          # cat("\n")

        }

        # calibrated posterior for all causes
        # print(head(MCMCout$p_calib))
        if(sum(donotcalib_asmat[k,])>0){

          # csmf
          p_calib_out_k = matrix(data = puncalib[k,],
                                 nrow = nMCMC, ncol = nCause,
                                 byrow = TRUE)
          # print(dim(p_calib_out[,donotcalib_touse, drop = F]))
          # p_calib_out[,donotcalib_touse] = matrix(data = puncalib[donotcalib_touse],
          #                                         nrow = nMCMC, ncol = sum(donotcalib_touse),
          #                                         byrow = T)
          p_calib_out_k[,!donotcalib_asmat[k,]] = (1 - sum(puncalib[k,donotcalib_asmat[k,]]))*
            (MCMCout_k$p_calib[,!donotcalib_asmat[k,]]/rowSums(MCMCout_k$p_calib[,!donotcalib_asmat[k,]]))
          MCMCout_k$p_calib = p_calib_out_k

        }
        colnames(MCMCout_k$p_calib) = causes
        # MCMCout$p_calib.postmean = colMeans(MCMCout$p_calib)
        # print(head(MCMCout_k$p_calib))

        # posterior summary of calibrated estimate
        MCMCout_k$p_calib_postsumm = apply(MCMCout_k$p_calib, 2,
                                           FUN = function(v){

                                             c(mean(v),
                                               quantile(x = v, probs = c(.025, .975)))

                                           })
        rownames(MCMCout_k$p_calib_postsumm) = c('postmean', 'lowcredI', 'upcredI')

      }

      # store MCMC output
      calibout[[k]] = MCMCout_k

      #### calibrated number of deaths ----
      va_deaths_calib_k = va_deaths[k,]
      if(calibrated_k){

        va_deaths_calib_k[!donotcalib_asmat[k,]] = round(colMeans(sum(va_deaths[k,!donotcalib_asmat[k,]])*
                                                                    (MCMCout_k$p_calib[,!donotcalib_asmat[k,]]/(1 - sum(puncalib[k,donotcalib_asmat[k,]])))))

        # exactly match total number of deaths
        if(!isTRUE(all.equal(sum(va_deaths_calib_k), sum(va_deaths[k,])))){

          (adjust.amnt_k = sum(va_deaths_calib_k) - sum(va_deaths[k,]))

          adjust.vec_k = rep(0, nCause)
          (adjust.vec_k[!donotcalib_asmat[k,]] = rep(floor(abs(adjust.amnt_k)/sum(!donotcalib_asmat[k,])),
                                                     sum(!donotcalib_asmat[k,])))

          id.extra.adjust_k = head(intersect(order(va_deaths_calib_k, decreasing = TRUE),
                                             which(!donotcalib_asmat[k,])),
                                   n = abs(adjust.amnt_k) %% sum(!donotcalib_asmat[k,]))
          adjust.vec_k[id.extra.adjust_k] = adjust.vec_k[id.extra.adjust_k] + 1
          adjust.vec_k

          if(adjust.amnt_k>0){

            va_deaths_calib_k = va_deaths_calib_k - adjust.vec_k

          }else{

            va_deaths_calib_k = va_deaths_calib_k + adjust.vec_k

          }

        }

      }

      va_deaths_calib[k,] = va_deaths_calib_k
      pcalib_asarray[k,,] = MCMCout_k$p_calib
      pcalib_postsumm_asarray[k,,] = MCMCout_k$p_calib_postsumm

    }# ending calibration for each algorithm


    ### ensemble calibration  ----
    if(K==1){

      if(is.null(ensemble)){

        ensemble = FALSE

      }else{

        if(ensemble){

          ensemble = FALSE
          if(verbose) message("Nothing to ensemble. Only one algorithm provided. Setting ensemble=F.")

        }

      }

    }else if(K>1){

      # defaulting to ensemble if multiple algorithms
      if(is.null(ensemble)) ensemble = TRUE

    }

    # checking if needs to run ensemble
    if(!ensemble){

      # output list
      input.list.now = as.list(environment())
      input.list = c(input.list.now[names(input.list)],
                     "K" = K, "nCause" = nCause, "causes" = causes)

      output = list("calib_MCMCout" = calibout,
                    "p_uncalib" = puncalib,
                    "p_calib" = pcalib_asarray,
                    "pcalib_postsumm" = pcalib_postsumm_asarray,
                    "va_deaths_uncalib" = va_deaths,
                    "va_deaths_calib_algo" = va_deaths_calib,
                    'Mmat.asDirich_input' = Mmat_input_asarray,
                    'Mmat.asDirich' = Mmat_asarray,
                    'donotcalib' = donotcalib_asmat,
                    'causes_notcalibrated' = causes_notcalibrated,
                    'input' = input.list)

    }else{

      # ensemble uncalibrated csmf
      (puncalib_ens = colSums(va_deaths)/sum(va_deaths))

      donotcalib_ens = donotcalib_asmat[1,]
      for(k in 2:K){

        donotcalib_ens = donotcalib_ens&donotcalib_asmat[k,]

      }

      idtocalib_ens = nTocalib_ens = NULL
      for(k in 1:K){

        idtocalib_ens_k = which(!donotcalib_asmat[k,])
        idtocalib_ens = c(idtocalib_ens, idtocalib_ens_k)
        nTocalib_ens = c(nTocalib_ens, length(idtocalib_ens_k))

      }

      calibrated_ens = sum(!donotcalib_ens)>1
      if(!calibrated_ens){

        MCMCout_ens = list('p_calib' = matrix(data = puncalib_ens,
                                              nrow = nMCMC, ncol = nCause,
                                              byrow = TRUE),
                           'loglik' = NULL,
                           "calibrated" = calibrated_ens,
                           'loo.out' = NULL, 'waic.out' = NULL,
                           'ic.df' = NULL, 'mcmc.diagnostic' = NULL,
                           'p_calib_postsumm' = rbind(puncalib_ens,
                                                      puncalib_ens,
                                                      puncalib_ens))

        # posterior summary of calibrated estimate
        colnames(MCMCout_ens$p_calib) = causes
        rownames(MCMCout_ens$p_calib_postsumm) = c('postmean', 'lowcredI', 'upcredI')

      }else{

        #### stan fit ----
        stanfit_ens = rstan::sampling(seqcalib_model,
                                      pars = c('Mmat', 'p_calib', 'loglik'),
                                      include = TRUE,
                                      data = list('nCause' = nCause,
                                                  'nAlgo' = K,
                                                  'aj' = va_deaths,
                                                  'nTocalib' = nTocalib_ens,
                                                  'idtocalib' = idtocalib_ens,
                                                  'nCumtocalib' = c(0, cumsum(nTocalib_ens)),
                                                  'p_uncalib' = puncalib_ens,
                                                  'Mmatprior_asDirich' = Mmat_asarray,
                                                  'pss' = pss
                                      ),
                                      chains = 1,
                                      iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                                      control = list('adapt_delta' = adapt_delta_stan),
                                      seed = seed,
                                      refresh = refresh.stan#,
                                      # init = list(list('Mmat' = Mmat.init,
                                      #                  'p_calib' = puncalib))
        )

        MCMCout_ens = rstan::extract(stanfit_ens) # MCMC output
        MCMCout_ens$calibrated = calibrated_ens

        #### mcmc diagnostic ----
        # max Rhat
        max_Rhat_ens = max(apply(X = MCMCout_ens$p_calib, 2,
                                 FUN = function(v){

                                   rstan::Rhat(v)

                                 }))

        # min bulk ESS
        min_ess_bulk_ens = min(apply(X = MCMCout_ens$p_calib, 2,
                                     FUN = function(v){

                                       rstan::ess_bulk(v)

                                     }))/nMCMC

        # loo ic
        MCMCout_ens$loo.out = loo::loo(MCMCout_ens$loglik,
                                       r_eff = loo::relative_eff(exp(MCMCout_ens$loglik),
                                                                 chain_id = rep(1, nrow(MCMCout_ens$loglik))),
                                       cores = 1)

        # waic
        MCMCout_ens$waic.out = loo::waic(MCMCout_ens$loglik)

        ic.df_ens = rbind(MCMCout_ens$waic.out$estimates,
                          MCMCout_ens$loo.out$estimates)

        ic.df.melt_ens = as.numeric(t(ic.df_ens))
        names(ic.df.melt_ens) = paste0(rep(rownames(ic.df_ens), each = ncol(ic.df_ens)),'_',rep(colnames(ic.df_ens), nrow(ic.df_ens)))

        MCMCout_ens$ic.df = ic.df.melt_ens

        # mcmc diagnostic summary
        MCMCout_ens$mcmc.diagnostic = c('max_Rhat' = max_Rhat_ens, 'min_ess_bulk' = min_ess_bulk_ens,
                                        'num_divergent' = rstan::get_num_divergent(stanfit_ens),
                                        'num_max_treedepth' = rstan::get_num_max_treedepth(stanfit_ens))

        if(verbose){

          # cat("\n")
          message(MCMCout_ens$mcmc.diagnostic)
          # cat("\n")

        }

        #### calibrated posterior for all causes ----
        # print(head(MCMCout$p_calib))
        if(sum(donotcalib_ens)>0){

          # csmf
          p_calib_out_ens = matrix(data = puncalib_ens,
                                   nrow = nMCMC, ncol = nCause,
                                   byrow = TRUE)
          # print(dim(p_calib_out[,donotcalib_touse, drop = F]))
          # p_calib_out[,donotcalib_touse] = matrix(data = puncalib[donotcalib_touse],
          #                                         nrow = nMCMC, ncol = sum(donotcalib_touse),
          #                                         byrow = T)
          p_calib_out_ens[,!donotcalib_ens] = (1 - sum(puncalib_ens[donotcalib_ens]))*
            (MCMCout_ens$p_calib[,!donotcalib_ens]/rowSums(MCMCout_ens$p_calib[,!donotcalib_ens]))
          MCMCout_ens$p_calib = p_calib_out_ens

        }
        colnames(MCMCout_ens$p_calib) = causes
        # MCMCout$p_calib.postmean = colMeans(MCMCout$p_calib)
        # print(head(MCMCout$p_calib))

        # posterior summary of calibrated estimate
        MCMCout_ens$p_calib_postsumm = apply(MCMCout_ens$p_calib, 2,
                                             FUN = function(v){

                                               c(mean(v),
                                                 quantile(x = v, probs = c(.025, .975)))

                                             })
        rownames(MCMCout_ens$p_calib_postsumm) = c('postmean', 'lowcredI', 'upcredI')

      }
      # print(1)

      #### calibrated number of deaths ----
      va_deaths_ens = va_deaths
      for(k in 1:K){

        va_deaths_ens_k = va_deaths[k,]

        if(calibrated_ens){

          # print(donotcalib_ens)
          # print(va_deaths[k,!donotcalib_ens])
          va_deaths_ens_k[!donotcalib_ens] = round(colMeans(sum(va_deaths[k,!donotcalib_ens])*
                                                              (MCMCout_ens$p_calib[,!donotcalib_ens]/(1 - sum(puncalib_ens[donotcalib_ens])))))

          # exactly match total number of deaths
          if(!isTRUE(all.equal(sum(va_deaths_ens_k), sum(va_deaths[k,])))){

            (adjust.amnt_ens_k = sum(va_deaths_ens_k) - sum(va_deaths[k,]))

            adjust.vec_ens_k = rep(0, nCause)
            (adjust.vec_ens_k[!donotcalib_ens] = rep(floor(abs(adjust.amnt_ens_k)/sum(!donotcalib_ens)),
                                                     sum(!donotcalib_ens)))

            id.extra.adjust_ens_k = head(intersect(order(va_deaths_ens_k, decreasing = TRUE),
                                                   which(!donotcalib_ens)),
                                         n = abs(adjust.amnt_ens_k) %% sum(!donotcalib_ens))
            adjust.vec_ens_k[id.extra.adjust_ens_k] = adjust.vec_ens_k[id.extra.adjust_ens_k] + 1
            adjust.vec_ens_k

            if(adjust.amnt_ens_k>0){

              va_deaths_ens_k = va_deaths_ens_k - adjust.vec_ens_k

            }else{

              va_deaths_ens_k = va_deaths_ens_k + adjust.vec_ens_k

            }

          }

        }

        va_deaths_ens[k,] = va_deaths_ens_k
        # print(k)

      }


      # ensemble calibration output
      calibout = c(calibout, list('ensemble' = MCMCout_ens))

      puncalib = rbind(puncalib, puncalib_ens)
      rownames(puncalib) = names(calibout)
      colnames(puncalib) = causes

      pcalib_asarray_wens = array(dim = c(K+1, nMCMC, nCause),
                                  dimnames = list(rownames(puncalib), NULL, causes))
      pcalib_asarray_wens[1:K,,] = pcalib_asarray
      pcalib_asarray_wens[K+1,,] = MCMCout_ens$p_calib

      pcalib_postsumm_asarray_wens = array(dim = c(K+1, 3, nCause),
                                           dimnames = list(rownames(puncalib),
                                                           dimnames(pcalib_postsumm_asarray)[[2]], causes))
      pcalib_postsumm_asarray_wens[1:K,,] = pcalib_postsumm_asarray
      pcalib_postsumm_asarray_wens[K+1,,] = MCMCout_ens$p_calib_postsumm

      donotcalib_wens = rbind(donotcalib_asmat, donotcalib_ens)
      rownames(donotcalib_wens) = names(calibout)
      colnames(donotcalib_wens) = causes


      # output list
      input.list.now = as.list(environment())
      input.list = input.list.now[names(input.list)]

      output = list("calib_MCMCout" = calibout,
                    "p_uncalib" = puncalib,
                    "p_calib" = pcalib_asarray_wens,
                    "pcalib_postsumm" = pcalib_postsumm_asarray_wens,
                    "va_deaths_uncalib" = va_deaths,
                    "va_deaths_calib_algo" = va_deaths_calib,
                    "va_deaths_calib_ensemble" = va_deaths_ens,
                    'Mmat.asDirich_input' = Mmat_input_asarray,
                    'Mmat.asDirich' = Mmat_asarray,
                    'donotcalib' = donotcalib_wens,
                    'causes_notcalibrated' = c(causes_notcalibrated,
                                               list("ensemble" = causes[donotcalib_ens])),
                    'input' = input.list)

      # print(dim(Mmat_input_asarray))
      # print(dim(Mmat_asarray))
      # print(dim(output$Mmat.asDirich_input))
      # print(dim(output$Mmat.asDirich))

      # end ensemble

    }

    # ending Mmat fixed

  }

  # print("calibrated succesfully")

  ## comparison plots ----
  if(plot_it){

    # heatmap of misclassification
    if(calibmodel.type=="Mmatfixed"){


      title_Mmat_input = title_Mmat = "Fixed Misclassification Matrix"

      if(stable){

        title_Mmat_input = paste0(title_Mmat_input, " (Estimated from CHAMPS Data)")
        title_Mmat = paste0(title_Mmat, " (Adjusted)")

      }


      # Mmat input
      # print(round(100*output$Mmat.fixed_input[1,,]))
      Mmat_input_toplot = round(100*output$Mmat.fixed_input)
      # print(dim(Mmat_input_toplot))
      for(k in 1:dim(Mmat_input_toplot)[1]){

        for(i in 1:dim(Mmat_input_toplot)[2]){

          if(!any(is.nan(Mmat_input_toplot[k,i,]))){

            nonzerocauseid_ki = which.max(Mmat_input_toplot[k,i,])
            Mmat_input_toplot[k,i,nonzerocauseid_ki] = 100 - sum(Mmat_input_toplot[k,i,-nonzerocauseid_ki])

          }

        }

        Mmat_input_toplot[k,donotcalib_asmat[k,],] =
          Mmat_input_toplot[k,,donotcalib_asmat[k,]] = NA

      }

      plotdf_Mmat_input = reshape2::melt(Mmat_input_toplot)
      head(plotdf_Mmat_input)
      value.labels_Mmat_input = plotdf_Mmat_input$value

      plotdf_Mmat_input$Var1 = factor(x = plotdf_Mmat_input$Var1, levels = dimnames(Mmat_input_toplot)[[1]])
      plotdf_Mmat_input$Var2 = factor(x = plotdf_Mmat_input$Var2, levels = rev(causes))
      plotdf_Mmat_input$Var3 = factor(x = plotdf_Mmat_input$Var3, levels = causes)

      plotdf_Mmat_input$diag = plotdf_Mmat_input$Var2==plotdf_Mmat_input$Var3
      plotdf_Mmat_input$diag[!plotdf_Mmat_input$diag] = NA
      # print(head(plotdf_Mmat_input))
      # print(K)

      if(K>1){

        ggplot2_Mmat_input =
          ggplot2::ggplot(plotdf_Mmat_input, ggplot2::aes(Var3, Var2, fill = value)) +
          ggplot2::geom_tile(color="white", linewidth=.5) +
          ggplot2::geom_tile(data = plotdf_Mmat_input[!is.na(plotdf_Mmat_input$diag), ],
                             ggplot2::aes(Var3, Var2, fill = value, color = diag), linewidth = .7) +
          ggplot2::scale_color_manual(guide = "none", values = c(`TRUE` = "blue3")) +
          ggplot2::geom_text(ggplot2::aes(Var3, Var2, label = value.labels_Mmat_input,
                                          size = value^(1/2)),
                             color = "black",
                             fontface = 'bold') +
          # ggplot2::scale_size(range = c(0, 8)) +
          ggplot2::scale_fill_gradient(low="white", high="red3",
                                       breaks = seq(0, 100, 20), limits = c(0,100),
                                       name = 'Classification Percentage') +
          ggplot2::facet_grid(.~Var1) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            # axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text.x = ggplot2::element_text(color = "black",
                                                angle = 30, hjust = 1, vjust = 1),
            axis.text.y = ggplot2::element_text(color = "black"),
            # axis.ticks.x = ggplot2::element_blank(),
            # axis.ticks.length.x = unit(.2, "cm"),
            # axis.ticks.y = ggplot2::element_line(linewidth = .5),
            # axis.ticks.length.y = ggplot2::unit(.2, "cm"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 1),
            panel.grid.major = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            panel.grid.minor = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            strip.text.x = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.text.y = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.background = ggplot2::element_rect(color="black", linewidth=1),
            legend.title = ggplot2::element_blank(),
            # legend.key.width = ggplot2::unit(.75, "cm"),
            # legend.key.height = ggplot2::unit(.75, "cm"),
            # legend.key.size = ggplot2::unit(.5, "cm"),
            # legend.spacing.x = ggplot2::unit(.5, 'cm'),
            # legend.text=ggplot2::element_text(size=22),
            legend.text=ggplot2::element_text(size=12),
            legend.position = 'none'
          ) +
          ggplot2::labs(title = title_Mmat_input,
                        x = 'VA Cause', y = 'CHAMPS Cause')
        # print(ggplot2_Mmat_input)

      }else{

        ggplot2_Mmat_input =
          ggplot2::ggplot(plotdf_Mmat_input, ggplot2::aes(Var3, Var2, fill = value)) +
          ggplot2::geom_tile(color="white", linewidth=.5) +
          ggplot2::geom_tile(data = plotdf_Mmat_input[!is.na(plotdf_Mmat_input$diag), ],
                             ggplot2::aes(Var3, Var2, fill = value, color = diag), linewidth = .7) +
          ggplot2::scale_color_manual(guide = "none", values = c(`TRUE` = "blue3")) +
          ggplot2::geom_text(ggplot2::aes(Var3, Var2, label = value.labels_Mmat_input,
                                          size = value^(1/2)),
                             color = "black",
                             fontface = 'bold') +
          # ggplot2::scale_size(range = c(0, 8)) +
          ggplot2::scale_fill_gradient(low="white", high="red3",
                                       breaks = seq(0, 100, 20), limits = c(0,100),
                                       name = 'Classification Percentage') +
          ggplot2::facet_grid(.~Var1) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            # axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text.x = ggplot2::element_text(color = "black",
                                                angle = 30, hjust = 1, vjust = 1),
            axis.text.y = ggplot2::element_text(color = "black"),
            # axis.ticks.x = ggplot2::element_blank(),
            # axis.ticks.length.x = unit(.2, "cm"),
            # axis.ticks.y = ggplot2::element_line(linewidth = .5),
            # axis.ticks.length.y = ggplot2::unit(.2, "cm"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 1),
            panel.grid.major = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            panel.grid.minor = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            strip.text.x = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.text.y = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.background = ggplot2::element_rect(color="black", linewidth=1),
            legend.title = ggplot2::element_blank(),
            # legend.key.width = ggplot2::unit(.75, "cm"),
            # legend.key.height = ggplot2::unit(.75, "cm"),
            # legend.key.size = ggplot2::unit(.5, "cm"),
            # legend.spacing.x = ggplot2::unit(.5, 'cm'),
            # legend.text=ggplot2::element_text(size=22),
            legend.text=ggplot2::element_text(size=12),
            legend.position = 'none'
          ) +
          ggplot2::labs(title = title_Mmat_input,
                        x = 'VA Cause', y = 'CHAMPS Cause')
        # print(ggplot2_Mmat_input)

      }


      # Mmat used for calibration
      # print(round(100*output$Mmat.fixed[1,,]))
      (Mmat_toplot = round(100*output$Mmat.fixed))
      for(k in 1:dim(Mmat_toplot)[1]){

        for(i in 1:dim(Mmat_toplot)[2]){

          if(!any(is.nan(Mmat_toplot[k,i,]))){

            nonzerocauseid_ki = which.max(Mmat_toplot[k,i,])
            Mmat_toplot[k,i,nonzerocauseid_ki] = 100 - sum(Mmat_toplot[k,i,-nonzerocauseid_ki])

          }

        }

        Mmat_toplot[k,donotcalib_asmat[k,],] =
          Mmat_toplot[k,,donotcalib_asmat[k,]] = NA

      }

      plotdf_Mmat = reshape2::melt(Mmat_toplot)
      head(plotdf_Mmat)
      value.labels_Mmat = plotdf_Mmat$value

      plotdf_Mmat$Var1 = factor(x = plotdf_Mmat$Var1, levels = dimnames(Mmat_toplot)[[1]])
      plotdf_Mmat$Var2 = factor(x = plotdf_Mmat$Var2, levels = rev(causes))
      plotdf_Mmat$Var3 = factor(x = plotdf_Mmat$Var3, levels = causes)

      plotdf_Mmat$diag = plotdf_Mmat$Var2==plotdf_Mmat$Var3
      plotdf_Mmat$diag[!plotdf_Mmat$diag] = NA
      head(plotdf_Mmat)

      if(K>1){

        ggplot2_Mmat =
          ggplot2::ggplot(plotdf_Mmat, ggplot2::aes(Var3, Var2, fill = value)) +
          ggplot2::geom_tile(color="white", linewidth=.5) +
          ggplot2::geom_tile(data = plotdf_Mmat[!is.na(plotdf_Mmat$diag), ],
                             ggplot2::aes(Var3, Var2, fill = value, color = diag), linewidth = .7) +
          ggplot2::scale_color_manual(guide = "none", values = c(`TRUE` = "blue3")) +
          ggplot2::geom_text(ggplot2::aes(Var3, Var2, label = value.labels_Mmat,
                                          size = value^(1/2)),
                             color = "black",
                             fontface = 'bold') +
          # ggplot2::scale_size(range = c(0, 8)) +
          ggplot2::scale_fill_gradient(low="white", high="red3",
                                       breaks = seq(0, 100, 20), limits = c(0,100),
                                       name = 'Classification Percentage') +
          ggplot2::facet_grid(.~Var1) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            # axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text.x = ggplot2::element_text(color = "black",
                                                angle = 30, hjust = 1, vjust = 1),
            axis.text.y = ggplot2::element_text(color = "black"),
            # axis.ticks.x = ggplot2::element_blank(),
            # axis.ticks.length.x = unit(.2, "cm"),
            # axis.ticks.y = ggplot2::element_line(linewidth = .5),
            # axis.ticks.length.y = ggplot2::unit(.2, "cm"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 1),
            panel.grid.major = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            panel.grid.minor = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            strip.text.x = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.text.y = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.background = ggplot2::element_rect(color="black", linewidth=1),
            legend.title = ggplot2::element_blank(),
            # legend.key.width = ggplot2::unit(.75, "cm"),
            # legend.key.height = ggplot2::unit(.75, "cm"),
            # legend.key.size = ggplot2::unit(.5, "cm"),
            # legend.spacing.x = ggplot2::unit(.5, 'cm'),
            # legend.text=ggplot2::element_text(size=22),
            legend.text=ggplot2::element_text(size=12),
            legend.position = 'none'
          ) +
          ggplot2::labs(title = title_Mmat,
                        x = 'VA Cause', y = 'CHAMPS Cause')
        # print(ggplot2_Mmat)

      }else{

        ggplot2_Mmat =
          ggplot2::ggplot(plotdf_Mmat, ggplot2::aes(Var3, Var2, fill = value)) +
          ggplot2::geom_tile(color="white", linewidth=.5) +
          ggplot2::geom_tile(data = plotdf_Mmat[!is.na(plotdf_Mmat$diag), ],
                             ggplot2::aes(Var3, Var2, fill = value, color = diag), linewidth = .7) +
          ggplot2::scale_color_manual(guide = "none", values = c(`TRUE` = "blue3")) +
          ggplot2::geom_text(ggplot2::aes(Var3, Var2, label = value.labels_Mmat,
                                          size = value^(1/2)),
                             color = "black",
                             fontface = 'bold') +
          # ggplot2::scale_size(range = c(0, 8)) +
          ggplot2::scale_fill_gradient(low="white", high="red3",
                                       breaks = seq(0, 100, 20), limits = c(0,100),
                                       name = 'Classification Percentage') +
          ggplot2::facet_grid(.~Var1) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            # axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text.x = ggplot2::element_text(color = "black",
                                                angle = 30, hjust = 1, vjust = 1),
            axis.text.y = ggplot2::element_text(color = "black"),
            # axis.ticks.x = ggplot2::element_blank(),
            # axis.ticks.length.x = unit(.2, "cm"),
            # axis.ticks.y = ggplot2::element_line(linewidth = .5),
            # axis.ticks.length.y = ggplot2::unit(.2, "cm"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 1),
            panel.grid.major = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            panel.grid.minor = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            strip.text.x = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.text.y = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.background = ggplot2::element_rect(color="black", linewidth=1),
            legend.title = ggplot2::element_blank(),
            # legend.key.width = ggplot2::unit(.75, "cm"),
            # legend.key.height = ggplot2::unit(.75, "cm"),
            # legend.key.size = ggplot2::unit(.5, "cm"),
            # legend.spacing.x = ggplot2::unit(.5, 'cm'),
            # legend.text=ggplot2::element_text(size=22),
            legend.text=ggplot2::element_text(size=12),
            legend.position = 'none'
          ) +
          ggplot2::labs(title = title_Mmat,
                        x = 'VA Cause', y = 'CHAMPS Cause')
        # print(ggplot2_Mmat)

      }


    }else if(calibmodel.type=="Mmatprior"){


      title_Mmat_input = title_Mmat = "Prior Mean of Misclassification Matrix"

      if(stable){

        title_Mmat_input = paste0(title_Mmat_input, " (Estimated from CHAMPS Data)")
        title_Mmat = paste0(title_Mmat, " (Adjusted)")

      }


      # Mmat input
      Mmat_input_toplot = output$Mmat.asDirich_input
      # Mmat_input_toplot = array(dim = dim(output$Mmat.asDirich_input),
      #                           dimnames = dimnames(output$Mmat.asDirich_input))
      for(k in 1:dim(Mmat_input_toplot)[1]){

        Mmat_input_toplot[k,,] = diag(nCause)
        Mmat_input_toplot[k,!donotcalib_asmat[k,],!donotcalib_asmat[k,]] =
          output$Mmat.asDirich_input[k,!donotcalib_asmat[k,],!donotcalib_asmat[k,]]/rowSums(output$Mmat.asDirich_input[k,!donotcalib_asmat[k,],!donotcalib_asmat[k,]])
        Mmat_input_toplot[k,,] = round(100*(Mmat_input_toplot[k,,]/rowSums(Mmat_input_toplot[k,,], na.rm = TRUE)))

        for(i in 1:dim(Mmat_input_toplot)[2]){

          if(!any(is.nan(Mmat_input_toplot[k,i,]))){

            nonzerocauseid_ki = which.max(Mmat_input_toplot[k,i,])
            Mmat_input_toplot[k,i,nonzerocauseid_ki] = 100 - sum(Mmat_input_toplot[k,i,-nonzerocauseid_ki])

          }

        }

        Mmat_input_toplot[k,donotcalib_asmat[k,],] =
          Mmat_input_toplot[k,,donotcalib_asmat[k,]] = NA

      }

      plotdf_Mmat_input = reshape2::melt(Mmat_input_toplot)
      head(plotdf_Mmat_input)
      value.labels_Mmat_input = plotdf_Mmat_input$value

      plotdf_Mmat_input$Var1 = factor(x = plotdf_Mmat_input$Var1, levels = dimnames(Mmat_input_toplot)[[1]])
      plotdf_Mmat_input$Var2 = factor(x = plotdf_Mmat_input$Var2, levels = rev(causes))
      plotdf_Mmat_input$Var3 = factor(x = plotdf_Mmat_input$Var3, levels = causes)

      plotdf_Mmat_input$diag = plotdf_Mmat_input$Var2==plotdf_Mmat_input$Var3
      plotdf_Mmat_input$diag[!plotdf_Mmat_input$diag] = NA
      # print(head(plotdf_Mmat_input))

      if(K>1){

        ggplot2_Mmat_input =
          ggplot2::ggplot(plotdf_Mmat_input, ggplot2::aes(Var3, Var2, fill = value)) +
          ggplot2::geom_tile(color="white", linewidth=.5) +
          ggplot2::geom_tile(data = plotdf_Mmat_input[!is.na(plotdf_Mmat_input$diag), ],
                             ggplot2::aes(Var3, Var2, fill = value, color = diag), linewidth = .7) +
          ggplot2::scale_color_manual(guide = "none", values = c(`TRUE` = "blue3")) +
          ggplot2::geom_text(ggplot2::aes(Var3, Var2, label = value.labels_Mmat_input,
                                          size = value^(1/2)),
                             color = "black",
                             fontface = 'bold') +
          # ggplot2::scale_size(range = c(0, 8)) +
          ggplot2::scale_fill_gradient(low="white", high="red3",
                                       breaks = seq(0, 100, 20), limits = c(0,100),
                                       name = 'Classification Percentage') +
          ggplot2::facet_grid(.~Var1) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            # axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text.x = ggplot2::element_text(color = "black",
                                                angle = 30, hjust = 1, vjust = 1),
            axis.text.y = ggplot2::element_text(color = "black"),
            # axis.ticks.x = ggplot2::element_blank(),
            # axis.ticks.length.x = unit(.2, "cm"),
            # axis.ticks.y = ggplot2::element_line(linewidth = .5),
            # axis.ticks.length.y = ggplot2::unit(.2, "cm"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 1),
            panel.grid.major = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            panel.grid.minor = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            strip.text.x = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.text.y = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.background = ggplot2::element_rect(color="black", linewidth=1),
            legend.title = ggplot2::element_blank(),
            # legend.key.width = ggplot2::unit(.75, "cm"),
            # legend.key.height = ggplot2::unit(.75, "cm"),
            # legend.key.size = ggplot2::unit(.5, "cm"),
            # legend.spacing.x = ggplot2::unit(.5, 'cm'),
            # legend.text=ggplot2::element_text(size=22),
            legend.text=ggplot2::element_text(size=12),
            legend.position = 'none'
          ) +
          ggplot2::labs(title = title_Mmat_input,
                        x = 'VA Cause', y = 'CHAMPS Cause')
        # print(ggplot2_Mmat_input)

      }else{

        ggplot2_Mmat_input =
          ggplot2::ggplot(plotdf_Mmat_input, ggplot2::aes(Var3, Var2, fill = value)) +
          ggplot2::geom_tile(color="white", linewidth=.5) +
          ggplot2::geom_tile(data = plotdf_Mmat_input[!is.na(plotdf_Mmat_input$diag), ],
                             ggplot2::aes(Var3, Var2, fill = value, color = diag), linewidth = .7) +
          ggplot2::scale_color_manual(guide = "none", values = c(`TRUE` = "blue3")) +
          ggplot2::geom_text(ggplot2::aes(Var3, Var2, label = value.labels_Mmat_input,
                                          size = value^(1/2)),
                             color = "black",
                             fontface = 'bold') +
          # ggplot2::scale_size(range = c(0, 8)) +
          ggplot2::scale_fill_gradient(low="white", high="red3",
                                       breaks = seq(0, 100, 20), limits = c(0,100),
                                       name = 'Classification Percentage') +
          ggplot2::facet_grid(.~Var1) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            # axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text.x = ggplot2::element_text(color = "black",
                                                angle = 30, hjust = 1, vjust = 1),
            axis.text.y = ggplot2::element_text(color = "black"),
            # axis.ticks.x = ggplot2::element_blank(),
            # axis.ticks.length.x = unit(.2, "cm"),
            # axis.ticks.y = ggplot2::element_line(linewidth = .5),
            # axis.ticks.length.y = ggplot2::unit(.2, "cm"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 1),
            panel.grid.major = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            panel.grid.minor = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            strip.text.x = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.text.y = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.background = ggplot2::element_rect(color="black", linewidth=1),
            legend.title = ggplot2::element_blank(),
            # legend.key.width = ggplot2::unit(.75, "cm"),
            # legend.key.height = ggplot2::unit(.75, "cm"),
            # legend.key.size = ggplot2::unit(.5, "cm"),
            # legend.spacing.x = ggplot2::unit(.5, 'cm'),
            # legend.text=ggplot2::element_text(size=22),
            legend.text=ggplot2::element_text(size=12),
            legend.position = 'none'
          ) +
          ggplot2::labs(title = title_Mmat_input,
                        x = 'VA Cause', y = 'CHAMPS Cause')
        # print(ggplot2_Mmat_input)

      }


      # Mmat used for calibration
      Mmat_toplot = output$Mmat.asDirich
      # Mmat_toplot = array(dim = dim(output$Mmat.asDirich),
      #                     dimnames = dimnames(output$Mmat.asDirich))
      for(k in 1:dim(Mmat_toplot)[1]){

        Mmat_toplot[k,,] = diag(nCause)
        Mmat_toplot[k,!donotcalib_asmat[k,],!donotcalib_asmat[k,]] =
          output$Mmat.asDirich[k,!donotcalib_asmat[k,],!donotcalib_asmat[k,]]/rowSums(output$Mmat.asDirich[k,!donotcalib_asmat[k,],!donotcalib_asmat[k,]])
        Mmat_toplot[k,,] = round(100*(Mmat_toplot[k,,]/rowSums(Mmat_toplot[k,,], na.rm = TRUE)))

        for(i in 1:dim(Mmat_toplot)[2]){

          if(!any(is.nan(Mmat_toplot[k,i,]))){

            nonzerocauseid_ki = which.max(Mmat_toplot[k,i,])
            Mmat_toplot[k,i,nonzerocauseid_ki] = 100 - sum(Mmat_toplot[k,i,-nonzerocauseid_ki])

          }

        }

        Mmat_toplot[k,donotcalib_asmat[k,],] =
          Mmat_toplot[k,,donotcalib_asmat[k,]] = NA

      }

      plotdf_Mmat = reshape2::melt(Mmat_toplot)
      head(plotdf_Mmat)
      value.labels_Mmat = plotdf_Mmat$value

      plotdf_Mmat$Var1 = factor(x = plotdf_Mmat$Var1, levels = dimnames(Mmat_toplot)[[1]])
      plotdf_Mmat$Var2 = factor(x = plotdf_Mmat$Var2, levels = rev(causes))
      plotdf_Mmat$Var3 = factor(x = plotdf_Mmat$Var3, levels = causes)

      plotdf_Mmat$diag = plotdf_Mmat$Var2==plotdf_Mmat$Var3
      plotdf_Mmat$diag[!plotdf_Mmat$diag] = NA
      head(plotdf_Mmat)

      if(K>1){

        ggplot2_Mmat =
          ggplot2::ggplot(plotdf_Mmat, ggplot2::aes(Var3, Var2, fill = value)) +
          ggplot2::geom_tile(color="white", linewidth=.5) +
          ggplot2::geom_tile(data = plotdf_Mmat[!is.na(plotdf_Mmat$diag), ],
                             ggplot2::aes(Var3, Var2, fill = value, color = diag), linewidth = .7) +
          ggplot2::scale_color_manual(guide = "none", values = c(`TRUE` = "blue3")) +
          ggplot2::geom_text(ggplot2::aes(Var3, Var2, label = value.labels_Mmat,
                                          size = value^(1/2)),
                             color = "black",
                             fontface = 'bold') +
          # ggplot2::scale_size(range = c(0, 8)) +
          ggplot2::scale_fill_gradient(low="white", high="red3",
                                       breaks = seq(0, 100, 20), limits = c(0,100),
                                       name = 'Classification Percentage') +
          ggplot2::facet_grid(.~Var1) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            # axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text.x = ggplot2::element_text(color = "black",
                                                angle = 30, hjust = 1, vjust = 1),
            axis.text.y = ggplot2::element_text(color = "black"),
            # axis.ticks.x = ggplot2::element_blank(),
            # axis.ticks.length.x = unit(.2, "cm"),
            # axis.ticks.y = ggplot2::element_line(linewidth = .5),
            # axis.ticks.length.y = ggplot2::unit(.2, "cm"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 1),
            panel.grid.major = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            panel.grid.minor = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            strip.text.x = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.text.y = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.background = ggplot2::element_rect(color="black", linewidth=1),
            legend.title = ggplot2::element_blank(),
            # legend.key.width = ggplot2::unit(.75, "cm"),
            # legend.key.height = ggplot2::unit(.75, "cm"),
            # legend.key.size = ggplot2::unit(.5, "cm"),
            # legend.spacing.x = ggplot2::unit(.5, 'cm'),
            # legend.text=ggplot2::element_text(size=22),
            legend.text=ggplot2::element_text(size=12),
            legend.position = 'none'
          ) +
          ggplot2::labs(title = title_Mmat,
                        x = 'VA Cause', y = 'CHAMPS Cause')
        # print(ggplot2_Mmat)

      }else{

        ggplot2_Mmat =
          ggplot2::ggplot(plotdf_Mmat, ggplot2::aes(Var3, Var2, fill = value)) +
          ggplot2::geom_tile(color="white", linewidth=.5) +
          ggplot2::geom_tile(data = plotdf_Mmat[!is.na(plotdf_Mmat$diag), ],
                             ggplot2::aes(Var3, Var2, fill = value, color = diag), linewidth = .7) +
          ggplot2::scale_color_manual(guide = "none", values = c(`TRUE` = "blue3")) +
          ggplot2::geom_text(ggplot2::aes(Var3, Var2, label = value.labels_Mmat,
                                          size = value^(1/2)),
                             color = "black",
                             fontface = 'bold') +
          # ggplot2::scale_size(range = c(0, 8)) +
          ggplot2::scale_fill_gradient(low="white", high="red3",
                                       breaks = seq(0, 100, 20), limits = c(0,100),
                                       name = 'Classification Percentage') +
          ggplot2::facet_grid(.~Var1) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            # axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text.x = ggplot2::element_text(color = "black",
                                                angle = 30, hjust = 1, vjust = 1),
            axis.text.y = ggplot2::element_text(color = "black"),
            # axis.ticks.x = ggplot2::element_blank(),
            # axis.ticks.length.x = unit(.2, "cm"),
            # axis.ticks.y = ggplot2::element_line(linewidth = .5),
            # axis.ticks.length.y = ggplot2::unit(.2, "cm"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 1),
            panel.grid.major = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            panel.grid.minor = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                                     colour = "grey90"),
            strip.text.x = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.text.y = ggplot2::element_text(
              size = 13,
              face = "bold"
            ),
            strip.background = ggplot2::element_rect(color="black", linewidth=1),
            legend.title = ggplot2::element_blank(),
            # legend.key.width = ggplot2::unit(.75, "cm"),
            # legend.key.height = ggplot2::unit(.75, "cm"),
            # legend.key.size = ggplot2::unit(.5, "cm"),
            # legend.spacing.x = ggplot2::unit(.5, 'cm'),
            # legend.text=ggplot2::element_text(size=22),
            legend.text=ggplot2::element_text(size=12),
            legend.position = 'none'
          ) +
          ggplot2::labs(title = title_Mmat,
                        x = 'VA Cause', y = 'CHAMPS Cause')
        # print(ggplot2_Mmat)

      }


    }


    # calibrated vs uncalibrated csmf
    plotdf_pcalib = NULL
    for(k in 1:dim(output$pcalib_postsumm)[1]){

      plotdf_pcalib = rbind.data.frame(plotdf_pcalib,
                                       cbind.data.frame(rbind.data.frame(data.frame('causes' = dimnames(output$pcalib_postsumm)[[3]],
                                                                                    'value' = unname(output$pcalib_postsumm[k,'postmean',]),
                                                                                    'llim' = unname(output$pcalib_postsumm[k,'lowcredI',]),
                                                                                    'ulim' = unname(output$pcalib_postsumm[k,'upcredI',]),
                                                                                    'calib_type' = 'Calibrated'),
                                                                         data.frame('causes' = colnames(output$p_uncalib),
                                                                                    'value' = unname(output$p_uncalib[k,]),
                                                                                    'llim' = NA,
                                                                                    'ulim' = NA,
                                                                                    'calib_type' = 'Uncalibrated')),
                                                        'vaalgo' = (rownames(output$p_uncalib))[k]))

    }
    head(plotdf_pcalib)

    plotdf_pcalib$causes = factor(x = plotdf_pcalib$causes,
                                  levels = colnames(output$p_uncalib))
    plotdf_pcalib$calib_type = factor(x = plotdf_pcalib$calib_type,
                                      levels = c("Uncalibrated", "Calibrated"))
    plotdf_pcalib$vaalgo = factor(x = plotdf_pcalib$vaalgo,
                                  levels = rownames(output$p_uncalib))
    head(plotdf_pcalib)

    if(K>1){

      ggplot2_pcalib =
        ggplot2::ggplot(data = plotdf_pcalib) +
        ggplot2::facet_grid(.~vaalgo) +
        ggplot2::coord_cartesian(ylim = c(0,1), expand = TRUE, default = FALSE, clip = 'on') +
        ggplot2::geom_col(ggplot2::aes(x = causes, y = value,
                                       fill = calib_type, color = calib_type),
                          linewidth = .3, alpha = .5, #color = "black",
                          position = ggplot2::position_dodge(width = .6), width = 0.5) +
        ggplot2::geom_errorbar(ggplot2::aes(x = causes, ymin = llim, ymax = ulim,
                                            color = calib_type),
                               # color = "black",
                               width = .3, linewidth = 1,
                               position = ggplot2::position_dodge(width = .6)) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.x = ggplot2::element_text(color = "black",
                                              angle = 30, hjust = 1, vjust = 1),
          axis.text.y = ggplot2::element_text(color = "black"),
          # axis.ticks.x = ggplot2::element_line(linewidth = .5),
          # axis.ticks.length.x = ggplot2::unit(.2, "cm"),
          # axis.ticks.y = ggplot2::element_line(linewidth = .5),
          # axis.ticks.length.y = ggplot2::unit(.2, "cm"),
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                               fill = NA, linewidth = 1),
          panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                   colour = "grey90"),
          panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                   colour = "grey90"),
          strip.text.x = ggplot2::element_text(
            size = 13,
            face = "bold"
          ),
          strip.text.y = ggplot2::element_text(
            size = 13,
            face = "bold"
          ),
          strip.background = ggplot2::element_rect(color="black", linewidth=1),
          legend.title = ggplot2::element_blank(),
          # legend.key.width = ggplot2::unit(1.5, "cm"),
          # legend.key.height = ggplot2::unit(.75, "cm"),
          # legend.key.spacing.x = ggplot2::unit(1, 'cm'),
          # legend.text=ggplot2::element_text(size=20),
          legend.text=ggplot2::element_text(size=12),
          legend.position = 'bottom'
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow=FALSE),
                        color = "none") +
        ggplot2::labs(title = ifelse(is.null(age_group),
                                     paste0('Cause-Specific Mortality Fractions (CSMF)'),
                                     paste0(age_group, ': Cause-Specific Mortality Fractions (CSMF)')),
                      x = "Cause",
                      y = 'Estimate')
      # ggplot2_pcalib

    }else{

      ggplot2_pcalib =
        ggplot2::ggplot(data = plotdf_pcalib) +
        ggplot2::facet_grid(.~vaalgo) +
        ggplot2::coord_cartesian(ylim = c(0,1), expand = TRUE, default = FALSE, clip = 'on') +
        ggplot2::geom_col(ggplot2::aes(x = causes, y = value,
                                       fill = calib_type, color = calib_type),
                          linewidth = .3, alpha = .5, #color = "black",
                          position = ggplot2::position_dodge(width = .6), width = 0.5) +
        ggplot2::geom_errorbar(ggplot2::aes(x = causes, ymin = llim, ymax = ulim,
                                            color = calib_type),
                               # color = "black",
                               width = .3, linewidth = 1,
                               position = ggplot2::position_dodge(width = .6)) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.x = ggplot2::element_text(color = "black",
                                              angle = 30, hjust = 1, vjust = 1),
          axis.text.y = ggplot2::element_text(color = "black"),
          # axis.ticks.x = ggplot2::element_line(linewidth = .5),
          # axis.ticks.length.x = ggplot2::unit(.2, "cm"),
          # axis.ticks.y = ggplot2::element_line(linewidth = .5),
          # axis.ticks.length.y = ggplot2::unit(.2, "cm"),
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                               fill = NA, linewidth = 1),
          panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                   colour = "grey90"),
          panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                   colour = "grey90"),
          strip.text.x = ggplot2::element_text(
            size = 13,
            face = "bold"
          ),
          strip.text.y = ggplot2::element_text(
            size = 13,
            face = "bold"
          ),
          strip.background = ggplot2::element_rect(color="black", linewidth=1),
          legend.title = ggplot2::element_blank(),
          # legend.key.width = ggplot2::unit(1.5, "cm"),
          # legend.key.height = ggplot2::unit(.75, "cm"),
          # legend.key.spacing.x = ggplot2::unit(1, 'cm'),
          # legend.text=ggplot2::element_text(size=20),
          legend.text=ggplot2::element_text(size=12),
          legend.position = 'bottom'
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow=FALSE),
                        color = "none") +
        ggplot2::labs(title = ifelse(is.null(age_group),
                                     paste0('Cause-Specific Mortality Fractions (CSMF)'),
                                     paste0(age_group, ': Cause-Specific Mortality Fractions (CSMF)')),
                      x = "Cause",
                      y = 'Estimate')
      # ggplot2_pcalib

    }


    if(K==1){

      if(stable){

        # print(ggplot2_Mmat_input)
        # print(ggplot2_Mmat)
        combinedplot = (ggplot2_Mmat_input | ggplot2_Mmat | ggplot2_pcalib)
        print(combinedplot)

      }else{

        combinedplot = (ggplot2_Mmat | ggplot2_pcalib)
        print(ggplot2_Mmat | ggplot2_pcalib)

      }

    }else if(K>1){

      if(stable){

        if(ensemble){

          print(((ggplot2_Mmat_input | NULL) + plot_layout(widths = c(K, 1))) / ((ggplot2_Mmat | NULL) + plot_layout(widths = c(K, 1))) / ggplot2_pcalib)

        }else{

          print(ggplot2_Mmat_input / ggplot2_Mmat / ggplot2_pcalib)

        }

      }else{

        if(ensemble){

          print(((ggplot2_Mmat | NULL) + plot_layout(widths = c(K, 1))) / ggplot2_pcalib)

        }else{

          print(ggplot2_Mmat / ggplot2_pcalib)

        }

      }

      # print(((ggplot2_Mmat | NULL) + plot_layout(widths = c(K, 1))) / (ggplot2_pcalib))

    }

  }


  if(saveoutput){

    if(is.null(output_filename)) output_filename = paste0("calibratedva_", calibmodel.type)

    saveRDS(output, file.path(tempdir(), output_filename))

  }else{

    return(output)

  }

}


