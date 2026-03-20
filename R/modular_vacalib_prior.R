#' Modular VA-Calibration using Dirichlet Prior on Misclassification Matrix
#'
#' This is a utility function. Please use \link{vacalibration}.
#'
#' @param va_unlabeled Same as \code{va_unlabeled} in \code{vacalibration()}
#' @param Mmat_calib Same as \code{missmat} in \code{vacalibration()}
#' @param studycause_map Same as \code{studycause_map} in \code{vacalibration()}
#' @param donotcalib Same as \code{donotcalib} in \code{vacalibration()}
#' @param donotcalib_type Same as \code{donotcalib_type} in \code{vacalibration()}
#' @param nocalib.threshold Same as \code{nocalib.threshold} in \code{vacalibration()}
#' @param path_correction Same as \code{path_correction} in \code{vacalibration()}
#' @param ensemble Same as \code{ensemble} in \code{vacalibration()}
#' @param pshrink_strength Same as \code{pshrink_strength} in \code{vacalibration()}
#' @param nMCMC,nBurn,nThin Same as \code{nMCMC}, \code{nBurn}, and \code{nThin} in \code{vacalibration()}
#' @param nChain Same as \code{nChain} in \code{vacalibration()}
#' @param nCore Same as \code{nCore} in \code{vacalibration()}
#' @param adapt_delta_stan Same as \code{adapt_delta_stan} in \code{vacalibration()}
#' @param refresh_stan Same as \code{refresh_stan} in \code{vacalibration()}
#' @param seed Same as \code{seed} in \code{vacalibration()}
#' @param verbose Same as \code{verbose} in \code{vacalibration()}
#' @param input_vacalib List of inputs in \code{vacalibration()}
#' @return Similar to the list returned in \code{vacalibration()}
#'
#' @aliases modular_vacalib_prior
#'
#' @import rstan
#'
#' @importFrom stats quantile
#' @importFrom utils head
#' @importFrom MASS ginv
#'
modular_vacalib_prior <- function(va_unlabeled,
                                  Mmat_calib, studycause_map,
                                  donotcalib, donotcalib_type, nocalib.threshold,
                                  path_correction,
                                  ensemble,
                                  pshrink_strength,
                                  nMCMC, nBurn, nThin,
                                  nChain, nCore,
                                  adapt_delta_stan, refresh_stan,
                                  seed, verbose, input_vacalib){

  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = nCore)


  # va_data preparation for STAN implementation ----
  K = length(va_unlabeled) # number of algorithms

  for(k in 1:K){

    if(k==1){

      causes = names(va_unlabeled[[k]])
      nCause = length(causes)

      va_deaths = array(dim = c(K, nCause),
                        dimnames = list(names(va_unlabeled), causes))

    }

    va_deaths[k,] = va_unlabeled[[k]]

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


  # preparing misclassification matrix ----

  if(verbose){

    # cat("\n")
    message("* Preparing misclassification for calibration")

  }

  Mmat_calib_input = Mmat_calib_study = Mmat_calib
  # print(Mmat_calib_input)

  eps = 0.001

  if(!is.null(studycause_map)){

    ## study causes outside CHAMPS ----

    # if(verbose){
    #
    #   # cat("\n")
    #   message("Mapping study causes to CHAMPS causes ... ")
    #
    # }

    map.cause = data.frame('study' = names(va_unlabeled[[1]]),
                           'champs' = unname(studycause_map[names(va_unlabeled[[1]])]))

    # champs cause frequencies to distribute false positives or maintain sensitivity
    champs_cause.freq = table(map.cause$champs)
    multiplied_champs_cause = names(champs_cause.freq)[champs_cause.freq>1]
    unique_champscause = unique(map.cause$champs)

    Mmat_calib_study = lapply(1:length(va_unlabeled),
                              FUN = function(k){

                                asDirich_champs_k = Mmat_calib_input[[k]][unique_champscause,unique_champscause]
                                rownames(asDirich_champs_k) = colnames(asDirich_champs_k) = unique_champscause
                                shapesumDirich_champs_k = rowSums(asDirich_champs_k)
                                mean_champs_k = asDirich_champs_k/shapesumDirich_champs_k

                                ## cause matching
                                Mmat_mean_k = mean_champs_k[map.cause$champs,map.cause$champs]
                                rownames(Mmat_mean_k) = colnames(Mmat_mean_k) = map.cause$study

                                ## expanding classification rates to study-specific causes
                                if(length(multiplied_champs_cause)>0){

                                  for(l in 1:length(multiplied_champs_cause)){

                                    (id_l = which(map.cause$champs %in% multiplied_champs_cause[l]))
                                    (idothers_l = which(!(map.cause$champs %in% multiplied_champs_cause[l])))

                                    tempmat = matrix(data = eps/(champs_cause.freq[multiplied_champs_cause[l]] - 1),
                                                     nrow = length(id_l), ncol = length(id_l))
                                    diag(tempmat) = 1 - eps
                                    Mmat_mean_k[id_l, id_l] = tempmat*Mmat_mean_k[id_l, id_l]

                                    Mmat_mean_k[idothers_l, id_l] =
                                      Mmat_mean_k[idothers_l, id_l]/champs_cause.freq[multiplied_champs_cause[l]]

                                  }

                                }

                                ## misclassification prior to use
                                Mmat_calib_study_k = shapesumDirich_champs_k[map.cause$champs]*Mmat_mean_k
                                rownames(Mmat_calib_study_k) = colnames(Mmat_calib_study_k) = map.cause$study

                                Mmat_calib_study_k

                              })
    names(Mmat_calib_study) = names(va_unlabeled)

    # if(verbose){
    #
    #   # cat("\n")
    #   message("Mapping study causes to CHAMPS causes ... Done.")
    #   message("")
    #
    # }

  }
  # print(Mmat_calib_study)

  if(verbose){

    # cat("\n")
    # message("Preparing misclassification for calibration ... Done.")
    message("")

  }


  # causes not calibrated ----
  # donotcalib_input = donotcalib
  # if(is.null(studycause_map)){
  #
  #   donotcalib_study = donotcalib_input
  #
  # }else{
  #
  #   donotcalib_study = lapply(1:length(va_unlabeled),
  #                             FUN = function(k){
  #
  #                               if(is.null(donotcalib[[k]])){
  #
  #                                 return(donotcalib[[k]])
  #
  #                               }else if((!all(donotcalib[[k]] %in% names(va_data_tomodel[[k]])))&(!all(donotcalib[[k]] %in% colnames(va_data_tomodel[[k]])))){
  #
  #                                 message(paste0(donotcalib[[k]][!(donotcalib[[k]] %in% names(va_data_tomodel[[k]]))], collapse = ', '))
  #                                 message(paste0("The above causes provided in 'donotcalib' for ", names(va_data_tomodel)[k], " algorithm are not present in its 'va_data'."))
  #                                 stop("")
  #
  #                               }
  #
  #                               names(studycause_map)[studycause_map %in% donotcalib_input[[k]]]
  #
  #                             })
  #   names(donotcalib_study) = names(va_unlabeled)
  #
  # }
  donotcalib_study = donotcalib


  # calibration for each algorithm ----

  # storage
  donotcalib_study_asmat = donotcalib_tomodel_asmat = matrix(nrow = K, ncol = nCause)

  calibout = vector(mode = "list", length = K)

  calibrated = lambda_calibpath = rep(NA, K)

  rownames(donotcalib_study_asmat) = rownames(donotcalib_tomodel_asmat) =
    names(calibout) = names(Mmat_calib_study)
  colnames(donotcalib_study_asmat) = colnames(donotcalib_tomodel_asmat) = causes

  Mmat_calib_input_asarray =
    array(dim = c(K, dim(Mmat_calib_input[[1]])),
          dimnames = list(names(Mmat_calib_input),
                          rownames(Mmat_calib_input[[1]]),
                          colnames(Mmat_calib_input[[1]])))

  Mmat_calib_study_asarray =
    Mmat_calib_tomodel_asarray =
    array(dim = c(K, nCause, nCause),
          dimnames = list(names(Mmat_calib_study), causes, causes))

  va_deaths_calib = va_deaths

  pcalib_asarray = array(dim = c(K, nMCMC, nCause),
                         dimnames = list(names(Mmat_calib_study), NULL, causes))
  pcalib_postsumm_asarray = array(dim = c(K, 3, nCause),
                                  dimnames = list(names(Mmat_calib_study), c('postmean', 'lowcredI', 'upcredI'), causes))

  for(k in 1:K){

    donotcalib_study_asmat[k,] =
      donotcalib_tomodel_asmat[k,] =
      causes %in% donotcalib_study[[k]]

    # learn more causes not to calibrate
    if(donotcalib_type=="learn"){

      # print(Mmat_calib_study[[k]])
      donotcalib_study_learn_k = apply(Mmat_calib_study[[k]]/rowSums(Mmat_calib_study[[k]]), 2,
                                       FUN = function(v){

                                         # mean(abs(v - mean(v)))
                                         diff(range(v))<=nocalib.threshold

                                       })
      donotcalib_tomodel_asmat[k,] = (donotcalib_study_asmat[k,]|donotcalib_study_learn_k)

    }


    ## calibration path correction ----
    Mmat_calib_input_asarray[k,,] = Mmat_calib_input[[k]]
    Mmat_calib_study_asarray[k,,] = Mmat_calib_study[[k]]

    if(sum(!donotcalib_tomodel_asmat[k,])<=1){

      ### no calibration ----

      if(verbose){

        message(paste0("** Cannot calibrate ", names(va_unlabeled)[k]))
        message("")

      }

      calibrated[k] = FALSE

      # Mmat_calib_tomodel_asarray[k,,] = Mmat_calib_study_asarray[k,,]

      #### calibration output ----
      MCMCout_k = list('p_calib' = matrix(data = puncalib[k,],
                                          nrow = 1, ncol = nCause,
                                          byrow = TRUE),
                       # 'loglik' = NULL,
                       # 'mcmc.diagnostic' = NULL,
                       'p_calib_postsumm' = rbind(puncalib[k,],
                                                  puncalib[k,],
                                                  puncalib[k,]))

      # posterior summary of calibrated estimate
      colnames(MCMCout_k$p_calib) = causes
      rownames(MCMCout_k$p_calib_postsumm) = c('postmean', 'lowcredI', 'upcredI')

    }else{

      ### calibration ----

      if(verbose){

        message(paste0("* Calibrating ", names(va_unlabeled)[k]))

      }

      calibrated[k] = TRUE

      idtocalib_k = which(!donotcalib_tomodel_asmat[k,])
      nTocalib_k = length(idtocalib_k)


      Mmat_calib_tomodel_asarray[k,,] = Mmat_calib_study_asarray[k,,]
      puncalib_k_sub = puncalib[k,!donotcalib_tomodel_asmat[k,]]/sum(puncalib[k,!donotcalib_tomodel_asmat[k,]])

      if(!path_correction){

        #### no path correction ----

        lambda_calibpath[k] = 0

      }else{

        #### path correction ----
        Mmat_calib_tomodel_sub_k = Mmat_calib_study_asarray[k,!donotcalib_tomodel_asmat[k,],!donotcalib_tomodel_asmat[k,]]
        nSamp_approx = 1000
        Imat = diag(nrow(Mmat_calib_tomodel_sub_k))
        Imat_asarray = Mmat_calib_tomodel_sub_samples_k =
          array(dim = c(nSamp_approx, nrow(Mmat_calib_tomodel_sub_k),
                        nrow(Mmat_calib_tomodel_sub_k)))
        for(i in 1:nrow(Imat)){

          Imat_asarray[,i,] = matrix(data = Imat[i,], nrow = nSamp_approx,
                                     ncol = ncol(Imat),
                                     byrow = T)
          Mmat_calib_tomodel_sub_samples_k[,i,] = LaplacesDemon::rdirichlet(n = nSamp_approx,
                                                                            alpha = Mmat_calib_tomodel_sub_k[i,])

        }
        # print(round(apply(Mmat_calib_tomodel_sub_samples_k, 2:3, mean), 2))
        # print(round(apply(Mmat_calib_tomodel_sub_samples_k, 2:3, var), 2))

        count_lambda_calibpath_k = 0
        got_it_k = FALSE
        while(!got_it_k){

          count_lambda_calibpath_k = count_lambda_calibpath_k + 1
          if(count_lambda_calibpath_k==1){

            lambda_calibpath[k] = .99

          }else{

            lambda_calibpath[k] = lambda_calibpath[k] - .01

          }

          # shrink towards identity
          Mmat_calib_tomodel_sub_samples_k_lambda =
            lambda_calibpath[k]*Imat_asarray +
            (1-lambda_calibpath[k])*Mmat_calib_tomodel_sub_samples_k
          Mmat_calib_tomodel_sub_samples_avg_k_lambda = apply(Mmat_calib_tomodel_sub_samples_k_lambda,
                                                              2:3, mean)

          # solving calibration eq
          puncalib_k_sub_lambda = as.numeric(MASS::ginv(t(Mmat_calib_tomodel_sub_samples_avg_k_lambda)) %*%
                                               as.matrix(puncalib_k_sub))

          # pcalib_lambda = puncalib
          # pcalib_lambda[!donotcalib_tomodel_asmat[k,]] = (1-sum(puncalib[donotcalib_tomodel_asmat[k,]]))*pcalib_lambda

          # checking if within simplex
          outside_simplex_k = any(puncalib_k_sub_lambda<0)|any(puncalib_k_sub_lambda>1)
          got_it_k = outside_simplex_k|isTRUE(all.equal(lambda_calibpath[k], 0))

        }

        if(outside_simplex_k){

          lambda_calibpath[k] = min(lambda_calibpath[k] + .01, .99)

        }
        # print(1-lambda_calibpath[k])

        # shrink towards identity
        Mmat_calib_tomodel_sub_samples_k_lambda =
          lambda_calibpath[k]*Imat_asarray +
          (1-lambda_calibpath[k])*Mmat_calib_tomodel_sub_samples_k
        Mmat_calib_tomodel_sub_samples_avg_k_lambda = apply(Mmat_calib_tomodel_sub_samples_k_lambda,
                                                            2:3, mean)

        # print(round(Mmat_calib_tomodel_sub_samples_avg_k_lambda,2))
        # print(round(apply(Mmat_calib_tomodel_sub_samples_k_lambda,
        #                   2:3, var), 2))

        Mmat_calib_tomodel_asarray[k,!donotcalib_tomodel_asmat[k,],!donotcalib_tomodel_asmat[k,]] =
          do.call("rbind",
                  lapply(1:dim(Mmat_calib_tomodel_sub_samples_k_lambda)[2],
                         FUN = function(i){

                           # print(Mmat_calib_tomodel_sub_samples_avg_k_lambda[i,])
                           # print(apply(Mmat_calib_tomodel_sub_samples_k_lambda[,i,],
                           #             2, mean))
                           # print(apply(Mmat_calib_tomodel_sub_samples_k_lambda[,i,],
                           #             2, var))
                           # print(apply(Mmat_calib_tomodel_sub_samples_k_lambda[,i,],
                           #             2, range))

                           mmat_nudge = (Mmat_calib_tomodel_sub_samples_k_lambda[,i,] + 1e-5)/rowSums(Mmat_calib_tomodel_sub_samples_k_lambda[,i,] + 1e-5)
                           # print(apply(mmat_nudge, 2, range))

                           mle.out = nlminb(start = 1,
                                            lower = 0,
                                            upper = Inf,
                                            objective = function(lambda_mmat){

                                              # -sum(LaplacesDemon::ddirichlet(x = Mmat_calib_tomodel_sub_samples_k_lambda[,i,],
                                              #                                alpha = dim(Mmat_calib_tomodel_sub_samples_k_lambda)[2]*lambda_mmat*Mmat_calib_tomodel_sub_samples_avg_k_lambda[i,],
                                              #                                log = T))
                                              -sum(LaplacesDemon::ddirichlet(x = mmat_nudge,
                                                                             alpha = dim(Mmat_calib_tomodel_sub_samples_k_lambda)[2]*lambda_mmat*Mmat_calib_tomodel_sub_samples_avg_k_lambda[i,],
                                                                             log = T))

                                            })
                           # print(paste0(i, " ", mle.out$par))

                           dim(Mmat_calib_tomodel_sub_samples_k_lambda)[2]*mle.out$par*Mmat_calib_tomodel_sub_samples_avg_k_lambda[i,]

                         }))

      }

      ### calibration output ----

      # if(shrink_towards=="calib"){
      #
      #   p0 = project_simplex(MASS::ginv(t(Mmat_calib_study_asarray[k,,]/rowSums(Mmat_calib_study_asarray[k,,]))) %*% puncalib[k,])
      #
      # }else if(shrink_towards=="uncalib"){
      #
      #   p0 = puncalib[k,]
      #
      # }
      p0 = puncalib[k,]

      if((sum(donotcalib_tomodel_asmat[k,])>0)&verbose){

        message(paste0("** Not calibrating: ", paste0(causes[donotcalib_tomodel_asmat[k,]], collapse = ', ')))

      }

      ### stan fit ----
      # print(list('nCause' = nCause,
      #            'nAlgo' = 1,
      #            'aj' = va_deaths[k,, drop = FALSE],
      #            'nTocalib' = as.array(nTocalib_k),
      #            'idtocalib' = idtocalib_k,
      #            'nCumtocalib' = c(0, cumsum(nTocalib_k)),
      #            'p0' = as.numeric(p0),
      #            'Mmatprior_asDirich' = Mmat_calib_tomodel_asarray[k,,, drop = FALSE],
      #            'pss' = pshrink_strength,
      #            'lambda' = lambda_eBayes,
      #            'Imat' = diag(nCause)
      # ))
      # print(lambda_calibpath)
      stanfit_k = rstan::sampling(get_stan_seqcalib(),
                                  # pars = c('p_calib', 'loglik'),
                                  # include = TRUE,
                                  data = list('nCause' = nCause,
                                              'nAlgo' = 1,
                                              'aj' = va_deaths[k,, drop = FALSE],
                                              'nTocalib' = as.array(nTocalib_k),
                                              'idtocalib' = idtocalib_k,
                                              'nCumtocalib' = c(0, cumsum(nTocalib_k)),
                                              'p0' = as.numeric(p0),
                                              'Mmatprior_asDirich' = Mmat_calib_tomodel_asarray[k,,, drop = FALSE],
                                              'pss' = pshrink_strength,
                                              'lambda' = 0,
                                              'Imat' = diag(nCause)
                                  ),
                                  chains = nChain, cores = nCore,
                                  iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                                  control = list('adapt_delta' = adapt_delta_stan),
                                  seed = seed,
                                  refresh = refresh_stan#,
                                  # init = list(list('Mmat' = Mmat.init,
                                  #                  'p_calib' = puncalib))
      )

      MCMCout_k = c(list("stanfit" = stanfit_k), rstan::extract(stanfit_k)) # MCMC output
      # MCMCout_k$calibrated = calibrated_k
      # MCMCout_k$lambda = lambda_k

      ### mcmc diagnostic ----
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

      # mcmc diagnostic summary
      MCMCout_k$mcmc.diagnostic = c('max_Rhat' = max_Rhat_k, 'min_ess_bulk' = min_ess_bulk_k,
                                    'num_divergent' = rstan::get_num_divergent(stanfit_k),
                                    'num_max_treedepth' = rstan::get_num_max_treedepth(stanfit_k))

      # if(verbose){
      #
      #   message("MCMC Diagnostic:")
      #   print(round(MCMCout_k$mcmc.diagnostic, 3))
      #   message("\n")
      #
      # }

      ### calibrated csmf for all causes ----
      # print(head(MCMCout$p_calib))
      if(sum(donotcalib_tomodel_asmat[k,])>0){

        # csmf
        p_calib_out_k = matrix(data = puncalib[k,],
                               nrow = nMCMC, ncol = nCause,
                               byrow = TRUE)
        # print(dim(p_calib_out[,donotcalib_touse, drop = F]))
        # p_calib_out[,donotcalib_touse] = matrix(data = puncalib[donotcalib_touse],
        #                                         nrow = nMCMC, ncol = sum(donotcalib_touse),
        #                                         byrow = T)
        p_calib_out_k[,!donotcalib_tomodel_asmat[k,]] =
          (1 - sum(puncalib[k,donotcalib_tomodel_asmat[k,]]))*
          (MCMCout_k$p_calib[,!donotcalib_tomodel_asmat[k,]]/rowSums(MCMCout_k$p_calib[,!donotcalib_tomodel_asmat[k,]]))
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

      if(verbose){

        # message(paste0("Calibrating ", names(va_unlabeled)[k], " ... Done."))
        message("")

      }

    }



    ## MCMC output ----
    calibout[[k]] = MCMCout_k


    ## calibrated number of deaths ----
    va_deaths_calib_k = va_deaths[k,]
    if(calibrated[k]){

      va_deaths_calib_k[!donotcalib_tomodel_asmat[k,]] =
        smart_round(x = colMeans(sum(va_deaths[k,!donotcalib_tomodel_asmat[k,]])*
                                   (MCMCout_k$p_calib[,!donotcalib_tomodel_asmat[k,]]/(1 - sum(puncalib[k,donotcalib_tomodel_asmat[k,]])))),
                    target_sum = sum(va_deaths[k,!donotcalib_tomodel_asmat[k,]]), digits = 0)
      # print(va_deaths_calib_k)
      # print(sum(va_deaths_calib_k))
      # print(sum(va_deaths[k,]))

    }

    va_deaths_calib[k,] = va_deaths_calib_k
    pcalib_asarray[k,,] = MCMCout_k$p_calib
    pcalib_postsumm_asarray[k,,] = MCMCout_k$p_calib_postsumm

  }# ending calibration for each algorithm


  # ensemble calibration  ----
  if(K==1){

    if(is.null(ensemble)){

      ensemble = FALSE

    }else{

      if(ensemble){

        ensemble = FALSE
        if(verbose){

          message("** Nothing to ensemble. 'va_data' is provided for one algorithm.")
          message("")

        }

      }

    }

  }else if(K>1){

    # defaulting to ensemble if multiple algorithms
    if(is.null(ensemble)) ensemble = TRUE

  }


  # performing ensemble calibration
  if(!ensemble){

    # uncalibrated csmf as %
    puncalib_percentage = do.call("rbind",
                                         lapply(1:nrow(puncalib),
                                                FUN = function(k){

                                                  smart_round(x = 100*puncalib[k,],
                                                              target_sum = 100, digits = 0)

                                                }))
    rownames(puncalib_percentage) = rownames(puncalib)
    colnames(puncalib_percentage) = colnames(puncalib)
    # print(dim(puncalib_percentage))

    # point estimate of calibrated csmf as %
    pcalib_postmean_percentage = do.call("rbind",
                                         lapply(1:dim(pcalib_postsumm_asarray)[1],
                                                FUN = function(k){

                                                  smart_round(x = 100*pcalib_postsumm_asarray[k,1,],
                                                              target_sum = 100, digits = 0)

                                                }))
    rownames(pcalib_postmean_percentage) = dimnames(pcalib_postsumm_asarray)[[1]]
    colnames(pcalib_postmean_percentage) = dimnames(pcalib_postsumm_asarray)[[3]]
    # print(dim(pcalib_postmean_percentage))

    # output list
    # print(names(input_vacalib))
    input_vacalib.now = as.list(environment())
    # print(names(input_vacalib.now))
    # print(names(input_vacalib.now)[names(input_vacalib) %in% names(input_vacalib.now)])
    input_vacalib[intersect(names(input_vacalib), names(input_vacalib.now))] = input_vacalib.now[intersect(names(input_vacalib), names(input_vacalib.now))]
    # print(names(input_vacalib))

    output = list("calib_MCMCout" = calibout,
                  "p_uncalib" = puncalib,
                  "p_calib" = pcalib_asarray,
                  "pcalib_postsumm" = pcalib_postsumm_asarray,
                  "puncalib_percentage" = puncalib_percentage,
                  "pcalib_postmean_percentage" = pcalib_postmean_percentage,
                  "va_deaths_uncalib" = va_deaths,
                  "va_deaths_calib_algo" = va_deaths_calib,
                  'Mmat_input' = Mmat_calib_input_asarray,
                  'Mmat_study' = Mmat_calib_study_asarray,
                  'Mmat_tomodel' = Mmat_calib_tomodel_asarray,
                  'donotcalib_study' = donotcalib_study_asmat,
                  'donotcalib_tomodel' = donotcalib_tomodel_asmat,
                  'calibrated' = calibrated,
                  'lambda_calibpath' = lambda_calibpath,
                  "K" = K, "nCause" = nCause, "causes" = causes,
                  'input' = input_vacalib)

  }else{

    if(verbose){

      message("* Ensemble calibration")

    }

    # ensemble uncalibrated csmf
    puncalib = rbind(puncalib, colSums(va_deaths)/sum(va_deaths))

    donotcalib_ens = donotcalib_tomodel_asmat[1,]
    for(k in 2:K){

      donotcalib_ens = donotcalib_ens&donotcalib_tomodel_asmat[k,]

    }
    donotcalib_tomodel_asmat = rbind(donotcalib_tomodel_asmat,
                                     donotcalib_ens)

    rownames(puncalib) =
      rownames(donotcalib_tomodel_asmat) = c(names(va_unlabeled), 'ensemble')

    idtocalib_ens = nTocalib_ens = NULL
    for(k in 1:K){

      idtocalib_ens_k = which(!donotcalib_tomodel_asmat[k,])
      idtocalib_ens = c(idtocalib_ens, idtocalib_ens_k)
      nTocalib_ens = c(nTocalib_ens, length(idtocalib_ens_k))

    }

    calibrated = c(calibrated,
                   'ensemble' = (sum(!donotcalib_tomodel_asmat[K+1,])>1))
    if(!calibrated[K+1]){

      MCMCout_ens = list('p_calib' = matrix(data = puncalib[K+1,],
                                            nrow = nMCMC, ncol = nCause,
                                            byrow = TRUE),
                         # 'loglik' = NULL,
                         # "calibrated" = calibrated_ens,
                         # 'mcmc.diagnostic' = NULL,
                         'p_calib_postsumm' = rbind(puncalib[K+1,],
                                                    puncalib[K+1,],
                                                    puncalib[K+1,]))

      # posterior summary of calibrated estimate
      colnames(MCMCout_ens$p_calib) = causes
      rownames(MCMCout_ens$p_calib_postsumm) = c('postmean', 'lowcredI', 'upcredI')

    }else{

      # if(shrink_towards=="calib"){
      #
      #   M_neq = matrix(0, K, K)
      #   b_neq = numeric(K)
      #   for(k in 1:K){
      #
      #     Mtemp = diag(nCause)
      #     Mtemp[!donotcalib_tomodel_asmat[K+1,],!donotcalib_tomodel_asmat[K+1,]] =
      #       Mmat_calib_tomodel_asarray[k,!donotcalib_tomodel_asmat[K+1,],!donotcalib_tomodel_asmat[K+1,]]/rowSums(Mmat_calib_tomodel_asarray[k,!donotcalib_tomodel_asmat[K+1,],!donotcalib_tomodel_asmat[K+1,]])
      #
      #     M_neq = M_neq + tcrossprod(Mtemp)
      #     b_neq = b_neq + (Mtemp%*%puncalib[k,])
      #
      #   }
      #
      #   p0 = project_simplex(MASS::ginv(M_neq) %*% b_neq)
      #
      # }else if(shrink_towards=="uncalib"){
      #
      #   p0 = puncalib[K+1,]
      #
      # }
      p0 = puncalib[K+1,]

      if((sum(donotcalib_tomodel_asmat[K+1,])>0)&verbose){

        message(paste0("** Not calibrating: ", paste0(causes[donotcalib_tomodel_asmat[K+1,]], collapse = ', ')))

      }

      ## stan fit ----
      stanfit_ens = rstan::sampling(get_stan_seqcalib(),
                                    # pars = c('p_calib', 'loglik'),
                                    # include = TRUE,
                                    data = list('nCause' = nCause,
                                                'nAlgo' = K,
                                                'aj' = va_deaths,
                                                'nTocalib' = nTocalib_ens,
                                                'idtocalib' = idtocalib_ens,
                                                'nCumtocalib' = c(0, cumsum(nTocalib_ens)),
                                                'p0' = as.numeric(p0),
                                                'Mmatprior_asDirich' = Mmat_calib_tomodel_asarray,
                                                'pss' = pshrink_strength,
                                                'lambda' = 0,
                                                'Imat' = diag(nCause)
                                    ),
                                    chains = nChain, cores = nCore,
                                    iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                                    control = list('adapt_delta' = adapt_delta_stan),
                                    seed = seed,
                                    refresh = refresh_stan#,
                                    # init = list(list('Mmat' = Mmat.init,
                                    #                  'p_calib' = puncalib))
      )

      MCMCout_ens = c(list("stanfit" = stanfit_ens),
                      rstan::extract(stanfit_ens)) # MCMC output
      # MCMCout_ens$calibrated = calibrated_ens

      ## mcmc diagnostic ----
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

      # mcmc diagnostic summary
      MCMCout_ens$mcmc.diagnostic = c('max_Rhat' = max_Rhat_ens, 'min_ess_bulk' = min_ess_bulk_ens,
                                      'num_divergent' = rstan::get_num_divergent(stanfit_ens),
                                      'num_max_treedepth' = rstan::get_num_max_treedepth(stanfit_ens))

      # if(verbose){
      #
      #   message("MCMC Diagnostic:")
      #   print(round(MCMCout_k$mcmc.diagnostic, 2))
      #   # message("\n")
      #
      # }

      ## calibrated posterior for all causes ----
      # print(head(MCMCout$p_calib))
      if(sum(donotcalib_tomodel_asmat[K+1,])>0){

        # csmf
        p_calib_out_ens = matrix(data = puncalib[K+1,],
                                 nrow = nMCMC, ncol = nCause,
                                 byrow = TRUE)
        # print(dim(p_calib_out[,donotcalib_touse, drop = F]))
        # p_calib_out[,donotcalib_touse] = matrix(data = puncalib[donotcalib_touse],
        #                                         nrow = nMCMC, ncol = sum(donotcalib_touse),
        #                                         byrow = T)
        p_calib_out_ens[,!donotcalib_tomodel_asmat[K+1,]] = (1 - sum(puncalib[K+1,donotcalib_tomodel_asmat[K+1,]]))*
          (MCMCout_ens$p_calib[,!donotcalib_tomodel_asmat[K+1,]]/rowSums(MCMCout_ens$p_calib[,!donotcalib_tomodel_asmat[K+1,]]))
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

    ## calibrated number of deaths ----
    va_deaths_ens = va_deaths
    for(k in 1:K){

      va_deaths_ens_k = va_deaths[k,]

      if(calibrated[K+1]){

        # print(donotcalib_tomodel_asmat[K+1,])
        # print(va_deaths[k,!donotcalib_tomodel_asmat[K+1,]])
        va_deaths_ens_k[!donotcalib_tomodel_asmat[K+1,]] =
          smart_round(x = colMeans(sum(va_deaths[k,!donotcalib_tomodel_asmat[K+1,]])*
                                     (MCMCout_ens$p_calib[,!donotcalib_tomodel_asmat[K+1,]]/(1 - sum(puncalib[K+1,donotcalib_tomodel_asmat[K+1,]])))),
                      target_sum = sum(va_deaths[k,!donotcalib_tomodel_asmat[K+1,]]), digits = 0)

      }

      va_deaths_ens[k,] = va_deaths_ens_k
      # print(k)

    }


    # ensemble calibration output
    calibout = c(calibout, list('ensemble' = MCMCout_ens))

    # puncalib = rbind(puncalib, puncalib[K+1,])
    # rownames(puncalib) = names(calibout)
    # colnames(puncalib) = causes


    pcalib_asarray_wens = array(dim = c(K+1, nMCMC, nCause),
                                dimnames = list(rownames(puncalib), NULL, causes))
    pcalib_asarray_wens[1:K,,] = pcalib_asarray
    pcalib_asarray_wens[K+1,,] = MCMCout_ens$p_calib

    pcalib_postsumm_asarray_wens = array(dim = c(K+1, 3, nCause),
                                         dimnames = list(rownames(puncalib),
                                                         dimnames(pcalib_postsumm_asarray)[[2]],
                                                         causes))
    pcalib_postsumm_asarray_wens[1:K,,] = pcalib_postsumm_asarray
    pcalib_postsumm_asarray_wens[K+1,,] = MCMCout_ens$p_calib_postsumm

    # donotcalib_wens = rbind(donotcalib_tomodel_asmat, donotcalib_tomodel_asmat[K+1,])
    # rownames(donotcalib_wens) = names(calibout)
    # colnames(donotcalib_wens) = causes

    # uncalibrated csmf as %
    puncalib_percentage = do.call("rbind",
                                  lapply(1:nrow(puncalib),
                                         FUN = function(k){

                                           smart_round(x = 100*puncalib[k,],
                                                       target_sum = 100, digits = 0)

                                         }))
    rownames(puncalib_percentage) = rownames(puncalib)
    colnames(puncalib_percentage) = colnames(puncalib)
    # print(dim(puncalib_percentage))

    # point estimate of calibrated csmf as %
    pcalib_postmean_percentage = do.call("rbind",
                                         lapply(1:dim(pcalib_postsumm_asarray_wens)[1],
                                                FUN = function(k){

                                                  smart_round(x = 100*pcalib_postsumm_asarray_wens[k,1,],
                                                              target_sum = 100, digits = 0)

                                                }))
    rownames(pcalib_postmean_percentage) = dimnames(pcalib_postsumm_asarray_wens)[[1]]
    colnames(pcalib_postmean_percentage) = dimnames(pcalib_postsumm_asarray_wens)[[3]]
    # print(dim(pcalib_postmean_percentage))


    # output list
    # print(names(input_vacalib))
    input_vacalib.now = as.list(environment())
    # print(names(input_vacalib.now))
    # print(names(input_vacalib.now)[names(input_vacalib) %in% names(input_vacalib.now)])
    input_vacalib[intersect(names(input_vacalib), names(input_vacalib.now))] = input_vacalib.now[intersect(names(input_vacalib), names(input_vacalib.now))]
    # print(names(input_vacalib))

    output = list("calib_MCMCout" = calibout,
                  "p_uncalib" = puncalib,
                  "p_calib" = pcalib_asarray_wens,
                  "pcalib_postsumm" = pcalib_postsumm_asarray_wens,
                  "puncalib_percentage" = puncalib_percentage,
                  "pcalib_postmean_percentage" = pcalib_postmean_percentage,
                  "va_deaths_uncalib" = va_deaths,
                  "va_deaths_calib_algo" = va_deaths_calib,
                  "va_deaths_calib_ensemble" = va_deaths_ens,
                  'Mmat_input' = Mmat_calib_input_asarray,
                  'Mmat_study' = Mmat_calib_study_asarray,
                  'Mmat_tomodel' = Mmat_calib_tomodel_asarray,
                  'donotcalib_study' = donotcalib_study_asmat,
                  'donotcalib_tomodel' = donotcalib_tomodel_asmat,
                  'calibrated' = calibrated,
                  'lambda_calibpath' = lambda_calibpath,
                  "K" = K, "nCause" = nCause, "causes" = causes,
                  'input' = input_vacalib)

    # print(dim(Mmat_input_asarray))
    # print(dim(Mmat_asarray))
    # print(dim(output$Mmat.fixed_input))
    # print(dim(output$Mmat.fixed))

    if(verbose){

      # message("Ensemble calibration ... Done.")
      message("")

    }

    # end ensemble

  }

  # print("calibrated succesfully")

  return(output)

}
