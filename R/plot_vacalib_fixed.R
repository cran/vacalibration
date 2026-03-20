#' Summary Plots of VA-Calibration Using Fixed Misclassification Matrix
#'
#' This is a utility function. Please use \link{plot_vacalib}.
#'
#' @param vacalib_fit Fitted object from \code{vacalibration()}
#' @param toplot Character. Same as \code{toplot} in \code{plot_vacalib_fixed()}
#' @return Plots misclassification matrices and/or cause-specific mortality fractions
#'
#' @aliases plot_vacalib_fixed
#'
#' @import patchwork
#' @import ggplot2
#'
#' @importFrom stats quantile
#' @importFrom utils head
#'
plot_vacalib_fixed <- function(vacalib_fit, toplot){

  if(toplot=="both"){

    # plot both misclassification and csmf ----
    ## input misclassification ----
    # print(round(100*vacalib_fit$Mmat.fixed_input[1,,]))
    # Mmat_toplot_input = round(100*vacalib_fit$Mmat_input)
    # print(dim(Mmat_toplot_input))
    Mmat_toplot_input = vacalib_fit$Mmat_input
    for(k in 1:dim(Mmat_toplot_input)[1]){

      Mmat_toplot_input[k,,] = do.call("rbind",
                                       lapply(1:dim(Mmat_toplot_input)[2],
                                              FUN = function(i){

                                                smart_round(x = 100*vacalib_fit$Mmat_input[k,i,],
                                                            target_sum = 100, digits = 0)

                                              }))

    }

    plotdf_Mmat_input = reshape2::melt(Mmat_toplot_input)
    head(plotdf_Mmat_input)
    value.labels_Mmat_input = plotdf_Mmat_input$value

    plotdf_Mmat_input$Var1 = factor(x = plotdf_Mmat_input$Var1, levels = dimnames(Mmat_toplot_input)[[1]])
    plotdf_Mmat_input$Var2 = factor(x = plotdf_Mmat_input$Var2, levels = rev(dimnames(vacalib_fit$Mmat_input)[[2]]))
    plotdf_Mmat_input$Var3 = factor(x = plotdf_Mmat_input$Var3, levels = dimnames(vacalib_fit$Mmat_input)[[2]])

    plotdf_Mmat_input$diag = plotdf_Mmat_input$Var2==plotdf_Mmat_input$Var3
    plotdf_Mmat_input$diag[!plotdf_Mmat_input$diag] = NA
    # print(head(plotdf_Mmat_input))
    # print(K)

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
        plot.subtitle = ggplot2::element_text(),
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
      )
    # print(ggplot2_Mmat_input)


    ## study-specific misclassification ----
    Mmat_toplot_study = vacalib_fit$Mmat_study
    # print(dim(Mmat_toplot_study))
    for(k in 1:dim(Mmat_toplot_study)[1]){

      if(!vacalib_fit$calibrated[k]){

        Mmat_toplot_study[k,,] = NA*Mmat_toplot_study[k,,]

      }else{

        Mmat_toplot_study[k,,] = diag(dim(Mmat_toplot_study)[2])
        Mmat_toplot_study[k,!vacalib_fit$donotcalib_tomodel[k,],!vacalib_fit$donotcalib_tomodel[k,]] =
          do.call("rbind",
                  lapply(which(!vacalib_fit$donotcalib_tomodel[k,]),
                         FUN = function(i){

                           smart_round(x = 100*vacalib_fit$Mmat_study[k,i,!vacalib_fit$donotcalib_tomodel[k,]],
                                       target_sum = 100, digits = 0)

                         }))

        Mmat_toplot_study[k,vacalib_fit$donotcalib_tomodel[k,],] =
          Mmat_toplot_study[k,,vacalib_fit$donotcalib_tomodel[k,]] = NA

      }

    }

    plotdf_Mmat_study = reshape2::melt(Mmat_toplot_study)
    head(plotdf_Mmat_study)
    value.labels_Mmat_study = plotdf_Mmat_study$value

    plotdf_Mmat_study$Var1 = factor(x = plotdf_Mmat_study$Var1, levels = dimnames(Mmat_toplot_study)[[1]])
    plotdf_Mmat_study$Var2 = factor(x = plotdf_Mmat_study$Var2, levels = rev(dimnames(Mmat_toplot_study)[[2]]))
    plotdf_Mmat_study$Var3 = factor(x = plotdf_Mmat_study$Var3, levels = dimnames(Mmat_toplot_study)[[2]])

    plotdf_Mmat_study$diag = plotdf_Mmat_study$Var2==plotdf_Mmat_study$Var3
    plotdf_Mmat_study$diag[!plotdf_Mmat_study$diag] = NA
    # print(head(plotdf_Mmat_study))
    # print(K)

    ggplot2_Mmat_study =
      ggplot2::ggplot(plotdf_Mmat_study, ggplot2::aes(Var3, Var2, fill = value)) +
      ggplot2::geom_tile(color="white", linewidth=.5) +
      ggplot2::geom_tile(data = plotdf_Mmat_study[!is.na(plotdf_Mmat_study$diag), ],
                         ggplot2::aes(Var3, Var2, fill = value, color = diag), linewidth = .7) +
      ggplot2::scale_color_manual(guide = "none", values = c(`TRUE` = "blue3")) +
      ggplot2::geom_text(ggplot2::aes(Var3, Var2, label = value.labels_Mmat_study,
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
        plot.subtitle = ggplot2::element_text(),
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
      )
    # print(ggplot2_Mmat_study)


    ## misclassification used for calibration ----
    Mmat_toplot_tomodel = vacalib_fit$Mmat_tomodel
    # print(dim(Mmat_toplot_tomodel))
    for(k in 1:dim(Mmat_toplot_tomodel)[1]){

      if(!vacalib_fit$calibrated[k]){

        Mmat_toplot_tomodel[k,,] = NA*Mmat_toplot_tomodel[k,,]

      }else{

        Mmat_toplot_tomodel[k,,] = diag(dim(Mmat_toplot_tomodel)[2])
        Mmat_toplot_tomodel[k,!vacalib_fit$donotcalib_tomodel[k,],!vacalib_fit$donotcalib_tomodel[k,]] =
          do.call("rbind",
                  lapply(which(!vacalib_fit$donotcalib_tomodel[k,]),
                         FUN = function(i){

                           smart_round(x = 100*vacalib_fit$Mmat_tomodel[k,i,!vacalib_fit$donotcalib_tomodel[k,]],
                                       target_sum = 100, digits = 0)

                         }))

      }

      Mmat_toplot_tomodel[k,vacalib_fit$donotcalib_tomodel[k,],] =
        Mmat_toplot_tomodel[k,,vacalib_fit$donotcalib_tomodel[k,]] = NA

    }

    plotdf_Mmat_tomodel = reshape2::melt(Mmat_toplot_tomodel)
    head(plotdf_Mmat_tomodel)
    value.labels_Mmat_tomodel = plotdf_Mmat_tomodel$value

    plotdf_Mmat_tomodel$Var1 = factor(x = plotdf_Mmat_tomodel$Var1, levels = dimnames(Mmat_toplot_tomodel)[[1]])
    plotdf_Mmat_tomodel$Var2 = factor(x = plotdf_Mmat_tomodel$Var2, levels = rev(dimnames(Mmat_toplot_tomodel)[[2]]))
    plotdf_Mmat_tomodel$Var3 = factor(x = plotdf_Mmat_tomodel$Var3, levels = dimnames(Mmat_toplot_tomodel)[[2]])

    plotdf_Mmat_tomodel$diag = plotdf_Mmat_tomodel$Var2==plotdf_Mmat_tomodel$Var3
    plotdf_Mmat_tomodel$diag[!plotdf_Mmat_tomodel$diag] = NA
    # print(head(plotdf_Mmat_tomodel))
    # print(K)

    ggplot2_Mmat_tomodel =
      ggplot2::ggplot(plotdf_Mmat_tomodel, ggplot2::aes(Var3, Var2, fill = value)) +
      ggplot2::geom_tile(color="white", linewidth=.5) +
      ggplot2::geom_tile(data = plotdf_Mmat_tomodel[!is.na(plotdf_Mmat_tomodel$diag), ],
                         ggplot2::aes(Var3, Var2, fill = value, color = diag), linewidth = .7) +
      ggplot2::scale_color_manual(guide = "none", values = c(`TRUE` = "blue3")) +
      ggplot2::geom_text(ggplot2::aes(Var3, Var2, label = value.labels_Mmat_tomodel,
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
        plot.subtitle = ggplot2::element_text(),
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
      )
    # print(ggplot2_Mmat_tomodel)



    ## calibrated vs uncalibrated csmf ----
    plotdf_pcalib = NULL
    for(k in 1:dim(vacalib_fit$pcalib_postsumm)[1]){

      plotdf_pcalib = rbind.data.frame(plotdf_pcalib,
                                       cbind.data.frame(rbind.data.frame(data.frame('causes' = dimnames(vacalib_fit$pcalib_postsumm)[[3]],
                                                                                    'value' = unname(vacalib_fit$pcalib_postsumm[k,'postmean',]),
                                                                                    'llim' = unname(vacalib_fit$pcalib_postsumm[k,'lowcredI',]),
                                                                                    'ulim' = unname(vacalib_fit$pcalib_postsumm[k,'upcredI',]),
                                                                                    'calib_type' = 'Calibrated'),
                                                                         data.frame('causes' = colnames(vacalib_fit$p_uncalib),
                                                                                    'value' = unname(vacalib_fit$p_uncalib[k,]),
                                                                                    'llim' = NA,
                                                                                    'ulim' = NA,
                                                                                    'calib_type' = 'Uncalibrated')),
                                                        'vaalgo' = (rownames(vacalib_fit$p_uncalib))[k]))

    }
    head(plotdf_pcalib)

    plotdf_pcalib$causes = factor(x = plotdf_pcalib$causes,
                                  levels = colnames(vacalib_fit$p_uncalib))
    plotdf_pcalib$calib_type = factor(x = plotdf_pcalib$calib_type,
                                      levels = c("Uncalibrated", "Calibrated"))
    plotdf_pcalib$vaalgo = factor(x = plotdf_pcalib$vaalgo,
                                  levels = rownames(vacalib_fit$p_uncalib))
    head(plotdf_pcalib)

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
      ggplot2::labs(title = 'Cause-Specific Mortality Fractions (CSMF)',
                    x = "Cause",
                    y = 'Estimate')
    # ggplot2_pcalib


    ## printing plots ----
    if(isTRUE(all.equal(vacalib_fit$K,1))){

      if(!is.null(vacalib_fit$input$studycause_map)){

        if(vacalib_fit$input$path_correction){

          ggplot2_Mmat_input = ggplot2_Mmat_input +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_study = ggplot2_Mmat_study +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Study-Specific",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          print((ggplot2_Mmat_input | ggplot2_Mmat_study | ggplot2_Mmat_tomodel | ggplot2_pcalib))

        }else{

          ggplot2_Mmat_input = ggplot2_Mmat_input +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          print((ggplot2_Mmat_input | ggplot2_Mmat_tomodel | ggplot2_pcalib))

        }

      }else{

        if(vacalib_fit$input$path_correction){

          ggplot2_Mmat_input = ggplot2_Mmat_input +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          print((ggplot2_Mmat_input | ggplot2_Mmat_tomodel | ggplot2_pcalib))

        }else{

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input and Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          print((ggplot2_Mmat_tomodel | ggplot2_pcalib))

        }

      }

    }else if(vacalib_fit$K>1){

      if(!is.null(vacalib_fit$input$studycause_map)){

        if(vacalib_fit$input$path_correction){

          ggplot2_Mmat_input = ggplot2_Mmat_input +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_study = ggplot2_Mmat_study +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Study-Specific",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          if(vacalib_fit$input$ensemble){

            print(((ggplot2_Mmat_input | NULL) + plot_layout(widths = c(vacalib_fit$K, 1))) / ((ggplot2_Mmat_study | NULL) + plot_layout(widths = c(vacalib_fit$K, 1))) / ((ggplot2_Mmat_tomodel | NULL) + plot_layout(widths = c(vacalib_fit$K, 1))) / ggplot2_pcalib)

          }else{

            print(ggplot2_Mmat_input / ggplot2_Mmat_study / ggplot2_Mmat_tomodel / ggplot2_pcalib)

          }

        }else{

          ggplot2_Mmat_input = ggplot2_Mmat_input +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          if(vacalib_fit$input$ensemble){

            print(((ggplot2_Mmat_input | NULL) + plot_layout(widths = c(vacalib_fit$K, 1))) / ((ggplot2_Mmat_tomodel | NULL) + plot_layout(widths = c(vacalib_fit$K, 1))) / ggplot2_pcalib)

          }else{

            print(ggplot2_Mmat_input / ggplot2_Mmat_tomodel / ggplot2_pcalib)

          }

        }

      }else{

        if(vacalib_fit$input$path_correction){

          ggplot2_Mmat_input = ggplot2_Mmat_input +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          if(vacalib_fit$input$ensemble){

            print(((ggplot2_Mmat_input | NULL) + plot_layout(widths = c(vacalib_fit$K, 1))) / ((ggplot2_Mmat_tomodel | NULL) + plot_layout(widths = c(vacalib_fit$K, 1))) / ggplot2_pcalib)

          }else{

            print(ggplot2_Mmat_input / ggplot2_Mmat_tomodel/ ggplot2_pcalib)

          }

        }else{

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input and Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          if(vacalib_fit$input$ensemble){

            print(((ggplot2_Mmat_tomodel | NULL) + plot_layout(widths = c(vacalib_fit$K, 1))) / ggplot2_pcalib)

          }else{

            print(ggplot2_Mmat_tomodel / ggplot2_pcalib)

          }

        }

      }

    }

  }else if(toplot=="missmat"){

    # plot only misclassification ----
    ## input misclassification ----
    # print(round(100*vacalib_fit$Mmat.fixed_input[1,,]))
    # Mmat_toplot_input = round(100*vacalib_fit$Mmat_input)
    # print(dim(Mmat_toplot_input))
    Mmat_toplot_input = vacalib_fit$Mmat_input
    for(k in 1:dim(Mmat_toplot_input)[1]){

      Mmat_toplot_input[k,,] = do.call("rbind",
                                       lapply(1:dim(Mmat_toplot_input)[2],
                                              FUN = function(i){

                                                smart_round(x = 100*vacalib_fit$Mmat_input[k,i,],
                                                            target_sum = 100, digits = 0)

                                              }))

    }

    plotdf_Mmat_input = reshape2::melt(Mmat_toplot_input)
    head(plotdf_Mmat_input)
    value.labels_Mmat_input = plotdf_Mmat_input$value

    plotdf_Mmat_input$Var1 = factor(x = plotdf_Mmat_input$Var1, levels = dimnames(Mmat_toplot_input)[[1]])
    plotdf_Mmat_input$Var2 = factor(x = plotdf_Mmat_input$Var2, levels = rev(dimnames(vacalib_fit$Mmat_input)[[2]]))
    plotdf_Mmat_input$Var3 = factor(x = plotdf_Mmat_input$Var3, levels = dimnames(vacalib_fit$Mmat_input)[[2]])

    plotdf_Mmat_input$diag = plotdf_Mmat_input$Var2==plotdf_Mmat_input$Var3
    plotdf_Mmat_input$diag[!plotdf_Mmat_input$diag] = NA
    # print(head(plotdf_Mmat_input))
    # print(K)

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
        plot.subtitle = ggplot2::element_text(),
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
      )
    # print(ggplot2_Mmat_input)


    ## study-specific misclassification ----
    Mmat_toplot_study = vacalib_fit$Mmat_study
    # print(dim(Mmat_toplot_study))
    for(k in 1:dim(Mmat_toplot_study)[1]){

      if(!vacalib_fit$calibrated[k]){

        Mmat_toplot_study[k,,] = NA*Mmat_toplot_study[k,,]

      }else{

        Mmat_toplot_study[k,,] = diag(dim(Mmat_toplot_study)[2])
        Mmat_toplot_study[k,!vacalib_fit$donotcalib_tomodel[k,],!vacalib_fit$donotcalib_tomodel[k,]] =
          do.call("rbind",
                  lapply(which(!vacalib_fit$donotcalib_tomodel[k,]),
                         FUN = function(i){

                           smart_round(x = 100*vacalib_fit$Mmat_study[k,i,!vacalib_fit$donotcalib_tomodel[k,]],
                                       target_sum = 100, digits = 0)

                         }))

        Mmat_toplot_study[k,vacalib_fit$donotcalib_tomodel[k,],] =
          Mmat_toplot_study[k,,vacalib_fit$donotcalib_tomodel[k,]] = NA

      }

    }

    plotdf_Mmat_study = reshape2::melt(Mmat_toplot_study)
    head(plotdf_Mmat_study)
    value.labels_Mmat_study = plotdf_Mmat_study$value

    plotdf_Mmat_study$Var1 = factor(x = plotdf_Mmat_study$Var1, levels = dimnames(Mmat_toplot_study)[[1]])
    plotdf_Mmat_study$Var2 = factor(x = plotdf_Mmat_study$Var2, levels = rev(dimnames(Mmat_toplot_study)[[2]]))
    plotdf_Mmat_study$Var3 = factor(x = plotdf_Mmat_study$Var3, levels = dimnames(Mmat_toplot_study)[[2]])

    plotdf_Mmat_study$diag = plotdf_Mmat_study$Var2==plotdf_Mmat_study$Var3
    plotdf_Mmat_study$diag[!plotdf_Mmat_study$diag] = NA
    # print(head(plotdf_Mmat_study))
    # print(K)

    ggplot2_Mmat_study =
      ggplot2::ggplot(plotdf_Mmat_study, ggplot2::aes(Var3, Var2, fill = value)) +
      ggplot2::geom_tile(color="white", linewidth=.5) +
      ggplot2::geom_tile(data = plotdf_Mmat_study[!is.na(plotdf_Mmat_study$diag), ],
                         ggplot2::aes(Var3, Var2, fill = value, color = diag), linewidth = .7) +
      ggplot2::scale_color_manual(guide = "none", values = c(`TRUE` = "blue3")) +
      ggplot2::geom_text(ggplot2::aes(Var3, Var2, label = value.labels_Mmat_study,
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
        plot.subtitle = ggplot2::element_text(),
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
      )
    # print(ggplot2_Mmat_study)


    ## misclassification used for calibration ----
    Mmat_toplot_tomodel = vacalib_fit$Mmat_tomodel
    # print(dim(Mmat_toplot_tomodel))
    for(k in 1:dim(Mmat_toplot_tomodel)[1]){

      if(!vacalib_fit$calibrated[k]){

        Mmat_toplot_tomodel[k,,] = NA*Mmat_toplot_tomodel[k,,]

      }else{

        Mmat_toplot_tomodel[k,,] = diag(dim(Mmat_toplot_tomodel)[2])
        Mmat_toplot_tomodel[k,!vacalib_fit$donotcalib_tomodel[k,],!vacalib_fit$donotcalib_tomodel[k,]] =
          do.call("rbind",
                  lapply(which(!vacalib_fit$donotcalib_tomodel[k,]),
                         FUN = function(i){

                           smart_round(x = 100*vacalib_fit$Mmat_tomodel[k,i,!vacalib_fit$donotcalib_tomodel[k,]],
                                       target_sum = 100, digits = 0)

                         }))

      }

      Mmat_toplot_tomodel[k,vacalib_fit$donotcalib_tomodel[k,],] =
        Mmat_toplot_tomodel[k,,vacalib_fit$donotcalib_tomodel[k,]] = NA

    }

    plotdf_Mmat_tomodel = reshape2::melt(Mmat_toplot_tomodel)
    head(plotdf_Mmat_tomodel)
    value.labels_Mmat_tomodel = plotdf_Mmat_tomodel$value

    plotdf_Mmat_tomodel$Var1 = factor(x = plotdf_Mmat_tomodel$Var1, levels = dimnames(Mmat_toplot_tomodel)[[1]])
    plotdf_Mmat_tomodel$Var2 = factor(x = plotdf_Mmat_tomodel$Var2, levels = rev(dimnames(Mmat_toplot_tomodel)[[2]]))
    plotdf_Mmat_tomodel$Var3 = factor(x = plotdf_Mmat_tomodel$Var3, levels = dimnames(Mmat_toplot_tomodel)[[2]])

    plotdf_Mmat_tomodel$diag = plotdf_Mmat_tomodel$Var2==plotdf_Mmat_tomodel$Var3
    plotdf_Mmat_tomodel$diag[!plotdf_Mmat_tomodel$diag] = NA
    # print(head(plotdf_Mmat_tomodel))
    # print(K)

    ggplot2_Mmat_tomodel =
      ggplot2::ggplot(plotdf_Mmat_tomodel, ggplot2::aes(Var3, Var2, fill = value)) +
      ggplot2::geom_tile(color="white", linewidth=.5) +
      ggplot2::geom_tile(data = plotdf_Mmat_tomodel[!is.na(plotdf_Mmat_tomodel$diag), ],
                         ggplot2::aes(Var3, Var2, fill = value, color = diag), linewidth = .7) +
      ggplot2::scale_color_manual(guide = "none", values = c(`TRUE` = "blue3")) +
      ggplot2::geom_text(ggplot2::aes(Var3, Var2, label = value.labels_Mmat_tomodel,
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
        plot.subtitle = ggplot2::element_text(),
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
      )
    # print(ggplot2_Mmat_tomodel)


    ## printing plots ----
    if(isTRUE(all.equal(vacalib_fit$K,1))){

      if(!is.null(vacalib_fit$input$studycause_map)){

        if(vacalib_fit$input$path_correction){

          ggplot2_Mmat_input = ggplot2_Mmat_input +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_study = ggplot2_Mmat_study +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Study-Specific",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          print((ggplot2_Mmat_input | ggplot2_Mmat_study | ggplot2_Mmat_tomodel))

        }else{

          ggplot2_Mmat_input = ggplot2_Mmat_input +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          print((ggplot2_Mmat_input | ggplot2_Mmat_tomodel))

        }

      }else{

        if(vacalib_fit$input$path_correction){

          ggplot2_Mmat_input = ggplot2_Mmat_input +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          print((ggplot2_Mmat_input | ggplot2_Mmat_tomodel))

        }else{

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input and Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          print(ggplot2_Mmat_tomodel)

        }

      }

    }else if(vacalib_fit$K>1){

      if(!is.null(vacalib_fit$input$studycause_map)){

        if(vacalib_fit$input$path_correction){

          ggplot2_Mmat_input = ggplot2_Mmat_input +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_study = ggplot2_Mmat_study +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Study-Specific",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          print(ggplot2_Mmat_input / ggplot2_Mmat_study / ggplot2_Mmat_tomodel / ggplot2_pcalib)

        }else{

          ggplot2_Mmat_input = ggplot2_Mmat_input +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          print(ggplot2_Mmat_input / ggplot2_Mmat_tomodel)

        }

      }else{

        if(vacalib_fit$input$path_correction){

          ggplot2_Mmat_input = ggplot2_Mmat_input +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          print(ggplot2_Mmat_input / ggplot2_Mmat_tomodel)

        }else{

          ggplot2_Mmat_tomodel = ggplot2_Mmat_tomodel +
            ggplot2::labs(title = "Fixed Misclassification Matrix",
                          subtitle = "Input and Used For Calibration",
                          x = 'VA Cause', y = 'CHAMPS Cause')

          print(ggplot2_Mmat_tomodel)

        }

      }

    }

  }else if(toplot=="csmf"){

    # plot only csmf ----
    ## calibrated vs uncalibrated csmf ----
    plotdf_pcalib = NULL
    for(k in 1:dim(vacalib_fit$pcalib_postsumm)[1]){

      plotdf_pcalib = rbind.data.frame(plotdf_pcalib,
                                       cbind.data.frame(rbind.data.frame(data.frame('causes' = dimnames(vacalib_fit$pcalib_postsumm)[[3]],
                                                                                    'value' = unname(vacalib_fit$pcalib_postsumm[k,'postmean',]),
                                                                                    'llim' = unname(vacalib_fit$pcalib_postsumm[k,'lowcredI',]),
                                                                                    'ulim' = unname(vacalib_fit$pcalib_postsumm[k,'upcredI',]),
                                                                                    'calib_type' = 'Calibrated'),
                                                                         data.frame('causes' = colnames(vacalib_fit$p_uncalib),
                                                                                    'value' = unname(vacalib_fit$p_uncalib[k,]),
                                                                                    'llim' = NA,
                                                                                    'ulim' = NA,
                                                                                    'calib_type' = 'Uncalibrated')),
                                                        'vaalgo' = (rownames(vacalib_fit$p_uncalib))[k]))

    }
    head(plotdf_pcalib)

    plotdf_pcalib$causes = factor(x = plotdf_pcalib$causes,
                                  levels = colnames(vacalib_fit$p_uncalib))
    plotdf_pcalib$calib_type = factor(x = plotdf_pcalib$calib_type,
                                      levels = c("Uncalibrated", "Calibrated"))
    plotdf_pcalib$vaalgo = factor(x = plotdf_pcalib$vaalgo,
                                  levels = rownames(vacalib_fit$p_uncalib))
    head(plotdf_pcalib)

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
      ggplot2::labs(title = 'Cause-Specific Mortality Fractions (CSMF)',
                    x = "Cause",
                    y = 'Estimate')
    # ggplot2_pcalib


    ## printing plots ----
    print(ggplot2_pcalib)

  }

}
