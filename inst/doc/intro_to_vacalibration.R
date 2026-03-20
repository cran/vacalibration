## ----results='asis', echo=FALSE-----------------------------------------------
cat('
<style>
  pre {
    position: relative;
  }
  .copy-btn {
    position: absolute;
    top: 6px;
    right: 6px;
    background-color: #f0f0f0;
    border: 1px solid #ccc;
    border-radius: 4px;
    padding: 4px 8px;
    font-size: 12px;
    cursor: pointer;
    opacity: 0.7;
    transition: opacity 0.3s ease;
  }
  .copy-btn:hover {
    opacity: 1;
  }
</style>

<script>
document.addEventListener("DOMContentLoaded", function() {
  document.querySelectorAll("pre").forEach(function(pre) {
    const button = document.createElement("button");
    button.className = "copy-btn";
    button.textContent = "Copy";
    pre.style.position = "relative";
    pre.appendChild(button);

    button.addEventListener("click", function() {
      const code = pre.querySelector("code");
      if (!code) return;
      navigator.clipboard.writeText(code.innerText).then(() => {
        button.textContent = "Copied!";
        setTimeout(() => {
          button.textContent = "Copy";
        }, 1500);
      });
    });
  });
});
</script>
')

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=F-------------------------------------------------------------------
# install.packages("vacalibration")
# library(vacalibration) # load

## ----eval=F-------------------------------------------------------------------
# # install "devtools" R package
# devtools::install_github("sandy-pramanik/vacalibration")
# library(vacalibration) # load

## ----eval=F-------------------------------------------------------------------
# data("comsamoz_CCVAoutput")
# 
# comsamoz_CCVAoutput$neonate$eava  # output from EAVA for neonates
# comsamoz_CCVAoutput$neonate$insilicova  # output from InSilicoVA for neonates
# comsamoz_CCVAoutput$neonate  # list of outputs for neonates from EAVA, InSilicoVA, and InterVA

## ----eval=F-------------------------------------------------------------------
# data("CCVA_missmat")

## ----eval=F-------------------------------------------------------------------
# CCVA_missmat$neonate$eava$postmean$Mozambique  # average
# CCVA_missmat$neonate$eava$asDirich$Mozambique  # Dirichlet approximation
# CCVA_missmat$neonate$eava$postsumm$Mozambique  # summary of distribution

## ----eval=F, results = 'hide', message = FALSE, warning = FALSE, fig.show = "hide"----
# vacalib_eava = vacalibration(va_data = list("eava" = comsamoz_CCVAoutput$neonate$eava),
#                              age_group = "neonate", country = "Mozambique")
# 
# # CSMF
# vacalib_eava$p_uncalib[1,]  # uncalibrated estimates
# vacalib_eava$p_calib[1,,]  # posterior of calibrated estimates
# vacalib_eava$pcalib_postsumm[1,,]  # posterior summary of calibrated estimates
# 
# # death counts
# vacalib_eava$va_deaths_uncalib[1,]  # uncalibrated
# vacalib_eava$va_deaths_calib_algo[1,]  # calibrated

## ----eval=F, results = 'hide', message = FALSE, warning = FALSE, fig.show = "hide"----
# vacalib_ensemble =
#   vacalibration(va_data = list("eava" = comsamoz_CCVAoutput$neonate$eava,
#                                "insilicova" = comsamoz_CCVAoutput$neonate$insilicova,
#                                "interva" = comsamoz_CCVAoutput$neonate$interva),
#                 age_group = "neonate", country = "Mozambique")
# 
# # CSMF
# vacalib_ensemble$p_uncalib  # uncalibrated estimates
# 
# # posterior of calibrated CSMF
# vacalib_ensemble$p_calib["eava",,]  # EAVA
# vacalib_ensemble$p_calib["insilicova",,]  # InSilicoVA
# vacalib_ensemble$p_calib["interva",,]  # InterVA
# vacalib_ensemble$p_calib["ensemble",,]  # ensemble
# 
# # posterior summary of calibrated CSMF
# vacalib_ensemble$pcalib_postsumm["eava",,]  # EAVA
# vacalib_ensemble$pcalib_postsumm["insilicova",,]  # InSilicoVA
# vacalib_ensemble$pcalib_postsumm["interva",,]  # InterVA
# vacalib_ensemble$pcalib_postsumm["ensemble",,]  # ensemble
# 
# # death counts
# vacalib_ensemble$va_deaths_uncalib  # uncalibrated
# vacalib_ensemble$va_deaths_calib_algo  # calibrated counts from algorithm-specific calibration
# vacalib_ensemble$va_deaths_calib_ensemble  # calibrated counts from ensemble calibration

## ----eval=F, results = 'hide', message = FALSE, warning = FALSE, fig.show = "hide"----
# plot_vacalib(vacalib_fit = vacalib_eava)

## ----echo=FALSE, out.width="100%"---------------------------------------------
knitr::include_graphics("figures/vacalib_eava_plot.png")

## ----eval=F, results = 'hide', message = FALSE, warning = FALSE, fig.show = "hide"----
# set_studycause_map = c("Intrapartum" = "ipre", "Congenital" = "congenital_malformation",
#                        "Diarrhoeal" = "sepsis_meningitis_inf", "LRI" = "pneumonia",
#                        "Sepsis" = "sepsis_meningitis_inf", "Preterm" = "prematurity",
#                        "Tetanus" = "sepsis_meningitis_inf", "Other" = "other")

## ----eval=F, results = 'hide', message = FALSE, warning = FALSE, fig.show = "hide"----
# 
# 
# vacalib_cacode = vacalibration(va_data = list("eava" = c("Intrapartum" = 82, "Congenital" = 17,
#                                                          "Diarrhoeal" = 6, "LRI" = 33,
#                                                          "Sepsis" = 108, "Preterm" = 35,
#                                                          "Tetanus" = 14, "Other" = 7)),
#                                age_group = "neonate", country = "Bangladesh",
#                                studycause_map = set_studycause_map)
# 
# # CSMF
# vacalib_cacode$p_uncalib[1,]  # uncalibrated estimates
# vacalib_cacode$p_calib[1,,]  # posterior of calibrated estimates
# vacalib_cacode$pcalib_postsumm[1,,]  # posterior summary of calibrated estimates
# 
# # death counts
# vacalib_cacode$va_deaths_uncalib[1,]  # uncalibrated
# vacalib_cacode$va_deaths_calib_algo[1,]  # calibrated

