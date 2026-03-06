#' PRECOG Prognostic Z-Scores
#'
#' Gene-level prognostic z-scores from the PRECOG (PREdiction of Clinical 
#' Outcomes from Genomic Profiles) meta-analysis.
#'
#' @format A data.frame with genes as rows and cancer types as columns.
#'   Values represent meta-z-scores indicating prognostic association.
#'
#' @source Gentles et al. (2015) Cell. DOI: 10.1016/j.cell.2015.11.025
#'
"precog"


#' TCGA Prognostic Z-Scores
#'
#' Gene-level prognostic z-scores derived from The Cancer Genome Atlas (TCGA) 
#' survival analyses.
#'
#' @format A data.frame with genes as rows and cancer types as columns.
#'   Values represent z-scores from Cox proportional hazards models.
#'
#' @source The Cancer Genome Atlas Research Network
#'
"tcga"


#' Pediatric Cancer Prognostic Z-Scores
#'
#' Gene-level prognostic z-scores for pediatric cancers.
#'
#' @format A data.frame with genes as rows and pediatric cancer types as columns.
#'
"pediatric"


#' Immune Checkpoint Inhibitor (ICI) Prognostic Z-Scores
#'
#' Gene-level prognostic z-scores for patients treated with immune checkpoint 
#' inhibitors, separated by primary and metastatic disease.
#'
#' @format A data.frame with genes as rows and columns in format:
#'   CANCER_THERAPY_STAGE_ENDPOINT (e.g., "MELANOMA_Anti-PD1_Primary_OS")
#'
"ici"
