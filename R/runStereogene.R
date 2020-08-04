#' @title runStereogene
#'
#' @description Writes a StereoGene script in the working directory
#'
#' @param track_files Vector of at least two track or interval file names to be
#' pairwise-correlated by StereoGene. Required.
#' @param name_config Name of corresponding configuration file. Required
#' @param pcorProfile Track for partial correlation. More information for this
#' can be found in the StereoGene README. Optional
#' @param confounder Confounder filename. More information for this can be found
#' in the StereoGene README. Optional
#'
#' @export

runStereogene <- function(track_files,
                            name_config,
                            pcorProfile = NULL,
                            confounder = NULL) {
    if (length(track_files) < 2) {
        stop("Must have at least two track or interval files for correlation.")
    }
    if (length(name_config) < 1) {
        stop("Must provide a name for the configuration file.")
    }
    tracks <- ""
    for (track in track_files) {
        tracks <- paste0(tracks, " ", track)
    }
    config <- paste0("-cfg ", name_config, " ")
    partial_corr <- ""
    if (length(pcorProfile) == 1) {
        partial_corr <- paste0("-pcorProfile ", pcorProfile, " ")
    }
    confound <- ""
    if (length(confounder) == 1) {
        confound <- paste0("-confounder ", confounder, " ")
    }
    system2("StereoGene", paste0(config, partial_corr, confound, tracks)
    )
}
