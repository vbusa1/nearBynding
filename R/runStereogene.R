#' @title runStereogene
#'
#' @description Writes a StereoGene script in the working directory
#'
#' @param track_files Vector of at least two track or interval file names to be
#' pairwise-correlated by StereoGene. Required.
#' @param name_config Name of corresponding configuration file; a string.
#' Required
#' @param pcorProfile Name of track file name for partial correlation; a string.
#' More information for this can be found in the StereoGene README. Optional
#' @param confounder Confounder file name; a string. More information for this
#' can be found in the StereoGene README. Optional
#' @param nShuffle Permutations used to estimate error. Default 5000.
#'
#' @return generates StereoGene output files in directory
#'
#' @examples
#' runStereogene(track_files = c("chr4and5_3UTR_stem_liftOver.bedGraph",
#'                              "chr4and5_liftOver.bedGraph"),
#'              name_config = "chr4and5_3UTR.cfg")
#'
#' @export

runStereogene <- function(track_files,
                            name_config,
                            pcorProfile = NULL,
                            confounder = NULL,
                            nShuffle = 5000) {
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
    if(.is_StereoGene_installed()){
        .StereoGene_run(paste0(config, partial_corr, confound, tracks))
    }else{return("Please install Stereogene and place in working PATH")}
}
