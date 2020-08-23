#' @title runStereogeneOnCapR
#'
#' @description Writes a configuration file and Stereogene script and runs
#' Stereogene for all CapR tracks
#'
#' @param dir_CapR_bg Directory of lifted-over CapR bedGraph files. Default
#' current directory
#' @param input_prefix Prefix string appended to input files; same as
#' input_prefix argument in processCapRout. Required
#' @param protein_file Name of protein file in bedGraph format. Required
#' @param output_prefix Prefix string to be appended to all output files.
#' Default to be same as input_prefix
#' @param chrom_size Name of chromosome size file. File must be in two-column
#' format without a header where first column is chromosome name and second
#' column is chromosome length, as from getChainChrSize. Required
#' @param name_config Name of output config file. Default config.cfg
#' @param ... includes all other parameters acceptable to write_config and
#' write_stereogene
#'
#' @return generates StereoGene output files, including *.dist files
#'
#' @examples
#' \donttest{
#'runStereogeneOnCapR(protein_file = "chr4and5_liftOver.bedGraph",
#'                    chrom_size = "chr4and5_3UTR.size",
#'                    name_config = "chr4and5_3UTR.cfg",
#'                    input_prefix = "chr4and5_3UTR")
#'}
#'
#' @importFrom R.utils doCall
#'
#' @export

runStereogeneOnCapR <- function(dir_CapR_bg = ".",
                                input_prefix,
                                protein_file,
                                output_prefix = input_prefix,
                                name_config = "config.cfg",
                                chrom_size,
                                ...) {
    if (missing(chrom_size)){stop("please provide a chrom_size file")}
    doCall(write_config, args = list(
        name_config = name_config,
        chrom_size = chrom_size, ...))
    track_files <- c(paste0(dir_CapR_bg, "/", input_prefix,
        "_hairpin_liftOver.bedGraph"), protein_file)
    name_sh <- paste0(output_prefix, "_stereogene_hairpin.sh")
    doCall(runStereogene, args = list(...,
        track_files = track_files,
        name_sh = name_sh,
        name_config = name_config))
    track_files <- c(paste0(dir_CapR_bg, "/", input_prefix,
        "_internal_liftOver.bedGraph"), protein_file)
    name_sh <- paste0(output_prefix, "_stereogene_internal.sh")
    doCall(runStereogene, args = list(...,
        track_files = track_files,
        name_sh = name_sh,
        name_config = name_config))
    track_files <- c(paste0(dir_CapR_bg, "/", input_prefix,
        "_exterior_liftOver.bedGraph"), protein_file)
    name_sh <- paste0(output_prefix, "_stereogene_exterior.sh")
    doCall(runStereogene, args = list(...,
        track_files = track_files,
        name_sh = name_sh,
        name_config = name_config))
    track_files <- c(paste0(dir_CapR_bg, "/", input_prefix,
        "_stem_liftOver.bedGraph"), protein_file)
    name_sh <- paste0(output_prefix, "_stereogene_stem.sh")
    doCall(runStereogene, args = list(...,
        track_files = track_files,
        name_sh = name_sh,
        name_config = name_config))
    track_files <- c(paste0(dir_CapR_bg, "/", input_prefix,
        "_multibranch_liftOver.bedGraph"), protein_file)
    name_sh <- paste0(output_prefix, "_stereogene_multibranch.sh")
    doCall(runStereogene, args = list(...,
        track_files = track_files,
        name_sh = name_sh,
        name_config = name_config))
    track_files <- c(paste0(dir_CapR_bg, "/", input_prefix,
        "_bulge_liftOver.bedGraph"), protein_file)
    name_sh <- paste0(output_prefix, "_stereogene_bulge.sh")
    doCall(runStereogene, args = list(...,
        track_files = track_files,
        name_sh = name_sh,
        name_config = name_config))
}
