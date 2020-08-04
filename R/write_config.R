#' @title write_config
#'
#' @description Writes a configuration file for use by Stereogenes in the
#' working directory.
#'
#' @param name_config Name of output config file. Default config.cfg
#' @param chrom_size Name of chromosome size file in two-column format without a
#' header where first column is chromosome name and second column is chromosome
#' length, as from GenomeMappingToChainFile. Required
#' @param Rscript Write R script for the result presentation. Equivalent to -r
#' argument in StereoGene. Default TRUE
#' @param verbose Provides a verbose output when Stereogene is run. Equivalent
#' to -v or -verbose argument in StereoGene.Default FALSE
#' @param na_noise Use NA values as unknown and fill them with noise. Equivalent
#' to -NA argument in StereoGene. Default FALSE
#' @param bin Bin size (integer) for input averaging. Default 1
#' @param threshold Threshold for input data to remove small values. An integer
#' between 0 and 250. Default 0
#' @param cross_width Width of cross-correlation plot output in Rscript.
#' Default 200.
#' @param wSize Window size. If windows are too small, cross correlations will
#' have a lot of noise; if they are too large, there may be too few windows for
#' robust statistical assessment. Default 10000
#' @param kernel_width Kernel span in nucleotides. Equivalent to KernelSigma in
#' StereoGene. Default 1000
#' @param outLC Write local kerneled correlations into a bedgraph file.
#' Default FALSE.
#' @param LCScale Local correlation scale: logarithmic ("LOG") or linear ("LIN")
#' scaling. Default "LOG".
#' @param LC_FDR Threshold for local kernel correlation FDR to be written into
#' the local correlation file. Default 0.5
#'
#' @note Not all StereoGene parameters are included in this function so refer to
#' the StereoGene manual and modify the output .cfg file manually if additional
#' parameters are desired.
#'
#' @export

write_config <- function(name_config = "config.cfg",
                            chrom_size,
                            Rscript = TRUE,
                            verbose = FALSE,
                            na_noise = FALSE,
                            bin = 1,
                            threshold = 0,
                            cross_width = 200,
                            wSize = 10000,
                            kernel_width = 1000,
                            outLC = FALSE,
                            LCScale = "LOG",
                            LC_FDR = .5) {
    if (missing(chrom_size)) {stop("please provide a chrom_size file")}
    if (Rscript == TRUE) {R <- 1} else {R <- 0}
    if (verbose == TRUE) {V <- 1} else {V <- 0}
    if (na_noise == TRUE) {N <- 1} else {N <- 0}
    if (outLC == TRUE) {L <- 1} else {L <- 0}
    conf <- paste("#!bash/bin", "",
        "profPath =./",
        "trackPath=./",
        "resPath=./",
        paste0("chrom=", chrom_size),
        paste0("Rscript=", R),
        paste0("verbose=", V),
        paste0("NA=", N),
        "outRes=TAB",
        "writeDistr=NONE",
        "Distances=0",
        paste0("bin=", bin),
        paste0("threshold=", threshold),
        paste0("wSize=", wSize),
        "maxZero=100",
        "maxNA=100",
        "nShuffle=100000",
        paste0("outLC=", L),
        paste0("LCScale=", LCScale),
        paste0("crossWidth=", cross_width),
        paste0("KernelSigma=", kernel_width),
        paste0("L_FDR =", LC_FDR),
        paste0("R_FDR =", LC_FDR),
        sep = "\n"
    )
    writeLines(conf, name_config)
}
