#' @title write_config
#'
#' @description Writes a configuration file for use by Stereogenes in the
#' working directory.
#'
#' @param name_config Name of output config file. Default config.cfg
#' @param chrom_size Name of chromosome size file. File must be in two-column
#' format without a header where first column is chromosome name and second
#' column is chromosome length, as from getChainChrSize. Required
#' @param Rscript Write R script for the result presentation. Equivalent to -r
#' argument in StereoGene. Default FALSE
#' @param silent Provides an output when Stereogene is run. Equivalent
#' to -s or -silent argument in StereoGene.Default TRUUE
#' @param na_noise Use NA values as unknown and fill them with noise. Equivalent
#' to -NA argument in StereoGene. Default FALSE
#' @param bin Bin size for input averaging; an integer. Default 1
#' @param threshold Threshold for input data to remove small values. An integer
#' between 0 and 250. Default 0
#' @param cross_width Width of cross-correlation plot output in Rscript; an
#' integer. Default 200.
#' @param wSize Window size; an integer. If windows are too small, cross
#' correlations will have a lot of noise; if they are too large, there may be
#' too few windows for robust statistical assessment. Default 10000
#' @param kernel_width Kernel span in nucleotides; an integer. Equivalent to
#' KernelSigma invStereoGene. Default 1000
#' @param resPath Folder to store results. Default is current directory.
#'
#' @return writes a configuration file into directory
#'
#' @note Not all StereoGene parameters are included in this function so refer to
#' the StereoGene manual and modify the output .cfg file manually if additional
#' parameters are desired.
#'
#' @examples
#' ## Write a config file named "test.cfg" with chromosome size file "test.size"
#' write_config(name_config = "test.cfg",
#'             chrom_size = "test.size")
#'
#' @export

write_config <- function(name_config = "config.cfg",
                            chrom_size,
                            Rscript = FALSE,
                            silent = TRUE,
                            na_noise = FALSE,
                            bin = 1,
                            threshold = 0,
                            cross_width = 200,
                            wSize = 10000,
                            kernel_width = 1000,
                            resPath = ".") {
    if (missing(chrom_size)) {stop("please provide a chrom_size file")}
    if (Rscript == TRUE) {R <- 1} else {R <- 0}
    if (silent == TRUE) {S <- 1} else {S <- 0}
    if (na_noise == TRUE) {N <- 1} else {N <- 0}
    conf <- paste("#!bash/bin", "",
        "profPath =./",
        "trackPath=./",
        paste0("resPath=", resPath, "/"),
        paste0("chrom=", chrom_size),
        paste0("Rscript=", R),
        paste0("silent=", S),
        paste0("NA=", N),
        "outRes=TAB",
        "writeDistr=NONE",
        "Distances=1",
        paste0("bin=", bin),
        paste0("threshold=", threshold),
        paste0("wSize=", format(wSize, scientific = F)),
        "maxZero=100",
        "maxNA=100",
        "outLC=0",
        paste0("crossWidth=", cross_width),
        paste0("KernelSigma=", kernel_width),
        sep = "\n"
    )
    writeLines(conf, name_config)
}
