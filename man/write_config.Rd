\name{write_config}
\alias{write_config}
\title{
Write configuration file
}
\description{
Writes a configuration file for use by Stereogenes in the working directory}
\usage{
write_config(name_config = "config.cfg",
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
            LC_FDR = .5)
}
\arguments{
    \item{name_config}{Name of output config file. Default "config.cfg"}
    \item{chrom_size}{Name of chromosome size file in two-column format without
    a header where first column is chromosome name and second column is
    chromosome length, as from GenomeMappingToChainFile. Required}
    \item{Rscript}{Rscript Write R script for the result presentation.
    Equivalent to -r argument in StereoGene. Default TRUE}
    \item{verbose}{Provides a verbose output when Stereogene is run. Equivalent
    to -v or -verbose argument in StereoGene.Default FALSE}
    \item{na_noise}{Use NA values as unknown and fill them with noise.
    Equivalent to -NA argument in StereoGene. Default FALSE}
    \item{bin}{Bin size (integer) for input averaging. Larger values provide
    faster run time but less precise protein binding information. Default 1}
    \item{threshold}{Threshold for input data to remove small values. An integer
    between 0 and 250. Default 0}
    \item{cross_width}{Width of cross-correlation plot output in Rscript.
    Irrelevant if using visualizatiion through other functions in this package.
    Default 200.}
    \item{wSize}{Window size. If windows are too small, cross correlations will
    have a lot of noise; if they are too large, there may be too few windows for
    robust statistical assessment. Default 10000}
    \item{kernel_width}{Kernel span in nucleotides. Equivalent to KernelSigma in
    StereoGene. Must be smaller than wSize. Default 1000}
    \item{outLC}{Write local kerneled correlations into a bedgraph file.
    Default FALSE.}
    \item{LCScale}{Local correlation scale: logarithmic ("LOG") or linear
    ("LIN") scaling. Default "LOG".}
    \item{LC_FDR}{Threshold for local kernel correlation FDR to be written into
    the local correlation file. Default 0.5}
}
\note{Not all StereoGene parameters are included in this function so refer to
the StereoGene manual and modify the output .cfg file manually if additional
parameters are desired.}
\value{
writes a configuration file into directory
}
\examples{
## Write a config file named "test.cfg" with chromosome size file "test.size"
write_config(name_config = "test.cfg",
                chrom_size = "test.size")
}
\references{
Stavrovskaya, Elena D., Tejasvi Niranjan, Elana J. Fertig, Sarah J. Wheelan,
Alexander V. Favorov, and Andrey A. Mironov. “StereoGene: Rapid Estimation
of Genome-Wide Correlation of Continuous or Interval Feature Data.”
Bioinformatics 33, no. 20 (October 15, 2017): 3158–65.
https://doi.org/10.1093/bioinformatics/btx379.
}
\author{Veronica Busa}

