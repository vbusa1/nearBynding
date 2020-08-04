\name{CleanBEDtoBG}
\alias{CleanBEDtoBG}
\title{
Write BED to bedGraph}
\description{
Writes a script to convert a BED file to a clean bedGraph file.}
\usage{
CleanBEDtoBG(in_bed,
            out_bedGraph = NA,
            unwanted_chromosomes = NULL,
            alignment = "hg19")
}

\arguments{
    \item{in_bed}{Name of sorted BAM file to be converted to a BedGraph file.
    Required}
    \item{out_bedGraph}{Name of bedGraph output file, including full directory
    path. Default in_bam prefix and directory.}
    \item{unwanted_chromosomes}{A vector of unwanted chromosomes that are
    present in the BAM file.}
    \item{alignment}{The human genome alignment used, either "hg19" or "hg38".
    Default "hg19"}
}
\value{
deposits bedGraph from BED in same directory
}
\examples{
bam <- system.file("extdata/chr4and5.bam", package="nearBynding")
out_bed <- "bamto.bed"
\dontrun{
system2("bedtools", paste0("bamtobed -i ", bam, " > ", out_bed))

CleanBEDtoBG(in_bed = out_bed,
            alignment = "hg38")
}
}
