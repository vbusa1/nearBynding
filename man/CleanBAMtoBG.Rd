\name{CleanBAMtoBG}
\alias{CleanBAMtoBG}
\title{
Write BAM to bedGraph}
\description{
Writes a script to convert a sorted BAM file to a clean bedGraph file.}
\usage{
CleanBAMtoBG(in_bam,
               out_bedGraph = NA,
               unwanted_chromosomes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{in_bam}{Name of sorted BAM file to be converted to a BedGraph file.
  Required}
  \item{out_bedGraph}{Name of bedGraph output file, including full directory
  path. Default in_bam prefix and directory.}
  \item{unwanted_chromosomes}{A vector of unwanted chromosomes that are present
  in the BAM file.}
}

\author{Veronica Busa}


\examples{
\dontrun{
CleanBAMtoBG(in_bam="ENCFF288LEG.bam",
unwanted_chromosomes="EBV")
}
}
