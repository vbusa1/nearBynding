\name{runStereogene}
\alias{runStereogene}
\title{
Write file to run Stereogene
}
\description{
Writes a file to the working directory to run StereoGene}
\usage{
runStereogene(track_files,
                name_config,
                pcorProfile = NULL,
                confounder = NULL)
}
\arguments{
    \item{track_files}{Vector of at least two track or interval file names to be
    pairwise-correlated by StereoGene. Required.}
    \item{name_config}{Name of corresponding configuration file. Required}
    \item{pcorProfile}{Track for partial correlation. More information for this
    can be found in the StereoGene README. Optional}
    \item{confounder}{Confounder filename. More information for this can be
    found in the StereoGene README. Optional}
}
\value{
generates StereoGene output files in directory
}

\examples{
\dontrun{
runStereogene(track_files = c("chr4and5_3UTR_stem_liftOver.bedGraph",
                                "chr4and5_liftOver.bedGraph"),
                name_config = "chr4and5_3UTR.cfg")
}
}
\references{
Stavrovskaya, Elena D., Tejasvi Niranjan, Elana J. Fertig, Sarah J. Wheelan,
Alexander V. Favorov, and Andrey A. Mironov. “StereoGene: Rapid Estimation
of Genome-Wide Correlation of Continuous or Interval Feature Data.”
Bioinformatics 33, no. 20 (October 15, 2017): 3158–65.
https://doi.org/10.1093/bioinformatics/btx379.
}
\author{Veronica Busa}

