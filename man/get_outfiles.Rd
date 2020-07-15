\name{get_outfiles}
\alias{get_outfiles}
\title{
Files for vignette}
\description{
Copy files necessary to complete the vignette onto the local machine in cases
where Stereogene, CapR, or bedtools are not available.}
\usage{
get_outfiles(dir = ".")
}
\arguments{
  \item{dir}{Directory into which files ought to be stored. Default current
  work directory.}
}
\value{
deposits six *.dist StereoGene output files into the selected directory
}
\examples{
## pull example StereoGene output files
get_outfiles()
}
