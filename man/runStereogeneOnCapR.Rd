\name{runStereogeneOnCapR}
\alias{runStereogeneOnCapR}
\title{
Input CapR contexts to Stereogene}
\description{
Writes a configuration file and Stereogene script and runs Stereogene for all
CapR tracks}
\usage{
runStereogeneOnCapR(dir_CapR_bg = ".",
                              input_prefix,
                              protein_file,
                              output_prefix = input_prefix,
                              name_config = "config.cfg",
                              chrom_size,
                              ...)
}
\arguments{
  \item{dir_CapR_bg}{Location of lifted-over CapR bedGraph files. Default
  current directory}
  \item{input_prefix}{Prefix appended to input files; same as input_prefix
  argument in processCapRout. Required}
  \item{protein_file}{Name of protein file in bedGraph format. Required}
  \item{output_prefix}{Prefix to be appended to all output files. Default to be
  same as input_prefix}
  \item{name_config}{Name of output config file. Default "config.cfg"}
  \item{chrom_size}{Name of chromosome size file in two-column format without a
  header where first column is chromosome name and second column is chromosome
  length, as from liftOverToExomicBG. Required}
  \item{...}{includes all other parameters acceptable to write_config and
  write_stereogene}
}
\value{
generates StereoGene output files, including *.dist files
}
\references{
Stavrovskaya, Elena D., Tejasvi Niranjan, Elana J. Fertig, Sarah J. Wheelan,
Alexander V. Favorov, and Andrey A. Mironov. “StereoGene: Rapid Estimation
of Genome-Wide Correlation of Continuous or Interval Feature Data.”
Bioinformatics 33, no. 20 (October 15, 2017): 3158–65.
https://doi.org/10.1093/bioinformatics/btx379.
}
\author{Veronica Busa}
