% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bindingContextDistanceCapR.R
\name{bindingContextDistanceCapR}
\alias{bindingContextDistanceCapR}
\title{bindingContextDistanceCapR}
\usage{
bindingContextDistanceCapR(
  dir_stereogene_output = ".",
  CapR_prefix = "",
  protein_file,
  protein_file_input = NULL,
  dir_stereogene_output_2 = NULL,
  CapR_prefix_2 = "",
  protein_file_2,
  protein_file_input_2 = NULL,
  context = "all",
  range = c(-200, 200)
)
}
\arguments{
\item{dir_stereogene_output}{Directory of Stereogene output for first
protein. Default current directory.}

\item{CapR_prefix}{The prefix common to CapR output files of protein_file, if
applicable. Equivalent to output_prefix from runStereogeneOnCapR. Default ""}

\item{protein_file}{A vector of strings with at least one protein file name
to be averaged for calculation of distance. File names must exclude
extensions such as ".bedGraph". All files in the list should be
experimental/biological replicates. Required.}

\item{protein_file_input}{A protein file name of background input to be
subtracted from protein_file signal. File name must exclude extension. Only
one input file is permitted. Optional.}

\item{dir_stereogene_output_2}{Directory of Stereogene output for second
protein. Default current directory.}

\item{CapR_prefix_2}{The prefix common to CapR output files of
protein_file_2, if applicable.Equivalent to output_prefix from r
unStereogeneOnCapR. Default ""}

\item{protein_file_2}{Similar to protein_file. A second vector of at least
one protein file name to be averaged for calculation of distance. File names
must exclude extensions such as ".bedGraph". All files in the list should be
experimental/biological replicates. Required.}

\item{protein_file_input_2}{Similar to protein_file_input. A second protein
file name of background input to be subtracted from protein_file_2 signal.
File name must exclude extension. Only one input file is permitted. Optional.}

\item{context}{The RNA structure context being compared for the two protein
file sets. Acceptable contexts include "all", which sums the distance of all
six contexts, or any of the contexts individually ("bulge", "hairpin",
"stem", "exterior", "multibranch", or "internal"). Default "all"}

\item{range}{A vector of two integers denoting the range upstream and
downstream of the center of protein binding to consider in the comparison.
Ranges that are too small miss the holistic binding context, while large
ranges amplify distal noise in the binding data. Cannot exceed wSize/2 from
write_config. Default c(-200, 200)}
}
\value{
Wasserstein distance between the two protein file sets provided for
the RNA structure context specified, minus the input binding signal if
applicable
}
\description{
Calculate the Wasserstein distance between two replicates' or
two proteins' binding contexts.
}
\note{
Wasserstein distance calculations are reciprocal, so it does not matter
which protein is first or second so long as replicates and input files
correspond to one another.
}
\examples{
## load example StereoGene output
get_outfiles()

## This boring example compares a protein's binding with itself for all
## contexts, therefore the distance is 0
bindingContextDistanceCapR(CapR_prefix = "chr4and5_3UTR",
                           protein_file = "chr4and5_liftOver",
                           CapR_prefix_2 = "chr4and5_3UTR",
                           protein_file_2 = "chr4and5_liftOver")

}
