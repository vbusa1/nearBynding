\name{bindingContextDistance}
\alias{bindingContextDistance}
\title{
Distance between binding contexts}
\description{
Calculate the Wasserstein distance between two replicates' or two proteins'
binding contexts}
\usage{
bindingContextDistance(dir_stereogene_output = ".",
                                 RNA_context,
                                 protein_file,
                                 protein_file_input = NULL,
                                 dir_stereogene_output_2 = NULL,
                                 RNA_context_2 = NULL,
                                 protein_file_2 = NULL,
                                 protein_file_input_2 = NULL,
                                 range = c(-200, 200))
}
\arguments{
  \item{dir_stereogene_output}{Directory of Stereogene output for first protein.
  Default current directory.}
  \item{RNA_context}{Name of the RNA context file input to Stereogene.File names
  must exclude extensions such as ".bedGraph". Requred}
  \item{protein_file}{A vector of at least one protein file name to be averaged
  for calculation of distance. File names must exclude extensions such as
  ".bedGraph". All files in the list should be experimental/biological
  replicates. Required.}
  \item{protein_file_input}{A protein file name of background input to be
  subtracted from protein_file signal. File name must exclude extension.
  Only one input file is permitted. Optional.}
  \item{dir_stereogene_output_2}{Directory of Stereogene output for second
  protein. Default dir_stereogene_output.}
  \item{RNA_context_2}{Name of the RNA context file input to Stereogene.File
  names must exclude extensions such as ".bedGraph". Default same as
  RNA_context.}
  \item{protein_file_2}{Similar to protein_file. A second vector of at least one
  protein file name to be averaged for calculation of distance. File names must
  exclude extensions such as ".bedGraph". All files in the list should be
  experimental/biological replicates. Default same as protein_file}
  \item{protein_file_input_2}{Similar to protein_file_input. A second protein
  file name of background input to be subtracted from protein_file_2 signal.
  File name must exclude extension. Only one input file is permitted. Optional.}
  \item{range}{The range upstream and downstream of the center of protein
  binding to consider in the comparison. Ranges that are too small miss the
  holistic binding context, while large ranges amplify distal noise in the
  binding data. Cannot exceed wSize/2 from write_config. Default c(-200, 200)}
}
