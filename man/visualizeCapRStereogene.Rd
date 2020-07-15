\name{visualizeCapRStereogene}
\alias{visualizeCapRStereogene}
\title{
Visualize binding for CapR contexts}
\description{
Creates a visual output of all CapR RNA structure contexts relative to protein
binding.}
\usage{
visualizeCapRStereogene(dir_stereogene_output = ".",
                                  CapR_prefix,
                                  protein_file,
                                  protein_file_input=NULL,
                                  x_lim = c(-100, 100),
                                  y_lim = NULL,
                                  out_file = "out_file",
                                  legend = TRUE,
                                  heatmap = FALSE)
}
\arguments{
  \item{dir_stereogene_output}{Directory of stereogene output. Default working
  directory}
  \item{CapR_prefix}{The prefix common to CapR output files of protein_file.
  Required.}
  \item{protein_file}{A vector of at least one protein file name to be averaged
  for visualization. File names must exclude extensions such as ".bedGraph".
  All files in the list should be experimental or biological replicates.
  Required.}
  \item{protein_file_input}{A protein file name of background input to be
  subtracted from protein_file signal. File name must exclude extension. Only
  one input file is permitted. Optional.}
  \item{x_lim}{A vector of two values denoting the lower and upper x axis
  limits. Cannot exceed wSize/2 from write_config. Default (-100, 100)}
  \item{y_lim}{A vector of two values denoting the lower and upper y axis
  limits. Optional}
  \item{out_file}{Name of output file, excluding extension. ".pdf" or ".jpeg"
  will be added as relevant to the output file name. Default "out_file"}
  \item{legend}{Whether a legend should be included with the output graph.
  Default TRUE}
  \item{heatmap}{Whether the output graph should be in the form of a heatmap
  (TRUE) or of a line graph (FALSE). Default FALSE}
}
\values{
heatmap (JPEG) or line graph (PDF) image file
}
\examples{
## pull example files
get_outfiles()
## heatmap
visualizeCapRStereogene(CapR_prefix = "chr4and5_3UTR",
                        protein_file = "chr4and5_liftOver",
                        heatmap = T,
                        out_file = "all_contexts_heatmap",
                        x_lim = c(-500, 500))
## line graph
visualizeCapRStereogene(CapR_prefix = "chr4and5_3UTR",
                        protein_file = "chr4and5_liftOver",
                        x_lim = c(-500, 500),
                        out_file = "all_contexts_line",
                        y_lim = c(-18, 22))
}
\author{Veronica Busa}
