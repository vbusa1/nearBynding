\name{liftOverToExomicBG}
\alias{liftOverToExomicBG}
\title{
Lift features to transcriptome}
\description{
Lifts features such as CLIP-seq reads or RNA structure annotations from genome
to transcriptome.}
\usage{
liftOverToExomicBG(input,
                     chain,
                     output_bg,
                     write_chr = TRUE,
                     out_chr,
                     format = "bedGraph")
}
\arguments{
  \item{input}{A single input file name or a vector of input file names in the
  format of c(forward_reads, reverse_reads) for strand-separated alignments.
  Files must be BED or bedGraph format. Required}
  \item{chain}{The name of the chain file to be used for liftOver. Format should
  be like chain files derived from GRangesMappingToChainFile function. Required}
  \item{output_bg}{The name of the lifted-over output bedGraph file. Required}
  \item{write_chr}{Whether to output a table of mapped chromosome names and
  lengths. Recommended TRUE for downstream pipeline functions. Default TRUE}
  \item{out_chr}{Name of the chromosome names and lengths table file. Required
  if write_chr is TRUE}
  \item{format}{File type of input file(s). Recommended "BED" or "bedGraph".
  Default "bedGraph"}
}
\values{
writes lifted-over bedGraph file
}
\author{Veronica Busa}
