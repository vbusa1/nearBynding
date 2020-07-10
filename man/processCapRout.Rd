\name{processCapRout}
\alias{processCapRout}
\title{
CapR outfile to bedGraphs}
\description{
Creates context-separated bedGraph files of CapR output for genome and
transcriptome alignments.}
\usage{
processCapRout(CapR_outfile,
               output_prefix,
               chrom_size,
               genome_gtf,
               RNA_fragment,
               chain)
}
\arguments{
  \item{CapR_outfile}{Name of CapR output file. Required}
  \item{output_prefix}{Prefix to be appended to all output files. Required}
  \item{chrom_size}{Name of chromosome size file in two-column format without a
  header where first column is chromosome name and second column is chromosome
  length, as from liftOverToExomicBG. Required}
  \item{genome_gtf}{The name of the GTF/GFF file that contains all exome
  annotations. Required}
  \item{RNA_fragment}{RNA component of interest. Options depend on the gtf file
  but often include "gene", "transcript", "exon", "CDS", "five_prime_utr",
  and/or "three_prime_utr". Default "exon" for the whole exome.}
  \item{chain}{The name of the chain file to be used. Format should be like
  chain files derived from GRangesMappingToChainFile. Required}
}
