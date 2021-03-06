% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processCapRout.R
\name{processCapRout}
\alias{processCapRout}
\title{processCapRout}
\usage{
processCapRout(
  CapR_outfile,
  output_prefix,
  chrom_size,
  genome_gtf,
  RNA_fragment,
  chain
)
}
\arguments{
\item{CapR_outfile}{Name of CapR output file. Required}

\item{output_prefix}{Prefix string to be appended to all output files.
Required.}

\item{chrom_size}{Name of chromosome size file. File must be in two-column
format without a header where first column is chromosome name and second
column is chromosome length, as from getChainChrSize. Required.}

\item{genome_gtf}{The name of the GTF/GFF file that contains all exome
annotations. Required}

\item{RNA_fragment}{RNA component of interest. Options depend on the gtf
file but often include "gene", "transcript", "exon", "CDS", "five_prime_utr",
and/or "three_prime_utr". Default "exon" for the whole exome.}

\item{chain}{The name of the chain file to be used. Format should be like
chain files derived from GRangesMappingToChainFile. Required}
}
\value{
writes bedGraph files of structure signal for each of the six
CapR contexts 1) mapped to the genome and 2) lifted-over to the transcriptome
}
\description{
Creates context-separated bedGraph files of CapR output for
genome and transcriptome alignments.
}
\examples{
## make chain file
load(system.file("extdata/transcript_list.Rda", package="nearBynding"))
gtf<-system.file("extdata/Homo_sapiens.GRCh38.chr4&5.gtf",
                package="nearBynding")
GenomeMappingToChainFile(genome_gtf = gtf,
                        out_chain_name = "test.chain",
                        RNA_fragment = "three_prime_utr",
                        transcript_list = transcript_list,
                        alignment = "hg38")
## get chromosome size file
getChainChrSize(chain = "test.chain",
               out_chr = "chr4and5_3UTR.size")

processCapRout(CapR_outfile = system.file("extdata/chr4and5_3UTR.out",
                                         package="nearBynding"),
              chain = "test.chain",
              output_prefix = "chr4and5_3UTR",
              chrom_size = "chr4and5_3UTR.size",
              genome_gtf = gtf,
              RNA_fragment = "three_prime_utr")

}
