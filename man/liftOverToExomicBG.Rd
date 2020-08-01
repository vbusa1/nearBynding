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
                     chrom_size,
                     output_bg,
                     format = "bedGraph")
}
\arguments{
  \item{input}{A single input file name or a vector of input file names in the
  format of c(forward_reads, reverse_reads) for strand-separated alignments.
  Files must be BED or bedGraph format. Required}
  \item{chain}{The name of the chain file to be used for liftOver. Format should
  be like chain files derived from GRangesMappingToChainFile function. Required}
  \item{chrom_size}{Name of chromosome size file in two-column format without a
  header where first column is chromosome name and second column is chromosome
  length, as from getChainChrSize. Required}
  \item{output_bg}{The name of the lifted-over output bedGraph file. Required}
  \item{format}{File type of input file(s). Recommended "BED" or "bedGraph".
  Default "bedGraph"}
}
\value{
writes lifted-over bedGraph file
}
\examples{
## first, get chain file
load(system.file("extdata/transcript_list.Rda", package="nearBynding"))
gtf<-system.file("extdata/Homo_sapiens.GRCh38.chr4&5.gtf",
                 package="nearBynding")
GenomeMappingToChainFile(genome_gtf = gtf,
                         out_chain_name = "test.chain",
                         RNA_fragment = "three_prime_utr",
                         transcript_list = transcript_list,
                         alignment = "hg38")
## and chain file chromosome sizes
getChainChrSize(chain = "test.chain",
                out_chr = "chr4and5_3UTR.size")

## get bedGraph file
chr4and5_sorted.bedGraph<-system.file("extdata/chr4and5_sorted.bedGraph",
                                        package="nearBynding")

liftOverToExomicBG(input = chr4and5_sorted.bedGraph,
                   chain = "test.chain",
                   chrom_size = "chr4and5_3UTR.size",
                   output_bg = "chr4and5_liftOver.bedGraph")
}
\author{Veronica Busa}
