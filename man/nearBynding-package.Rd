\name{nearBynding}
\alias{nearBynding}
\alias{nearBynding-package}
\docType{package}
\title{Discern RNA structure proximal to protein binding}

\description{
nearBynding is a package designed to discern annotated RNA structures proximal
to protein binding. It allows users to annotate RNA structure contexts
via CapR or input their own annotations in BED/bedGraph format and it
accomodates protein binding information from CLIP-seq experiments as either
aligned CLIP-seq reads or peak-called intervals.
}

\details{
\tabular{ll}{
Package: \tab nearBynding\cr
Type: \tab Package\cr
Title: \tab nearBynding package \cr
Version: \tab 0.99.4\cr
Date: \tab July 21, 2020`\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
URL: \tab http://github.com/vbusa1/nearBynding\cr
}}

\author{Veronica Busa \email{vbusa1@jhmi.edu}}

\references{
\tabular{ll}{
StereoGene: \tab Stavrovskaya, Elena D., Tejasvi Niranjan, Elana J. Fertig,
Sarah J. Wheelan, Alexander V. Favorov, and Andrey A. Mironov. “StereoGene:
Rapid Estimation of Genome-Wide Correlation of Continuous or Interval Feature
Data.” Bioinformatics 33, no. 20 (October 15, 2017): 3158–65.
https://doi.org/10.1093/bioinformatics/btx379.\cr
CapR: \tab Tsukasa Fukunaga, Haruka Ozaki, Goro Terai, Kiyoshi Asai, Wataru
Iwasaki, and Hisanori Kiryu. "CapR: revealing structural specificities of
RNA-binding protein target recognition using CLIP-seq data." Genome Biology,
15, R16. (2014)
}}

\examples{
\dontrun{

library(nearBynding)
library(Rsamtools)

# get transcript list
load(system.file("extdata/transcript_list.Rda", package="nearBynding"))
# get GTF file
gtf<-system.file("extdata/Homo_sapiens.GRCh38.chr4&5.gtf",
                 package="nearBynding")
# make chain file
GenomeMappingToChainFile(genome_gtf = gtf,
                         out_chain_name = "test.chain",
                         RNA_fragment = "three_prime_utr",
                         transcript_list = transcript_list,
                         alignment = "hg38")
# get size of chromosomes of chain file
getChainChrSize(chain = "test.chain",
                out_chr = "chr4and5_3UTR.size")

# get transcript sequences
ExtractTranscriptomeSequence(transcript_list = transcript_list,
                     ref_genome = "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                     genome_gtf = gtf,
                     RNA_fragment = "three_prime_utr",
                     exome_prefix = "chr4and5_3UTR")
# run CapR on extracted sequences
runCapR(in_file = "chr4and5_3UTR.fa")

# get BAM file of protein binding
bam <- system.file("extdata/chr4and5.bam", package="nearBynding")
# sort it and convert to bedGraph format
sorted_bam<-sortBam(bam, "chr4and5_sorted")
CleanBAMtoBG(in_bam = sorted_bam)

# lift over protein binding and RNA structure to chain
liftOverToExomicBG(input = "chr4and5_sorted.bedGraph",
                   chain = "test.chain",
                   chrom_size = "chr4and5_3UTR.size",
                   output_bg = "chr4and5_liftOver.bedGraph")
processCapRout(CapR_outfile = "chr4and5_3UTR.out",
               chain = "test.chain",
               output_prefix = "chr4and5_3UTR",
               chrom_size = "chr4and5_3UTR.size",
               genome_gtf = gtf,
               RNA_fragment = "three_prime_utr")

# input to StereoGene
runStereogeneOnCapR(protein_file = "chr4and5_liftOver.bedGraph",
                    chrom_size = "chr4and5_3UTR.size",
                    name_config = "chr4and5_3UTR.cfg",
                    input_prefix = "chr4and5_3UTR")

# visualize protein binding context
visualizeCapRStereogene(CapR_prefix = "chr4and5_3UTR",
                        protein_file = "chr4and5_liftOver",
                        heatmap = T,
                        out_file = "all_contexts_heatmap",
                        x_lim = c(-500, 500))
}}

\keyword{ package }
\seealso{ See the nearBynding package vignette.}
