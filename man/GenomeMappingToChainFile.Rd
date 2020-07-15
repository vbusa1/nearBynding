\name{GenomeMappingToChainFile}
\alias{GenomeMappingToChainFile}
\title{
Create chain file from genome
}
\description{
Writes a chain file of non-overlapping isoforms mapped from a genome annotation
(GTF/GFF) file.}
\usage{
GenomeMappingToChainFile(genome_gtf,
                        out_chain_name,
                        RNA_fragment="exon",
                        transcript_list,
                        chrom_suffix = "exome",
                        verbose=FALSE,
                        alignment="hg19",
                        check_overwrite = FALSE)
}
\arguments{
  \item{genome_gtf}{The name of the GTF/GFF file that will be converted to an
  exome chain file. Required}
  \item{out_chain_name}{The name of the chain file to be created. Required}
  \item{RNA_fragment}{RNA component of interest. Options depend on the gtf file
  but often include "gene", "transcript", "exon", "CDS", "five_prime_utr",
  and/or "three_prime_utr". Default "exon" for the whole exome.}
  \item{transcript_list}{A vector of transcript names that represent the most
  expressed isoform of their respective genes. Transcript names must correspond
  to gtf annotation names. Isoforms cannot overlap. Required}
  \item{chrom_suffix}{The suffix to be appended to all chromosome names created
  in the chain file. Default "exome"}
  \item{verbose}{Output updates while the function is running. Default FALSE}
  \item{alignment}{The human genome alignment used, either "hg19" or "hg38".
  Suggested naming format: 'chr#' or just # for chromosome name, e.g. "chr1"" or
  "1". Default alignment is "hg19"}
  \item{check_overwrite}{Check for file wth the same out_chain_name before
  writing new chain file. Default FALSE.}
}
\values{
writes a chain file into directory
}
\author{Veronica Busa}

