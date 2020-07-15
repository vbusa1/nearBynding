\name{ExtractTranscriptomeSequence}
\alias{ExtractTranscriptomeSequence}
\title{
Write FASTA of transcriptome}
\description{
Writes a FASTA file of transcript sequences from a list of transcripts.}
\usage{
ExtractTranscriptomeSequence(transcript_list,
                               ref_genome,
                               genome_gtf,
                               RNA_fragment="exon",
                               exome_prefix="exome")
}
\arguments{
  \item{transcript_list}{a vector of transcript names that represent the most
  expressed isoform of their respective genes and correspond to GTF annotation
  names. Required}
  \item{ref_genome}{The name of the reference genome FASTA from which exome
  sequences will be derived. Required}
  \item{genome_gtf}{The name of the GTF/GFF file that contains all exome
  annotations. Coordinates must match the file input for the ref_genome
  parameter. Required}
  \item{RNA_fragment}{RNA component of interest. Options depend on the GTF
  file but often include "gene", "transcript", "exon", "CDS", "five_prime_utr",
  and/or "three_prime_utr". Default "exon" for the whole exome.}
  \item{exome_prefix}{The prefix for all output files. Default "exome"}
}
\value{
writes FASTA file of transcriptome sequences into directory
}
\author{Veronica Busa}
