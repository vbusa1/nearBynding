#' @title ExtractTranscriptomeSequence
#'
#' @description Writes a FASTA file of transcript sequences from a list of
#' transcripts.
#'
#' @param transcript_list a vector of transcript names that represent the most
#' expressed isoform of their respective genes and correspond to GTF annotation
#' names. Required
#' @param ref_genome The name of the reference genome FASTA from which exome
#' sequences will be derived. Required
#' @param genome_gtf The name of the GTF/GFF file that contains all exome
#' annotations. Coordinates must match the file input for the ref_genome
#' parameter. Required
#' @param RNA_fragment RNA component of interest. Options depend on the gtf
#' file but often include "gene", "transcript", "exon", "CDS", "five_prime_utr",
#' and/or "three_prime_utr". Default "exon" for the whole exome.
#' @param exome_prefix The prefix for all output files. Default "exome"
#'
#' @note transcript_list, genome_gtf, and RNA_fragment arguments should be the
#' same as GenomeMappingToChainFile function arguments
#'
#' @export

ExtractTranscriptomeSequence <- function(transcript_list,
                                         ref_genome,
                                         genome_gtf,
                                         RNA_fragment = "exon",
                                         exome_prefix = "exome") {
  # make sure there is a gtf file available; if there is, load in the
  # information; if not, request one
  if (missing(genome_gtf)) {
    return("Please provide a genome GTF file via the genome_gtf parameter.")
  } else {
    gtf <- import(genome_gtf)
    gtf <- with(gtf, plyranges::filter(gtf, type == RNA_fragment))
    gtf_transcripts <- gtf [(elementMetadata(gtf)[, "transcript_id"] %in%
      transcript_list)]
    input_bed <- paste0(exome_prefix, ".bed")
    df <- data.frame(
      seqnames = seqnames(gtf_transcripts),
      starts = start(gtf_transcripts) - 1,
      ends = end(gtf_transcripts),
      names = elementMetadata(gtf_transcripts)[, "transcript_id"],
      strands = strand(gtf_transcripts)
    )
    # remove junky transcripts and write BED file
    df <- filter(df, nchar(as.character(seqnames)) < 6)
    write.table(df, input_bed,
      sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = FALSE
    )
  }

  # call bedtools on local system, get FASTA sequences
  system2("bedtools", paste0(
    "getfasta -s -fi ",
    ref_genome, " -bed ", input_bed, " -name -tab -fo ",
    exome_prefix, ".tmp"
  ))

  # FASTA sequences are fragmented;
  # must compile sequences into whole transcripts to input into
  # CapR folding algorithm
  edit_tbl <- read.table(paste0(exome_prefix, ".tmp"),
    stringsAsFactors = FALSE
  )
  system2("rm", paste0(exome_prefix, ".tmp"))
  edit_tbl$V1 <- substr(edit_tbl$V1, 1, nchar(edit_tbl$V1) - 2)
  edit_tbl$strand <- df$strand[match(df$names, edit_tbl$V1)] # get strand info

  # make matrix for name and sequences to be written into FASTA
  hold_matrix <- matrix(NA, ncol = 2, nrow = length(unique(edit_tbl$V1)))
  hold_matrix[, 1] <- unique(edit_tbl$V1)
  for (txpt in unique(edit_tbl$V1)) {
    filter(edit_tbl, .data$V1 == txpt) -> txpt_tbl
    if (unique(txpt_tbl$strand) == "-") { # write - strand seqs in reverse
      seq <- paste(rev(txpt_tbl$V2), collapse = "", sep = "")
    } else {
      seq <- paste(txpt_tbl$V2, collapse = "", sep = "")
    }
    hold_matrix[which(hold_matrix[, 1] == txpt), 2] <- seq
  }

  write_fasta(hold_matrix[, 2], hold_matrix[, 1], paste0(exome_prefix, ".fa"))
}
