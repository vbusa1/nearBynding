#' @title GenomeMappingToChainFile
#'
#' @description Writes a chain file mapped from a genome annotation file.
#'
#' @param genome_gtf The name of the GTF/GFF file that will be converted to an
#' exome chain file. Required
#' @param out_chain_name The name of the chain file to be created. Required
#' @param RNA_fragment RNA component of interest. Options depend on the gtf file
#' but often include "gene", "transcript", "exon", "CDS", "five_prime_utr",
#' and/or "three_prime_utr". Default "exon" for the whole exome.
#' @param transcript_list A vector of transcript names that represent the most
#' expressed isoform of their respective genes and correspond to gtf annotation
#' names. Isoforms cannot overlap. Required
#' @param chrom_suffix The suffix to be appended to all chromosome names created
#' in the chain file. Default "exome"
#' @param verbose Output updates while the function is running. Default FALSE
#' @param alignment The human genome alignment used, either "hg19" or "hg38".
#' Default "hg19"
#' @param check_overwrite Check for file wth the same out_chain_name
#' before writing new file. Default FALSE.
#'
#' @return writes a chain file into directory
#'
#' @importFrom S4Vectors elementMetadata Rle
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom methods as
#' @importFrom magrittr '%>%'
#' @importFrom utils write.table
#'
#' @examples
#' ## load transcript list
#'load(system.file("extdata/transcript_list.Rda", package="nearBynding"))
#'## get GTF file
#' gtf<-system.file("extdata/Homo_sapiens.GRCh38.chr4&5.gtf",
#'                 package="nearBynding")
#'
#' GenomeMappingToChainFile(genome_gtf = gtf,
#'                         out_chain_name = "test.chain",
#'                         RNA_fragment = "three_prime_utr",
#'                         transcript_list = transcript_list,
#'                         alignment = "hg38")
#'
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import rtracklayer
#' @import BiocGenerics
#'
#' @export

GenomeMappingToChainFile <- function(genome_gtf,
                                        out_chain_name,
                                        RNA_fragment = "exon",
                                        transcript_list,
                                        chrom_suffix = "exome",
                                        verbose = FALSE,
                                        alignment = "hg19",
                                        check_overwrite = FALSE) {
    if (alignment == "hg19") {
        seqinft <- as.data.frame(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))
    }
    if (alignment == "hg38") {
        seqinft <- as.data.frame(seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene))
    }
    if (!(alignment %in% c("hg19", "hg38"))) {
        stop("Acceptable alignments are 'hg19' and 'hg38'")
    }
    # create seqinfo object for all chromosomes to get lengths
    seqinf <- as(seqinft[nchar(rownames(seqinft)) < 6, ], "Seqinfo")
    if (check_overwrite == TRUE) {
        if (file.exists(out_chain_name)) {
            input <- readline(prompt = paste0(
                "file ", out_chain_name,
                " already exists. Do you want to replace it? y/n :"))
            if (input == "y" | input == "Y" | input == "yes") {
                system2("rm", out_chain_name)
            } else if (input == "n" | input == "N" | input == "no") {
                stop("Input a different out_chain_name to proceed.")
            } else {
                print("Please enter y or n")
                if (input == "y" | input == "Y" | input == "yes") {
                    system2("rm", out_chain_name)
                } else if (input == "n" | input == "N" | input == "no") {
                    stop("Input a different out_chain_name to proceed.")
                } else {
                    stop(paste0("Cannot proceed without indication of whether
                                        you wish to replace ", out_chain_name))
                }
            }
        }
    }
    gtf <- import(genome_gtf)
    gtf <- with(gtf, plyranges::filter(gtf, type == RNA_fragment))
    gtf_transcripts <- gtf[(elementMetadata(gtf)[,"transcript_id"] %in%
        transcript_list)]
    if (verbose == TRUE) {
        print("annotation data finished loading")
    }
    # write chain file chromosome-by-chromosome
    if (sum(seqnames(gtf_transcripts) %in% seqnames(seqinf)) == 0) {
        gtf_transcripts@seqnames <- paste0("chr", seqnames(gtf_transcripts)) %>%
            Rle()
        if (sum(seqnames(gtf_transcripts) %in% seqnames(seqinf)) == 0) {
            stop("GTF chromosomes do not resemble UCSC chromosome names.
                Suggested format: 'chr#' or just # for chromosome name,
                e.g. chr1 chr10 chrM")
        }
    }
    for (chr in seqnames(seqinf)) {
        for (str in c("-", "+")) {
            if (verbose == TRUE) {
                print(paste("Chromosome", chr,"strand", str, "starting"))
            }
            gtf_hold <- gtf_transcripts %>% plyranges::filter(
                strand == str,
                seqnames == chr)
            # in case of chromosomes without data on one strand
            if (length(ranges(gtf_hold)) == 0) {next}
            length_chr <- sum(width(ranges(gtf_hold)))
            if (str == "+") {sign <- "plus"} else {sign <- "minus"}
            transcripts <-
                elementMetadata(gtf_hold)@listData[["transcript_id"]] %>%
                unique()
            chrom_chains <- matrix(NA, ncol = 3, nrow = (length(gtf_hold) +
                (2 * length(transcripts))))
            row <- 1
            start_position <- 1
            # make a chain for each transcript in the chromosome
            for (t in transcripts) {
                gtf_t <-
                    gtf_hold[(elementMetadata(gtf_hold)[,"transcript_id"] == t)]
                length_t <- sum(width(ranges(gtf_t)))
                end_position <- start_position + length_t - 1
                first_line <- c(paste0("chain 42 ", chr, " ",
                    seqlengths(seqinf)[which(seqnames(seqinf) == chr)],
                    " ", str, " ", start(ranges(gtf_t))[1], " ",
                    start(ranges(gtf_t))[length(start(ranges(gtf_t)))] +
                    width(ranges(gtf_t))[length(width(ranges(gtf_t)))] - 1,
                    " ", chr, "_", sign, "_", chrom_suffix, " ",
                    length_chr, " ", str, " ", start_position, " ",
                    end_position, " ", t), "", "")
                txpt <- matrix(NA, nrow = length(start(ranges(gtf_t))), ncol= 3)
                txpt[, 1] <- width(ranges(gtf_t))
                # for cases where there is only 1 exon in a transcript
                if (length(width(ranges(gtf_t))) > 1) {
                    txpt[, 2] <-
                        c(start(ranges(gtf_t))[2:length(start(ranges(gtf_t)))] -
                        (width(ranges(gtf_t))[seq(length(width(ranges(gtf_t))) -
                                                                        1)] +
                        start(ranges(gtf_t))[seq(length(width(ranges(gtf_t))) -
                                                                1)]) - 1, "")
                    txpt[, 3] <- 0
                    txpt[nrow(txpt), 3] <- ""
                }
                chrom_chains[row, ] <- first_line
                chrom_chains[(row + 1):(row + nrow(txpt)), ] <- txpt
                chrom_chains[row + nrow(txpt) + 1, ] <- c("", "", "")
                start_position <- end_position
                row <- row + nrow(txpt) + 2
            }
            # append the chains for the chromosome to the chain file
            write.table(chrom_chains,
                file = out_chain_name, append = TRUE, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE, na = ""
            )
            if (verbose == TRUE) {
                print(paste("Chromosome", chr, "strand", str,"complete"))
            }
        }
    }
    # import doesn't work without a truly blank line between chains
    remove_blanks <- readLines(out_chain_name)
    remove_blanks <- gsub("\t\t", "", remove_blanks)
    write(remove_blanks, file = out_chain_name, sep = "")
}
