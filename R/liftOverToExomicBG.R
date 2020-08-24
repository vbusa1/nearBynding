#' @title liftOverToExomicBG
#'
#' @description Lifts features such as CLIP-seq reads or RNA structure
#' annotations from genome to transcriptome.
#'
#' @param input A single input file name or a vector of input file names
#' in the format of c(forward_reads, reverse_reads) for strand-separated
#' alignments. Files must be BED or bedGraph format. Required
#' @param chain The name of the chain file to be used for liftOver. Format
#' should be like chain files derived from getChainChrSize function.
#' Required
#' @param chrom_size Name of chromosome size file. File must be in two-column
#' format without a header where first column is chromosome name and second
#' column is chromosome length, as from liftOverToExomicBG. Required.
#' @param output_bg The name of the lifted-over output bedGraph file. Required.
#' @param format File type of input file(s). Recommended "BED" or "bedGraph".
#' Default "bedGraph"
#'
#' @return writes lifted-over bedGraph file
#'
#' @examples
#' ## first, get chain file
#' load(system.file("extdata/transcript_list.Rda", package="nearBynding"))
#' gtf<-system.file("extdata/Homo_sapiens.GRCh38.chr4&5.gtf",
#'                 package="nearBynding")
#' GenomeMappingToChainFile(genome_gtf = gtf,
#'                         out_chain_name = "test.chain",
#'                         RNA_fragment = "three_prime_utr",
#'                         transcript_list = transcript_list,
#'                         alignment = "hg38")
#' ## and chain file chromosome sizes
#' getChainChrSize(chain = "test.chain",
#'                out_chr = "chr4and5_3UTR.size")
#'
#' ## get bedGraph file
#' chr4and5_sorted.bedGraph<-system.file("extdata/chr4and5_sorted.bedGraph",
#'                                      package="nearBynding")
#'
#' liftOverToExomicBG(input = chr4and5_sorted.bedGraph,
#'                   chain = "test.chain",
#'                   chrom_size = "chr4and5_3UTR.size",
#'                   output_bg = "chr4and5_liftOver.bedGraph")
#'
#' @importFrom magrittr '%>%'
#' @importFrom S4Vectors elementNROWS
#' @importFrom plyranges bind_ranges
#' @importFrom GenomicRanges countOverlaps
#' @importFrom utils read.delim
#' @importFrom rtracklayer import
#'
#'
#' @export

liftOverToExomicBG <- function(input,
                                chain,
                                chrom_size,
                                output_bg,
                                format = "bedGraph") {
    # allow flexible naming systems for forward and reverse files
    if (length(input) == 2) {
        f_bg <- input[1]
        r_bg <- input[2]
        plus <- import(f_bg, format = format)
        minus <- import(r_bg, format = format)
    } else if (length(input) == 1) {
        f_bg <- input
        plus <- import(f_bg, format = format)
        minus <- plus
    } else {
        print("input must be a single file or a vector of file names in the
            format of c(forward_reads, reverse_reads) of type BED or
            bedGraph.")
    }
    BiocGenerics::strand(plus) <- "+"
    BiocGenerics::strand(minus) <- "-"
    chain_object <- import(chain, exclude = "junk")
    # Perform liftOver for + and - strands given an exome chain object;
    # reduce to only relevant intervals and combine into one bedGraph
    liftOver_plus <- liftOver(plus, chain_object)
    # remove any blank intervals (shouldn't be any > 1)
    liftOver_plus <- liftOver_plus[(elementNROWS(range(liftOver_plus))==1L)] %>%
        unlist()
    liftOver_minus <- liftOver(minus, chain_object)
    liftOver_minus<-liftOver_minus[(elementNROWS(range(liftOver_minus))==1L)]%>%
        unlist()
    liftOver_plus <- liftOver_plus[grepl("plus",
                                    as.character(liftOver_plus@seqnames))]
    liftOver_minus <- liftOver_minus[grepl("minus",
                                    as.character(liftOver_minus@seqnames))]
    liftOver <- bind_ranges(liftOver_minus, liftOver_plus)
    # must re-introduce chromosome lengths into seqinfo
    seq_info <- read.delim(chrom_size, header = FALSE)
    colnames(seq_info) <- c("chr", "length")
    liftOver@seqinfo@seqlengths <-
        seq_info[which(seq_info$chr %in% liftOver@seqinfo@seqnames), 2] %>%
        as.character() %>%
        as.integer()
    # remove duplicates
    liftOver <- unique(liftOver)
    liftOver <- liftOver[countOverlaps(liftOver, liftOver) <= 1L]
    # output exome bg file
    export(liftOver, con = output_bg, format = "bedGraph")
}
