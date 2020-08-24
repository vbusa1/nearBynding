#' @title getChainChrSize
#'
#' @description Output a table of mapped chromosome names and lengths from a
#' chain file.
#'
#' @param chain The name of the chain file for which chromosome sizes should be
#' determined and output; a string. Required.
#' @param out_chr Name of the chromosome names and lengths table file; a string.
#' Required.
#'
#' @return writes a two-column tab-delineated file without a header containing
#' chromosome names and lengths for a given chain file
#'
#' @examples
#' ## first, make the chain file
#' load(system.file("extdata/transcript_list.Rda", package="nearBynding"))
#' gtf<-system.file("extdata/Homo_sapiens.GRCh38.chr4&5.gtf",
#'                 package="nearBynding")
#' GenomeMappingToChainFile(genome_gtf = gtf,
#'                         out_chain_name = "test.chain",
#'                         RNA_fragment = "three_prime_utr",
#'                         transcript_list = transcript_list,
#'                         alignment = "hg38")
#'
#' getChainChrSize(chain = "test.chain",
#'                out_chr = "chr4and5_3UTR.size")
#'
#' @importFrom magrittr '%>%'
#' @importFrom utils write.table
#'
#' @export

getChainChrSize <- function(chain,
                            out_chr) {
    chain_object <- import(chain, exclude = "junk")
    seqnames_liftedOver <- c()
    for (chr in names(chain_object)) {
        seqnames_liftedOver <- append(
            seqnames_liftedOver,
            unique(space(chain_object@listData[[chr]])))
    }
    lines <- c()
    for (chr in seqnames_liftedOver) {
        con <- file(chain, "r")
        while (TRUE) {
            line <- readLines(con, 1)
            if (length(line) == 0) {
                break
            } else if (grepl(chr, line)) {
                lines <- c(lines, line)
                break
            }
        }
        close(con)
    }
    seq_info <- do.call(rbind.data.frame, lapply(lines,
        function(x) {
            strsplit(x, " ") %>%
                unlist()
        }))[, 8:9] %>%
        unique()
    write.table(seq_info, out_chr,
        sep = "\t",
        col.names = FALSE, row.names = FALSE, quote = FALSE
    )
}
