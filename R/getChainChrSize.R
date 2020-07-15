#' @title getChainChrSize
#'
#' @description Output a table of mapped chromosome names and lengths from a
#' chain file.
#'
#' @param chain The name of the chain file for which chromosome sizes should be
#' determined and output. Required.
#' @param out_chr Name of the chromosome names and lengths table file. Required
#'
#' @export

getChainChrSize <- function(chain,
                            out_chr) {
  chain_object <- import(chain, exclude = "junk")

  seqnames_liftedOver <- c()
  for (chr in names(chain_object)) {
    seqnames_liftedOver <- append(
      seqnames_liftedOver,
      unique(space(chain_object@listData[[chr]]))
    )
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
  seq_info <- do.call(rbind.data.frame, lapply(
    lines,
    function(x) {
      strsplit(x, " ") %>%
        unlist()
    }
  ))[, 8:9] %>%
    unique()
  write.table(seq_info, out_chr,
    sep = "\t",
    col.names = FALSE, row.names = FALSE, quote = FALSE
  )
}
