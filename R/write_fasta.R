#' @title write_fasta
#'
#' @description Writes a FASTA file from a vector of sequences
#'
#' @param sequences A vector of sequences
#' @param names A vector of names corresponding to the sequences
#' @param file.out Name of output FASTA file; a string
#'
#' @return writes FASTA file into directory
#'
#' @examples
#'sequences<-c(paste0(sample(c("A", "T", "G", "C"), 20, replace = TRUE),
#'                    collapse = ""),
#'             paste0(sample(c("A", "T", "G", "C"), 20, replace = TRUE),
#'                    collapse = ""),
#'             paste0(sample(c("A", "T", "G", "C"), 20, replace = TRUE),
#'                    collapse = ""))
#'write_fasta(sequences,
#'            c("one", "two", "three"),
#'            "test.fa")
#'
#' @export
write_fasta <- function(sequences,
                        names,
                        file.out) {
    # create output file
    outfile <- file(description = file.out, open = "w")
    invisible(lapply(
        seq(length(sequences)),
        function(x) {
            # write one sequence into output file
            sequence <- as.character(sequences[[x]])
            writeLines(paste(">", names[x], sep = ""), outfile)
            writeLines(sequence, outfile)
        }
    ))
    close(outfile)
}
