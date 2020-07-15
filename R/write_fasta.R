#' @title write_fasta
#'
#' @description Writes a FASTA file from a vector of sequences
#'
#' @param sequences A vector of sequences
#' @param names A vector of names corresponding to the sequences
#' @param file.out Name of output FASTA file
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
