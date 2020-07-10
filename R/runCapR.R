#' @title runCapR
#'
#' @description Runs CapR
#'
#' @param dir_CapR Directory in which CapR is located. Default current directory
#' @param in_file An .fa file like from ExtractTranscriptomeSequence that is a
#' list of fasta sequences to be folded. Required
#' @param out_file Name of the CapR output file of nucleotide folding
#' probabilities. Default is in_file prefix.out
#' @param max_dist Maximum distance between folded nucleotides in sequences.
#' Recommeded between 50 and 100, with higher values yielding potentially more
#' accurate results but dramatically increasing run time. Default 100.
#'
#' @export

runCapR <- function(dir_CapR = ".",
                    in_file,
                    out_file = NA,
                    max_dist = 100){
    if(is.na(out_file)){
        out_file = paste0(substr(in_file, 0, nchar(in_file)-2), "out")
    }
    system2(paste0(dir_CapR, "/CapR"), args= c(in_file, out_file, max_dist))
}
