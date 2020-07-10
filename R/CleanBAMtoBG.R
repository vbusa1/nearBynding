#' @title CleanBAMtoBG
#'
#' @description Writes a script to convert a BAM file to a clean bedGraph file.
#'
#' @param in_bam Name of sorted BAM file to be converted to a bedGraph file.
#' Required.
#' @param out_bedGraph Name of bedGraph output file, including full directory
#' path. Default in_bam prefix.
#' @param unwanted_chromosomes A vector of unwanted chromosomes that are present
#' in the BAM file.
#'
#' @notes Requires bedtools in the local environment to run.
#'
#' @export

CleanBAMtoBG<-function(in_bam,
                       out_bedGraph = NA,
                       unwanted_chromosomes = NULL){
    if(is.na(out_bedGraph)){
        out_bedGraph <- paste0(substr(in_bam, 0, nchar(in_bam)-3), "bedGraph")
    }
    indexBam(in_bam)
    if(length(unwanted_chromosomes) > 0){
        tmp <- tempfile(tmpdir=".")
        idxstatsBam(in_bam) %>% dplyr::select(1) %>%
            unlist() %>% as.vector() -> hold

        matches <- unique(grep(paste(unwanted_chromosomes, collapse="|"),
                                hold, value=TRUE, invert = TRUE))
        filter <- FilterRules(list(reads = function(x) x$rname %in% matches))
        filterBam(in_bam, tmp, param = ScanBamParam(what = "rname"),
                  filter = filter)
        system2("bedtools", paste0("genomecov -ibam ", tmp, " -bg > ",
                                  out_bedGraph))
        system2("rm", tmp)
        system2("rm", paste0(tmp, ".bai"))
    }else{
        system2("bedtools", paste0("genomecov -ibam ", in_bam, " -bg > ",
                                   out_bedGraph))
    }
}
