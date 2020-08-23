#' @title CleanBEDtoBG
#'
#' @description Writes a script to convert a BED file to a clean bedGraph file.
#'
#' @param in_bed Name of sorted BAM file to be converted to a bedGraph file.
#' Required.
#' @param out_bedGraph Name of bedGraph output file, including full directory
#' path. Default in_bam prefix.
#' @param unwanted_chromosomes A vector of unwanted chromosomes that are present
#' in the BAM file.
#' @param alignment The human genome alignment used, either "hg19" or "hg38".
#' Default "hg19"
#'
#' @return deposits bedGraph from BED in same directory
#'
#' @examples
#' bam <- system.file("extdata/chr4and5.bam", package="nearBynding")
#' out_bed <- "bamto.bed"
#' \dontrun{
#'     system2("bedtools", paste0("bamtobed -i ", bam, " > ", out_bed))
#'
#'     CleanBEDtoBG(in_bed = out_bed,
#'                  alignment = "hg38")
#' }
#'
#' @export

CleanBEDtoBG <- function(in_bed, out_bedGraph = NA,
                            unwanted_chromosomes = NULL,
                            alignment = "hg19") {
    if (alignment == "hg19") {
        genome <- system.file("extdata/hg19.genome", package = "nearBynding")
    }
    if (alignment == "hg38") {
        genome <- system.file("extdata/hg38.genome", package = "nearBynding")
    }
    if (!(alignment %in% c("hg19", "hg38"))) {
        stop("Acceptable alignments are 'hg19' and 'hg38'")
    }
    # sort bed
    out_bed_sorted <- paste0(substr(in_bed, 0, nchar(in_bed) - 4),"_sorted.bed")
    system2("sort", paste0("-k 1,1 -k2,2n ", in_bed, " > ", out_bed_sorted))
    if (is.na(out_bedGraph)) {
        out_bedGraph <- paste0(substr(in_bed, 0, nchar(in_bed) - 3), "bedGraph")
    }
    if (length(unwanted_chromosomes) > 0) {
        tmp <- tempfile(tmpdir = ".")
        system2("sort", paste0("'/", paste(unwanted_chromosomes,
                                collapse = "/d;/"),
                                "/d' ", out_bed_sorted, " > ", tmp))
        system2("bedtools", paste0("genomecov -i ",
                                tmp, " -g ", genome, " -bg > ",out_bedGraph))
        system2("rm", tmp)
    } else {
        system2("bedtools",paste0("genomecov -i ",out_bed_sorted, " -g ",genome,
                                " -bg > ", out_bedGraph))
    }
    system2("rm", out_bed_sorted)
}
