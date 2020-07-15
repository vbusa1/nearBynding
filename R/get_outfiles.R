#' @title get_outfiles
#'
#' @description Copy files necessary to complete the vignette onto the local
#' machine in cases where Stereogene, CapR, or bedtools are not available.
#'
#' @param dir Directory into which files ought to be stored. Default current
#' work directory.
#'
#' @export

get_outfiles<-function(dir = "."){
    file.copy(system.file("extdata/stem_liftOver_chr4and5.dist",
                          package="nearBynding"), dir)
    file.rename(paste0(dir, "/stem_liftOver_chr4and5.dist"),
            paste0(dir, "/chr4and5_3UTR_stem_liftOver~chr4and5_liftOver.dist"))
    file.copy(system.file("extdata/bulge_liftOver_chr4and5.dist",
                          package="nearBynding"), dir)
    file.rename(paste0(dir, "/bulge_liftOver_chr4and5.dist"),
                paste0(dir,
                       "/chr4and5_3UTR_bulge_liftOver~chr4and5_liftOver.dist"))
    file.copy(system.file("extdata/exterior_liftOver_chr4and5.dist",
                          package="nearBynding"), dir)
    file.rename(paste0(dir, "/exterior_liftOver_chr4and5.dist"),
                paste0(dir,
                     "/chr4and5_3UTR_exterior_liftOver~chr4and5_liftOver.dist"))
    file.copy(system.file("extdata/internal_liftOver_chr4and5.dist",
                          package="nearBynding"), dir)
    file.rename(paste0(dir, "/internal_liftOver_chr4and5.dist"),
                paste0(dir,
                     "/chr4and5_3UTR_internal_liftOver~chr4and5_liftOver.dist"))
    file.copy(system.file("extdata/hairpin_liftOver_chr4and5.dist",
                          package="nearBynding"), dir)
    file.rename(paste0(dir, "/hairpin_liftOver_chr4and5.dist"),
                paste0(dir,
                      "/chr4and5_3UTR_hairpin_liftOver~chr4and5_liftOver.dist"))
    file.copy(system.file("extdata/multibranch_liftOver_chr4and5.dist",
                          package="nearBynding"), dir)
    file.rename(paste0(dir, "/multibranch_liftOver_chr4and5.dist"),
                paste0(dir,
                "/chr4and5_3UTR_multibranch_liftOver~chr4and5_liftOver.dist"))
}
