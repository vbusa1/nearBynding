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
    file.copy(system.file("extdata/chr4and5_3UTR_stem_liftOver~chr4and5_liftOver.dist",
                          package="nearBynding"), dir)
    file.copy(system.file("extdata/chr4and5_3UTR_internal_liftOver~chr4and5_liftOver.dist",
                          package="nearBynding"), dir)
    file.copy(system.file("extdata/chr4and5_3UTR_hairpin_liftOver~chr4and5_liftOver.dist",
                          package="nearBynding"), dir)
    file.copy(system.file("extdata/chr4and5_3UTR_multibranch_liftOver~chr4and5_liftOver.dist",
                          package="nearBynding"), dir)
    file.copy(system.file("extdata/chr4and5_3UTR_exterior_liftOver~chr4and5_liftOver.dist",
                          package="nearBynding"), dir)
    file.copy(system.file("extdata/chr4and5_3UTR_bulge_liftOver~chr4and5_liftOver.dist",
                          package="nearBynding"), dir)
    file.copy(system.file("extdata/chr4and5_3UTR.out",
                          package="nearBynding"), dir)
}
