#' @title runStereogene
#'
#' @description Writes a StereoGene script in the working directory
#'
#' @param track_files Vector of at least two track or interval file names to be
#' pairwise-correlated by StereoGene. Required.
#' @param name_config Name of corresponding configuration file; a string.
#' Required
#' @param pcorProfile Name of track file name for partial correlation; a string.
#' More information for this can be found in the StereoGene README. Optional
#' @param confounder Confounder file name; a string. More information for this
#' can be found in the StereoGene README. Optional
#' @param nShuffle Permutations used to estimate error. Default 1000.
#' @param get_error Whether to calculate the standard error of background
#' permutations from nShuffle. FALSE will save calculation time. Default TRUE
#'
#' @return generates StereoGene output files in directory
#'
#' @examples
#' runStereogene(track_files = c("chr4and5_3UTR_stem_liftOver.bedGraph",
#'                              "chr4and5_liftOver.bedGraph"),
#'              name_config = "chr4and5_3UTR.cfg")
#'
#' @export

runStereogene <- function(track_files,
                            name_config,
                            pcorProfile = NULL,
                            confounder = NULL,
                            nShuffle = 1000,
                            get_error = TRUE) {
    if (length(track_files) != 2) {
        stop("Must have two track or interval files for correlation.")
    }
    if (length(name_config) < 1) {
        stop("Must provide a name for the configuration file.")
    }
    tracks <- paste(track_files, collapse = " ")
    config <- paste0("-cfg ", name_config, " ")
    partial_corr <- ""
    if (length(pcorProfile) == 1) {
        partial_corr <- paste0("-pcorProfile ", pcorProfile, " ")
    }
    confound <- ""
    if (length(confounder) == 1) {
        confound <- paste0("-confounder ", confounder, " ")
    }
    if(.is_StereoGene_installed()){
        window<- as.numeric(substring(grep("wSize=",
                                           readLines(name_config),
                                           value = TRUE), 7))
        if(get_error == TRUE){
        Bkg<-matrix(nrow = window, ncol = nShuffle)
        for(n in 1:nShuffle){
            .StereoGene_run(paste0("-nShuffle 1 ", config, partial_corr,
                                   confound, tracks))
            Bkg[,n]<-read.table(paste0(gsub("\\..*", "", track_files[1]), "~",
                                       gsub("\\..*", "", track_files[2]),
                                       ".dist"), header = TRUE)$Bkg
            if(n > 1){
                if(Bkg[1,n] == Bkg[1,(n-1)]){ # sometimes get repeats
                .StereoGene_run(paste0("-nShuffle 1 ", config, partial_corr,
                                       confound, tracks))
                Bkg[,n]<-read.table(paste0(gsub("\\..*", "",
                                                track_files[1]), "~",
                                           gsub("\\..*", "",
                                                track_files[2]),
                                           ".dist"), header = TRUE)$Bkg
                }
                }
        }
        dist<-data.frame(x = read.table(paste0(gsub("\\..*", "",
                                                    track_files[1]), "~",
                                               gsub("\\..*", "",
                                                    track_files[2]),
                                               ".dist"), header = TRUE)$x,
                         Fg = read.table(paste0(gsub("\\..*", "",
                                                     track_files[1]), "~",
                                                gsub("\\..*", "",
                                                     track_files[2]),
                                                ".dist"), header = TRUE)$Fg,
                         Bkg = apply(Bkg, 1, mean),
                         Bkg_se = apply(Bkg, 1, sd)/nShuffle)
        }else{
        .StereoGene_run(paste0("-nShuffle ", nShuffle, " ", config,
                               partial_corr, confound, tracks))
        dist<-data.frame(x = read.table(paste0(gsub("\\..*", "",
                                                    track_files[1]), "~",
                                               gsub("\\..*", "",
                                                    track_files[2]),
                                               ".dist"), header = TRUE)$x,
                         Fg = read.table(paste0(gsub("\\..*", "",
                                                     track_files[1]), "~",
                                                gsub("\\..*", "",
                                                     track_files[2]),
                                                ".dist"), header = TRUE)$Fg,
                         Bkg = read.table(paste0(gsub("\\..*", "",
                                                     track_files[1]), "~",
                                                gsub("\\..*", "",
                                                     track_files[2]),
                                                ".dist"), header = TRUE)$Bkg,
                         Bkg_se = 0)
        }
        write.table(dist, file = paste0(gsub("\\..*", "", track_files[1]), "~",
                                  gsub("\\..*", "", track_files[2]),
                                  ".dist"), quote = F, row.names = FALSE)
    }else{return("Please install Stereogene and place in working PATH")}
}
