#' @title processCapRout
#'
#' @description Creates context-separated bedGraph files of CapR output for
#' genome and transcriptome alignments.
#'
#' @param CapR_outfile Name of CapR output file. Required
#' @param output_prefix Prefix to be appended to all output files. Required
#' @param chrom_size Name of chromosome size file in two-column format without
#' a header where first column is chromosome name and second column is
#' chromosome length, as from liftOverToExomicBG Required
#' @param genome_gtf The name of the GTF/GFF file that contains all exome
#' annotations. Required
#' @param RNA_fragment RNA component of interest. Options depend on the gtf
#' file but often include "gene", "transcript", "exon", "CDS", "five_prime_utr",
#' and/or "three_prime_utr". Default "exon" for the whole exome.
#' @param chain The name of the chain file to be used. Format should be like
#' chain files derived from GRangesMappingToChainFile. Required
#'
#' @export

processCapRout<-function(CapR_outfile,
                       output_prefix,
                       chrom_size,
                       genome_gtf,
                       RNA_fragment,
                       chain){
    folded_transcripts_appended <- readLines(CapR_outfile)

    transcript_names <- folded_transcripts_appended[which(
        substr(folded_transcripts_appended, 1, 1) == ">")]
    transcript_names <- lapply(transcript_names, function(x){
        substr(x, 2, nchar(x))}) %>% unlist()

    bulge<-folded_transcripts_appended[which(
        substr(folded_transcripts_appended, 1, 5) == "Bulge")]
    exterior<-folded_transcripts_appended[which(
        substr(folded_transcripts_appended, 1, 8) == "Exterior")]
    hairpin<-folded_transcripts_appended[which(
        substr(folded_transcripts_appended, 1, 7) == "Hairpin")]
    internal<-folded_transcripts_appended[which(
        substr(folded_transcripts_appended, 1, 8) == "Internal")]
    multibranch<-folded_transcripts_appended[which(
        substr(folded_transcripts_appended, 1, 11) == "Multibranch")]
    stem<-folded_transcripts_appended[which(
        substr(folded_transcripts_appended, 1, 4) == "Stem")]

    cut_transcripts<-function(x){
        transcript<-strsplit(x, ' ') %>% unlist()
        transcript<-transcript[2:length(transcript)]
        return(transcript)
    }

    bulge_cut<-lapply(bulge, cut_transcripts) %>% unlist()
    exterior_cut<-lapply(exterior, cut_transcripts) %>% unlist()
    hairpin_cut<-lapply(hairpin, cut_transcripts) %>% unlist()
    internal_cut<-lapply(internal, cut_transcripts) %>% unlist()
    multibranch_cut<-lapply(multibranch, cut_transcripts) %>% unlist()
    stem_cut<-lapply(stem, cut_transcripts) %>% unlist()

    # double check that transcriptome length total (chr size file)
    # equals number of folded transcript positions
    chr_size<-read.delim(chrom_size, header=F)

    if(sum(chr_size$V2)-length(bulge_cut) != 0){
        stop("The length of the CapR output does not equal the length of
             chromosome size. Make sure that the input chromosome size
             corresponds to the folded CapR transcripts.")
    }

    gtf <- import(genome_gtf)
    gtf <- gtf %>% filter(type == RNA_fragment)

    # write bedgraph file for each transcript
    # create bedGraph for each structure separately
    bg_shell <- as.data.frame(matrix(NA, ncol = 4, nrow = sum(chr_size$V2)))
    colnames(bg_shell) <- c("chr", "start", "end", "score")
    n <- 1
    for(t in transcript_names){
        gtf_t<-gtf[(elementMetadata(gtf)[,"transcript_id"] == t)]
        bg_shell[n:(n + sum(gtf_t@ranges@width)-1),"chr"]<-
            gtf_t@seqnames@values %>% as.character()
        if(gtf_t@strand@values[1] == "+"){
            for(unit in 1:length(gtf_t)){
                gtf_t_unit<-gtf_t[unit]
                bg_shell[n:(n+gtf_t_unit@ranges@width-1),"start"]<-
                    gtf_t_unit@ranges@start:(gtf_t_unit@ranges@start +
                                                 gtf_t_unit@ranges@width-1)
                n<-n+gtf_t_unit@ranges@width
            }
        }else{
            for(unit in length(gtf_t):1){
                gtf_t_unit<-gtf_t[unit]
                bg_shell[n:(n+gtf_t_unit@ranges@width-1),"start"]<-
                    rev(gtf_t_unit@ranges@start:(gtf_t_unit@ranges@start +
                                                     gtf_t_unit@ranges@width-1))
                n<-n+gtf_t_unit@ranges@width
            }
        }
    }
    bg_shell$end <- bg_shell$start
    if(nchar(bg_shell[1,"chr"]) < 3){
        bg_shell$chr <- paste0("chr", bg_shell$chr)
    }

    # bulge
    df_bulge<-bg_shell
    df_bulge$score<-bulge_cut %>% as.numeric()
    GRanges_bulge<-makeGRangesFromDataFrame(df_bulge,
                                            keep.extra.columns = T)
    export.bedGraph(GRanges_bulge, paste0(output_prefix,
                                          "_bulge.bedGraph"))
    liftOverToExomicBG(input = paste0(output_prefix,
                                         "_bulge.bedGraph"),
                       chain = chain,
                       output_bg = paste0(output_prefix,
                                          "_bulge_liftOver.bedGraph"),
                       write_chr = F)
    # exterior
    df_exterior<-bg_shell
    df_exterior$score<-exterior_cut %>% as.numeric()
    GRanges_exterior<-makeGRangesFromDataFrame(df_exterior,
                                               keep.extra.columns = T)
    export.bedGraph(GRanges_exterior, paste0(output_prefix,
                                             "_exterior.bedGraph"))
    liftOverToExomicBG(input = paste0(output_prefix,
                                         "_exterior.bedGraph"),
                       chain = chain,
                       output_bg = paste0(output_prefix,
                                          "_exterior_liftOver.bedGraph"),
                       write_chr = F)
    # hairpin
    df_hairpin<-bg_shell
    df_hairpin$score<-hairpin_cut %>% as.numeric()
    GRanges_hairpin<-makeGRangesFromDataFrame(df_hairpin,
                                              keep.extra.columns = T)
    export.bedGraph(GRanges_hairpin, paste0(output_prefix,
                                            "_hairpin.bedGraph"))
    liftOverToExomicBG(input = paste0(output_prefix,
                                         "_hairpin.bedGraph"),
                       chain = chain,
                       output_bg = paste0(output_prefix,
                                          "_hairpin_liftOver.bedGraph"),
                       write_chr = F)
    # internal
    df_internal<-bg_shell
    df_internal$score<-internal_cut %>% as.numeric()
    GRanges_internal<-makeGRangesFromDataFrame(df_internal,
                                               keep.extra.columns = T)
    export.bedGraph(GRanges_internal, paste0(output_prefix,
                                             "_internal.bedGraph"))
    liftOverToExomicBG(input = paste0(output_prefix,
                                         "_internal.bedGraph"),
                       chain = chain,
                       output_bg = paste0(output_prefix,
                                          "_internal_liftOver.bedGraph"),
                       write_chr = F)
    # multibranch
    df_multibranch<-bg_shell
    df_multibranch$score<-multibranch_cut %>% as.numeric()
    GRanges_multibranch<-makeGRangesFromDataFrame(df_multibranch,
                                                  keep.extra.columns = T)
    export.bedGraph(GRanges_multibranch, paste0(output_prefix,
                                                "_multibranch.bedGraph"))
    liftOverToExomicBG(input = paste0(output_prefix,
                                         "_multibranch.bedGraph"),
                       chain = chain,
                       output_bg = paste0(output_prefix,
                                          "_multibranch_liftOver.bedGraph"),
                       write_chr = F)
    # stem
    df_stem<-bg_shell
    df_stem$score<-stem_cut %>% as.numeric()
    GRanges_stem<-makeGRangesFromDataFrame(df_stem,
                                           keep.extra.columns = T)
    export.bedGraph(GRanges_stem, paste0(output_prefix,
                                         "_stem.bedGraph"))
    liftOverToExomicBG(input = paste0(output_prefix, "_stem.bedGraph"),
                       chain = chain,
                       output_bg = paste0(output_prefix,
                                          "_stem_liftOver.bedGraph"),
                       write_chr = F)
}
