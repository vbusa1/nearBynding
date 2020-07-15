#' @title assessGrouping
#'
#' @description Assess grouping of samples assigned to the same category
#' relative to random.
#'
#' @param distances Data frame with at least three columns where the first three
#' columns are sample 1 name, sample 2 name, and the distance between them.
#' @param annotations Data frame with at least two columns where the first two
#' columns are sample name and the category of the sample for grouping. Sample
#' names must match sample 1 and sample 2 names in distances data frame.
#' @param measurement The measurement for comparison between cases and controls
#' and statistical analysis ("mean", "max", or "min). Default "mean"
#' @param output What information will be returned, either a list of test and
#' control measurement distances ("measurements"), the p-value of the
#' Kolmogorov-Smirnov test comparing test and control distributions
#' ("KS.pvalue"), or a ggplot object plotting the test and control distributions
#' ("plot"). Default "KS.pvalue"
#' @param ctrl_iterations The number of iterations to test for the control
#' distribution. Default 10000.
#'
#' @export

assessGrouping<-function(distances,
                         annotations,
                         measurement = "mean",
                         output = "KS.pvalue",
                         ctrl_iterations = 10000){
    if(measurement != "mean" & measurement != "max" & measurement != "min"){
        stop("measurement argument must equal 'mean', 'max', or 'min'
             for statistical analysis.")
    }
    if(output != "measurements" & output != "plot" & output != "KS.pvalue"){
        stop("output argument must equal 'measurements', 'plot, or 'KS.pvalue")
    }
    if(ctrl_iterations < 20){
        stop("ctrl_iterations argument is too small; consider increasing
             iteration count")
    }
    if(ncol(distances) < 3){
        stop("The first three columns of the distances data frame should be
             sample 1 name, sample 2 name, and the distance between them")
    }
    if(ncol(annotations) < 2){
        stop("The first two columns of the annotations data frame should be
             sample name and the category of the sample for grouping")
    }

    colnames(distances)[c(1,2,3)]<-c("Var1", "Var2", "dist")
    colnames(annotations)[c(1,2)]<-c("sample", "category")

    # filter for distances with proper annotations
    test_distances<-dplyr::filter(distances, Var1 %in% annotations$sample &
                               Var2 %in% annotations$sample)
    test_distances<-test_distances[!duplicated(t(apply(test_distances[,c(1,2)],
                                        1, sort))),] # remove duplicate pairs
    test_distances<-dplyr::filter(test_distances, Var1 != Var2)
    test_categories<-dplyr::filter(annotations, sample %in% distances$Var1 |
                                sample %in% distances$Var2)
    if(nrow(test_distances) < 2 | nrow(test_categories) < 2){
        stop("Sample names in the distance data frame do not match sample names
             in annotations data frame.")
    }

    # get categories for which there is at least two samples
    unique_categories <- test_categories %>%
        dplyr::group_by(category) %>%
        dplyr::summarise(no_rows = length(category), .groups = 'drop') %>%
        dplyr::filter(no_rows > 1)

    test_output<-lapply(unique_categories$category, function(i){
        subset_samples<-dplyr::filter(test_categories, category == i)$sample %>%
            as.character()
        subset_distances <- dplyr::filter(test_distances, Var1 %in%
                                subset_samples & Var2 %in% subset_samples)$dist

        if(measurement == "mean"){
            return(mean(subset_distances))
        } else if(measurement == "max"){
            return(max(subset_distances))
        } else{
            return(min(subset_distances))
        }
    }) %>% unlist()

    ctrl_output<- lapply(seq_len(ctrl_iterations), function(i){
        subset_distances<- sample(test_distances$dist ,
                                  sample(unique_categories$no_rows, 1))

        if(measurement == "mean"){
            return(mean(subset_distances))
        } else if(measurement == "max"){
            return(max(subset_distances))
        } else{
            return(min(subset_distances))
        }
    }) %>% unlist()

    if(output == "measurements"){
        return(list(test_output = test_output,
                    ctrl_output = ctrl_output))
    } else if(output == "KS.pvalue"){
        return(suppressWarnings(ks.test(test_output, ctrl_output,
                       alternative = "greater")$p.value))
    } else if(output == "plot"){
        group <- c(rep("ctrl", length(ctrl_output)),
                   rep("test", length(test_output)))
        dat_plot <- data.frame(KSD = c(ctrl_output, test_output), group = group)
        plot<-ggplot(dat_plot,
                     aes(x = KSD, group = group, color = group)) +
            stat_ecdf(size=1, na.rm = TRUE) +
            theme_classic() +
            theme(legend.position ="top") +
            ylab("cummulative distribution") +
            xlab("distance") +
            theme(legend.position = NULL)
        return(plot)
    }
}
