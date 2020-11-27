#==============================================================================
#
#                   Switch between nondetect representations
#
#==============================================================================


#' Switch between nondetect representations
#' @param x         matrix
#' @param verbose   logical(1)
#' @return Updated SummarizedExperiment
#' @examples
#' require(magrittr)
#'
#' # 0 -> NA (proteingroups LFQ intensities)
#'
#' # NaN -> NA (proteingroups ratios)
#'     file <- load_data('billing16.proteingroups.txt')
#'     object <- read_proteingroups(file)
#'     exprs(object) %<>% nan_to_na(verbose=TRUE)
#'
#' # -Inf -> NA (log2 transformed proteingroups LFQ intensity)
#'
#' # NA -> 0
#'     file <- load_data('halama18.metabolon.xlsx')
#'     object <- exprs(read_metabolon(file)
#'     exprs(object) %<>% na_to_zero(verbose=TRUE)
#' @noRd
zero_to_na <- function(x, verbose = FALSE){
    selector <- x == 0
    if (any(c(selector), na.rm = TRUE)){
        if (verbose) cmessage(
                    paste0('\t\tReplace 0 -> NA for %d/%d values ',
                            '(in %d/%d features and %d/%d samples)'),
                    sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                    sum(rowAnys(selector)), nrow(x),
                    sum(colAnys(selector)), ncol(x))
        x[selector] <- NA_real_
    }
    x
}


#' @noRd
nan_to_na <- function(x, verbose = FALSE){
    selector <- is.nan(x)
    if (any(c(selector), na.rm = TRUE)){
        if (verbose) cmessage(
                    paste0( '\t\tReplace NaN -> NA for %d/%d values ',
                            '(in %d/%d features and %d/%d samples)'),
                    sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                    sum(rowAnys(selector)), nrow(x),
                    sum(colAnys(selector)), ncol(x))
        x[selector] <- NA_real_
    }
    x
}


na_to_zero <- function(x, verbose = FALSE){
    selector <- is.na(x)
    if (any(selector)){
        if (verbose) cmessage(
                        paste0( '\t\tReplace NA -> 0 for %d/%d values ',
                                '(in %d/%d features and %d/%d samples)'),
                        sum(selector), nrow(selector)*ncol(selector),
                        sum(rowAnys(selector)), nrow(x),
                        sum(colAnys(selector)), ncol(x))
        x[selector] <- 0
    }
    x
}


inf_to_na <- function(x, verbose){
    selector <- is.infinite(x)
    if (any(c(selector), na.rm = TRUE)){
        if (verbose) cmessage(
                    paste0(
                        '\t\tReplace -Inf -> NA for %d/%d values ',
                        '(in %d/%d features and %d/%d samples)'),
                    sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                    sum(rowAnys(selector)), nrow(x),
                    sum(colAnys(selector)), ncol(x))
        x[selector] <- NA_real_
    }
    x
}


minusinf_to_na <- function(x, verbose = FALSE){
    selector <- x==-Inf
    if (any(c(selector), na.rm = TRUE)){
        if (verbose) cmessage(
                    paste0( '\t\tReplace -Inf -> NA for %d/%d values ',
                            '(in %d/%d features and %d/%d samples)'),
                    sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                    sum(rowAnys(selector)), nrow(x),
                    sum(colAnys(selector)), ncol(x))
        x[selector] <- NA_real_
    }
    x
}


na_to_string <- function(x){
    x[is.na(x)] <- ''
    x
}


#=============================================================================
#
#                   normimpute
#                   halfnormimpute
#                   zeroimpute
#
#=============================================================================

#' @rdname halfnormimpute
#' @export
normimpute <- function(x, selector = is.na(x)){
    x[selector] <- rnorm(length(x[selector]), sd = sd(x[!is.na(x)]))
    x
}


#' Impute from half-normal distribution around 0
#' @param x          NA-containing numeric vector
#' @param selector   which values to impute
#' @return numeric vector of same length
#' @examples
#' require(data.table)
#' x <- rnorm(1e5)
#' idx <- runif(length(x))>0.9
#' x[idx] <- NA
#' dt1 <- data.table(value = normimpute(x), distr = 'norm')
#'
#' x <- abs(rnorm(1e5)); x[idx] <- NA
#' dt2 <- data.table(value = halfnormimpute(x), distr = 'halfnorm')
#'
#' x <- abs(rnorm(1e5)); x[idx] <- NA
#' dt3 <- data.table(value = zeroimpute(x), distr = 'zero')
#'
#' require(ggplot2)
#' ggplot(rbind(dt1,dt2,dt3), aes(x=value, fill=distr)) +
#' geom_density(alpha=0.5)
#' @export
halfnormimpute <- function(x, selector = is.na(x)){
    x[selector] <- abs(rnorm(length(x[selector]), sd = 2*sd(x[!is.na(x)])))
    x
}


#' @rdname halfnormimpute
#' @export
zeroimpute <- function(x, selector = is.na(x)){
    x[selector] <- 0
    x
}

#=============================================================================
#
#                        split_by_svar
#
#==============================================================================

#' Split by svar
#' @param object SummarizedExperiment
#' @param svar   svar to split on
#' @return list of SummarizedExperiment
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute=FALSE, plot = FALSE)
#' split_by_svar(object)
#' @export
split_by_svar <- function(object, svar = subgroup){
    svar <- enquo(svar)
    svarstr <- as_name(svar)

    if (is.null(svar)) return(list(object))
    extract_samples  <- function(sg){
                            idx <- sdata(object)[[svarstr]] == sg
                            object[, idx]
                        }
    Map(extract_samples, slevels(object, svarstr))
}



#=============================================================================
#
#                     is_systematic_detect
#                     is_random_detect
#                     is_full_detect
#
#=============================================================================


is_systematic_detect <- function(object, group = subgroup){
    group <- enquo(group)
    split_by_svar(object, !!group) %>%
    lapply(function(x) rowAlls(is.na(exprs(x)))) %>%
    Reduce("|", .)
}

is_random_detect <- function(object, group = subgroup){
    group <- enquo(group)
    rowAnys(is.na(exprs(object))) & !is_systematic_detect(object, !!group)
}


is_full_detect <- function(object){
    rowAlls(!is.na(exprs(object)))
}


#=============================================================================
#
#                     venn_detects
#
#=============================================================================


#' Venn detects
#'
#' Venn diagram full/systematic/random detects
#'
#' @param object SummarizedExperiment
#' @return NULL
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute=FALSE, plot = FALSE)
#' venn_detects(object)
#' @export
venn_detects <- function(object){
    limma::vennDiagram(data.matrix(cbind(
        systematic = is_systematic_detect(object),
        random     = is_random_detect(object),
        full       = is_full_detect(object))))
}


#=============================================================================
#
#                     impute_systematic_nondetects
#
#==============================================================================


#' Impute systematic nondetects
#' @param object SummarizedExperiment
#' @param group  group svar
#' @param fun    imputation function
#' @param plot   TRUE or FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute = FALSE, plot = FALSE)
#' impute_systematic_nondetects(object)
#' @export
impute_systematic_nondetects <- function(
    object, group = subgroup, fun = halfnormimpute, plot = TRUE
){
# Process
    absent <- replicated <- systematic <- NULL
    group <- enquo(group)
    groupstr <- as_name(group)
# Filter
    object %<>% filter_exprs_replicated_in_some_subgroup()
# Impute
    dt  <-  sumexp_to_long_dt(object)[
                , absent     := all(is.na(value)),    c('feature_id', groupstr)
            ][  , replicated := sum(!is.na(value))>1, c('feature_id', groupstr)
            ][  , systematic := any(absent) & any(replicated), 'feature_id']
    set.seed(39)
    dt[, is_imputed := systematic & absent]
    dt[, value := fun(value, is_imputed), by='feature_id']
# Update object
    ff <- fnames(object)
    ss <- snames(object)
    assays(object)$exprs      <- dt2exprs(dt)[ff, ss]
    assays(object)$is_imputed <- dt2mat(data.table::dcast(
                    dt, feature_id ~ sample_id, value.var = 'is_imputed'))
# Plot
    if (plot) print(plot_detects(object, group = !!group))
# Return
    object
}


#==============================================================================
#
#                      cluster_order_features
#                      detect_order_features
#
#==============================================================================

cluster_order_features <- function(object){
    if (nrow(object) < 3) return(object)
    idx <- is.na(exprs(object))
    exprs(object)[idx] <- 0
    order <- hclust(dist(exprs(object)))$order
    exprs(object)[idx] <- NA
    object %<>% extract(order, )
    object
}

detect_order_features <- function(object){
    x <- object
    exprs(x)[is_imputed(x)] <- NA
    idx1 <- fnames(cluster_order_features(x[is_systematic_detect(x),     ]))
    idx2 <- fnames(cluster_order_features(x[is_random_detect(x), ]))
    idx3 <- fnames(cluster_order_features(x[is_full_detect(x),       ]))
    SummarizedExperiment::rbind(object[idx1,], object[idx2,], object[idx3,])
}


#==============================================================================
#
#                           plot_detects
#
#==============================================================================


#' Plot detects
#' @param object   SummarizedExperiment
#' @param group    subgroup
#' @param fill     svar mapped to plot fill
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute=FALSE, plot = FALSE)
#' plot_detects(object)
#' object %<>% impute_systematic_nondetects(plot = FALSE)
#' plot_detects(object)
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, impute = FALSE, plot = FALSE)
#' plot_detects(object)
#' @export
plot_detects <- function(object, group = subgroup, fill = subgroup){
# Process
    detection <- feature_id <- NULL
    set.seed(39)
    group <- enquo(group);         fill     <- enquo(fill)
    groupstr <- as_name(group);    fillstr  <- as_name(fill)
# Reorder samples
    sdata(object)[[groupstr]] %<>% factor()
    object %<>% extract(, order(sdata(.)[[groupstr]]))
# Reorder/block features
    object %<>% detect_order_features()
    y <- object; exprs(y)[is_imputed(y)] <- NA
    nfull       <- sum(is_full_detect(y))
    nsystematic <- sum(is_systematic_detect(y))
    nrandom     <- sum(is_random_detect(y))
# Melt
    plotdt  <-  sumexp_to_long_dt(object)
    alpha <- NULL
    plotdt[,             detection := 'detect']
    plotdt[is.na(value), detection := 'nondetect']
    if ('is_imputed' %in% names(assays(object))){
        plotdt %<>% cbind(is_imputed =
                          sumexp_to_long_dt(object, assay = 'is_imputed')$value)
        plotdt[is_imputed==TRUE, detection := 'impute']
    }
    plotdt[, detection := factor(detection, c('nondetect', 'impute', 'detect'))]
    plotdt[, sample_id  := factor( sample_id, unique(snames(object)))]
    plotdt[, feature_id := factor(feature_id, rev(unique(fnames(object))))]
    colors <- make_colors(unique(sdata(object)[[fillstr]]))
    colors %<>% c(nondetect = "#FFFFFF")
# Plot
    ggplot(plotdt) +
    geom_tile(aes(x=sample_id, y=feature_id, fill=!!fill, alpha=detection)) +
    scale_fill_manual(values = colors) +
    scale_alpha_manual(values = c(nondetect=0, impute=0.3, detect = 1)) +
    ylab('features') +
    xlab('samples') +
    ggtitle(sprintf('detects: %d full, %d random, %d systematic',
                    nfull, nrandom, nsystematic)) +
    theme_bw() +
    theme(  axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
            axis.text.x  = element_blank(), legend.title = element_blank()) +
    geom_hline(yintercept = cumsum(c(nfull, nrandom, nsystematic))) +
    guides(alpha=FALSE)
}


#==============================================================================
#
#                           explore_imputations
#
#==============================================================================


#' Explore imputations
#' @param object SummarizedExperiment
#' @param xbiplot biplot x axis. Default pca1 (symbol)
#' @param ybiplot biplot y axis. Default pca2 (symbol)
#' @param ... aesthetic mappings
#' @return ggplot object
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute = FALSE, plot = FALSE)
#' explore_imputations(object)
#' explore_transformations(object)
#' @export
explore_imputations <- function(object, xbiplot = pca1, ybiplot = pca2, ...){
    imputed     <- impute_systematic_nondetects(object, plot=FALSE)
    zeroed <- impute_systematic_nondetects(object, fun = zeroimpute, plot = FALSE)
    legend <- gglegend(biplot(object))

    do_plot_detects <- function(obj, ...){
        plot_detects(obj, ...) + guides(color=FALSE, fill=FALSE)
    }

    do_biplot <- function(obj, ...){
        biplot(obj, x=!!enquo(xbiplot), y=!!enquo(ybiplot), nloadings=0,...) +
        guides(color=FALSE, fill=FALSE) +
        ggtitle(NULL)
    }

    do_plot_sample_densities <- function(obj, ...){
        plot_sample_densities(obj, ...) +
        guides(color=FALSE, fill=FALSE) +
        ggtitle(NULL)
    }

    p1 <- do_plot_detects(object, ...)  + ggtitle('Original')
    p2 <- do_plot_detects(imputed, ...) + ggtitle('Halfnorm imputed')
    p3 <- do_plot_detects(zeroed, ...)  + ggtitle('Zero imputed')
    p4 <- do_plot_sample_densities(object,  ...)
    p5 <- do_plot_sample_densities(imputed, ...)
    p6 <- do_plot_sample_densities(zeroed,  ...)
    p7 <- do_biplot(object, ...)
    p8 <- do_biplot(imputed, ...)
    p9 <- do_biplot(zeroed, ...)
    grid.arrange(arrangeGrob(p1,p2,p3,p4,p5,p6,p7,p8,p9, nrow=3),
                            legend, ncol=2, widths = c(8, 1))
}


#==============================================================================
#'@title Get/set is_imputed
#'@description Get/Set is_imputed
#'@param object SummarizedExperiment
#'@param value matrix
#'@return matrix (get) or updated object (set)
#'@examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' sum(is_imputed(object))
#' @rdname is_imputed
#' @export
setGeneric("is_imputed",  function(object) standardGeneric("is_imputed") )

#' @rdname is_imputed
setMethod("is_imputed", signature("SummarizedExperiment"),  function(object){
    if ('is_imputed' %in% names(assays(object))){
        assays(object)$is_imputed
    } else {
        matrix(FALSE, nrow = nrow(object), ncol = ncol(object),
                dimnames = dimnames(object))
    }
})

#' @rdname is_imputed
#' @export
setGeneric(
    "is_imputed<-",
    function(object, value)  standardGeneric("is_imputed<-") )

#' @rdname is_imputed
setReplaceMethod(
    "is_imputed",
    signature("SummarizedExperiment", "matrix"),
    function(object, value){
        assays(object)$is_imputed <- value; object})

#' @rdname is_imputed
setReplaceMethod(
    "is_imputed",
    signature("SummarizedExperiment", "NULL"),
    function(object, value){object})

