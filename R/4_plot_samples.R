#==============================================================================
#
#                  add_highlights
#
#==============================================================================
add_highlights <- function(p, hl, geom = geom_point, fixed_color = "black") {
    hlstr <- as_name(hl)
    hl_df <- p$data[get(hlstr)==TRUE]
    args <- list(data = hl_df)
    if (identical(geom, geom_point)) {
        many_hl <- length(unique(args$data$feature_name)) > 6
        if (many_hl) args$data$feature_name <- hlstr
        args %<>% c(
            list(aes(shape = feature_name), size = rel(3), color = fixed_color))
    }
    p <- p + do.call(geom, args)
    if (identical(geom, geom_point)) p <- p +
        labs(shape = if (many_hl) NULL else hlstr) +
        guides(fill = guide_legend(override.aes = list(shape = NA)))
    p
}


#=============================================================================
#
#                 plot_sample_densities()
#                 plot_sample_boxplots()
#                 plot_sample_violins()
#
#=============================================================================

#' Plot sample densities/boxplots/violins
#'
#' @param object      SummarizedExperiment
#' @param fill        svar mapped to fill
#' @param color       svar mapped to color
#' @param group       svar mapped to group
#' @param ...         mapped aesthetics
#' @param highlight   fvar expressing which feature should be highlighted
#' @param fixed       fixed aesthetics
#' @return  ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% extract(, order(.$subgroup))
#' plot_sample_densities(object)
#' plot_sample_boxplots(object)
#' plot_sample_violins(object)
#'
#' object %<>% extract(, 1:10)
#' plot_sample_densities(object)
#' plot_sample_boxplots(object)
#' plot_sample_violins(object)
#'
#' fdata(object) %<>% cbind(control=.$feature_name %in% c('biotin','phosphate'))
#' plot_sample_boxplots(object, highlight = control)
#' plot_sample_violins(object,  highlight = control)
#' @export
plot_sample_densities <- function(object, fill = subgroup, color = NULL,
    group = sample_id, ..., fixed = list(alpha = 0.5, na.rm = TRUE)
){
    dt <- sumexp_to_long_dt(object, svars = svars(object))
    fill  <- enquo(fill)
    color <- enquo(color)
    group <- enquo(group)
    value <- NULL
    plot_data(dt, geom = geom_density, x = value, fill = !!fill,
        color = !!color, group = !!group, ..., fixed = fixed) +
    ggtitle('sample densities')
}


#' @rdname plot_sample_densities
#' @export
plot_sample_violins <- function(object, fill = subgroup, color = NULL, ...,
    highlight = NULL, fixed = list(na.rm=TRUE)
){
# Process
    sample_id <- value <- NULL
    fill      <- enquo(fill)
    color     <- enquo(color)
    highlight <- enquo(highlight)
    fillstr      <- if (quo_is_null(fill))      NULL else as_name(fill)
    highlightstr <- if (quo_is_null(highlight)) NULL else as_name(highlight)
# Prepare
    show_highlights <- !is.null(highlightstr) &&
                        is.element(highlightstr, fvars(object)) &&
                        any(fdata(object)[[highlightstr]], na.rm = TRUE)
    fvars <- character(0)
    if (show_highlights) fvars %<>% c(highlightstr, "feature_name")
    dt <- sumexp_to_long_dt(object, svars = svars(object), fvars = fvars)
# Plot
    p <- plot_data(dt, geom = geom_violin, x = sample_id, y = value,
                   fill = !!fill, color= !!color, ..., fixed = fixed)
    if (show_highlights) p %<>% add_highlights(highlight, geom = geom_point)
# Finish
    breaks <- levels(dt$sample_id)
    if (length(breaks)>50) breaks <- dt[, .SD[1], by = fillstr][, sample_id]
    p <- p + xlab(NULL) + scale_x_discrete(breaks = breaks) +
        theme(axis.text.x = element_text(angle=90, hjust=1))
# Return
    p
}


#' @rdname plot_sample_densities
#' @export
plot_sample_boxplots <- function(object, fill = subgroup, color = NULL, ...,
    highlight = NULL, fixed = list(na.rm=TRUE)
){
# Process
    sample_id <- value <- NULL
    fill      <- enquo(fill)
    color     <- enquo(color)
    highlight <- enquo(highlight)
    fillstr      <- if (quo_is_null(fill))      NULL else as_name(fill)
    highlightstr <- if (quo_is_null(highlight)) NULL else as_name(highlight)
# Prepare
    show_highlights <- !is.null(highlightstr) &&
                        is.element(highlightstr, fvars(object)) &&
                        any(fdata(object)[[highlightstr]], na.rm = TRUE)
    fvars <- character(0)
    if (show_highlights) fvars %<>% c(highlightstr, "feature_name")
    dt <- sumexp_to_long_dt(object, svars = svars(object), fvars = fvars)
# Plot
    p <- plot_data(dt, geom = geom_boxplot, x = sample_id, y = value,
                   fill = !!fill, color = !!color, ..., fixed = fixed)
    if (show_highlights) p %<>% add_highlights( highlight, geom = geom_point,
                                                fixed_color = "darkred")
# Finish
    breaks <- levels(dt$sample_id)
    if (length(breaks)>50) breaks <- dt[, .SD[1], by = fillstr][, sample_id]
    p <- p + xlab(NULL) + scale_x_discrete(breaks = breaks) +
        theme(axis.text.x = element_text(angle=90, hjust=1))
# Return
    p
}


#=============================================================================
#
#                       biplot()
#
#==============================================================================


#' Biplot
#' @param object         SummarizedExperiment
#' @param x              pca1, etc.
#' @param y              pca2, etc.
#' @param color          svar mapped to color (symbol)
#' @param label          svar mapped to label (symbol)
#' @param ...            additional svars mapped to aesthetics
#' @param feature_label  fvar mapped to (loadings) label
#' @param fixed          fixed plot aesthetics
#' @param nloadings      number of loadings per half-axis to plot
#' @return ggplot object
#' @examples
#' # Read data
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot = FALSE)
#' # biplot.SummarizedExperiment
#'     biplot(object)
#'     biplot(object, x=pca1, y=pca2)
#'     biplot(object, x=pls1, y=pls2)
#'     biplot(object, x=lda1, y=lda2)
#'     biplot(object, x=pca3, y=pca4)
#'     biplot(object, nloadings = 1)
#'     biplot(object, color = TIME_POINT)
#'     biplot(object, color = NULL)
#'# biplot.matrix
#'     object %<>% exprs()
#' @export
setGeneric('biplot',
function(object, ...)  standardGeneric("biplot"), signature = "object")

add_scores <- function(
    p, sdata, x = pca1, y = pca2, color = subgroup, ...,
    fixed = list(shape=15, size=3)
){
    x     <- enquo(x)
    y     <- enquo(y)
    color <- enquo(color)

    p + layer(  geom = 'point',
                mapping = aes(x = !!x, y = !!y, color = !!color, ...),
                stat    = "identity",
                data    = sdata,
                params  = fixed,
                position= 'identity')

}

add_loadings <- function(
    p, fdata, sdata, x = pca1, y = pca2, label = feature_name, nloadings = 1
){
# Process args
    if (nloadings==0) return(p)
    x     <- enquo(x);  xstr <- rlang::as_name(x)
    y     <- enquo(y);  ystr <- rlang::as_name(y)
    label <- enquo(label)
# Loadings
    xloadings <- fdata[[xstr]]
    yloadings <- fdata[[ystr]]
    idx <- unique(c(headtail(order(xloadings, na.last=NA), nloadings),
                    headtail(order(yloadings, na.last=NA), nloadings)))
# Scale loadings to scoreplot
    xscores <- sdata[[xstr]]
    yscores <- sdata[[ystr]]
    maxscore <- min(abs(min(c(xscores, yscores, na.rm=TRUE))),
                    abs(max(c(xscores, yscores, na.rm=TRUE))), na.rm=TRUE)
    scorefactor <- maxscore/max(abs(c(xloadings, yloadings)),  na.rm=TRUE)
    plotdt <- fdata
    plotdt[[xstr]] %<>% multiply_by(scorefactor)
    plotdt[[ystr]] %<>% multiply_by(scorefactor)
    plotdt %<>% extract(idx, )
# Plot
    feature_name <- NULL
    if (!'feature_name' %in% names(plotdt)){
        setnames(plotdt, 'feature_id', 'feature_name')}
    p + layer(  geom     = 'segment',
                mapping  = aes(x=0, y=0, xend=!!x, yend=!!y),
                stat     = "identity",
                data     = plotdt, # list(alpha = 0.05, size=3),
                params   = list(alpha = 0.1, size=1, na.rm = TRUE),
                position = "identity") +
        layer(  geom     = "text",
                mapping  = aes(x = !!x, y = !!y, label = !!label),
                stat     = "identity",
                data     = plotdt,
                params   = list(alpha = 0.5, na.rm = TRUE),
                position ='identity')
}


#' @rdname biplot
#' @export
setMethod("biplot", signature("matrix"),
function(object, x=pca1, y=pca2, color=subgroup, group=subgroup, label = NULL,
         feature_label = feature_name, ..., sampledata, featuredata,
         fixed = list(shape=15, size=3), nloadings = 0
){
    x     <- enquo(x);     xstr <- as_name(x)
    y     <- enquo(y);     ystr <- as_name(y)
    group <- enquo(group)
    color <- enquo(color)
    label <- enquo(label)
    dots  <- enquos(...)
    fixed %<>% extract(setdiff(names(fixed), names(dots)))
    feature_label <- enquo(feature_label)
    methodx <- substr(xstr, 1, 3)
    methody <- substr(ystr, 1, 3)
    xdim <- xstr %>% substr(4, nchar(.)) %>% as.numeric()
    ydim <- ystr %>% substr(4, nchar(.)) %>% as.numeric()
    groups <- rlang::eval_tidy(group, data=sampledata)
    outx <- get(methodx)(object, groups, ndim=xdim, verbose = FALSE)
    outy <- get(methody)(object, groups, ndim=ydim, verbose = FALSE)

    sampledata  %<>% data.table()
    featuredata %<>% data.table()
    sampledata[, (xstr) := outx$samples[, xstr]]
    sampledata[, (ystr) := outy$samples[, ystr]]
    featuredata[,(xstr) := outx$features[,xstr]]
    featuredata[,(ystr) := outy$features[,ystr]]

    xlab  <- paste0(xstr, ' : ', outx$variances[[xstr]],'% ')
    ylab  <- paste0(ystr, ' : ', outy$variances[[ystr]],'% ')

    p <- ggplot() + theme_bw() + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
    p <- p + ggtitle(paste0(unique(c(methodx, methody)), collapse = '/'))
    p %<>% add_scores(sampledata, x=!!x, y=!!y, color = !!color, !!!dots,
                      fixed = fixed)
    p %<>% add_loadings(fdata=featuredata, sdata=sampledata, x = !!x, y = !!y,
                        label = !!feature_label, nloadings = nloadings)
    if (!quo_is_null(label)) p <- p + geom_text_repel(
            aes(x=!!x, y=!!y, label=!!label), data=sampledata)
    p %<>% add_color_scale(!!color, sampledata)
    p

})

#' @rdname biplot
#' @export
setMethod("biplot", signature("SummarizedExperiment"),
function(object, x=pca1, y=pca2, group=subgroup, color=subgroup, label = NULL,
    feature_label = feature_name, ...,
    fixed = list(shape=15, size=3), nloadings = 0
){
    biplot(exprs(object), x=!!enquo(x), y=!!enquo(y), group=!!enquo(group),
           color=!!enquo(color), label=!!enquo(label),
           feature_label=!!enquo(feature_label), ...,
           sampledata = sdata(object), featuredata = fdata(object),
           fixed = fixed)
})





headtail <- function(x, n){
    c(x[seq(1, n)], x[seq(length(x)+1-n, length(x))])
}

.filter_minvar <- function(object, method, minvar) {
    variances <- metadata(object)[[method]]
    discard_components <- variances[variances < minvar] %>% names()

    sdata(object)[discard_components] <- NULL
    fdata(object)[discard_components] <- NULL
    metadata(object)[[method]] <-
        variances[!names(variances) %in% discard_components]
    object
}

#=============================================================================
#
#                       biplot_corrections()
#                       biplot_covariates()
#
#==============================================================================

#' @export
#' @rdname biplot_corrections
plot_corrections <- function(...){
    .Deprecated("biplot_corrections")
    biplot_corrections(...)
}

#' Biplot batch corrections
#'
#' @param object      SummarizedExperiment
#' @param method      'pca', 'pls', 'lda', or 'sma'
#' @param color       variable mapped to color (symbol)
#' @param covariates  covariates to be batch-corrected
#' @param varcols     number of covariate columns
#' @param plot        TRUE/FALSE: plot?
#' @param ...         used to maintain deprecated functions
#' @return  grid object
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' biplot_corrections(object,  covariates = c('SEX', 'T2D', 'SUB', 'SET'))
#' @seealso biplot_covariates
#' @export
biplot_corrections <- function(
    object, method = 'pca', color = subgroup, covariates = character(0),
    varcols = ceiling(sqrt(1+length(covariates))), plot = TRUE
){
    p <- biplot(object, pca1, pca2, color = !!enquo(color), nloadings=0)
    p <- p + ggtitle('INPUT')
    p <- p + guides(color=FALSE, fill=FALSE)
    plotlist <- list(p)
    for (ibatch in covariates){
        exprs(object) %<>% removeBatchEffect(batch=sdata(object)[[ibatch]])
        p <- biplot(object, pca1, pca2, color = !!enquo(color), nloadings=0)
        p <- p + ggtitle(paste0(' - ', ibatch))
        p <- p + guides(color=FALSE, fill=FALSE)
        plotlist %<>% c(list(p))
    }
    pp <- arrangeGrob(grobs = plotlist, ncol = varcols)
    if (plot) grid::grid.draw(pp)
    invisible(pp)
}


#' @rdname biplot_covariates
#' @export
plot_covariates <- function(...){
    .Deprecated('biplot_covariates')
    biplot_covariates(...)
}


#' Biplot covariates
#'
#' @param object     SummarizedExperiment
#' @param method     'pca', 'pls', 'lda', or 'sma'
#' @param covariates  covariates: mapped to color or batch-corrected
#' @param ndim        number of dimensions to plot
#' @param dimcols     number of dimension columns
#' @param varcols     number of covariate columns
#' @param plot        TRUE or FALSE: whether to plot
#' @param ...         used to maintain deprecated functions
#' @return  ggplot object
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' biplot_covariates(object, ndim = 12, dimcols = 3)
#' biplot_covariates(object, covariates = c('SEX', 'T2D', 'SUB', 'SET'))
#' biplot_covariates(object, covariates = c('SEX', 'T2D', 'SUB', 'SET'), ndim=2)
#' biplot_covariates(object, covariates = c('subgroup'), dimcols = 3)
#' @seealso biplot_corrections
#' @export
biplot_covariates <- function(
    object, method = 'pca', covariates = 'subgroup', ndim = 6,
    dimcols = 1, varcols = length(covariates), plot = TRUE
){
    x <- y <- NULL
    plotdt <- prep_covariates(object, method = 'pca', ndim=ndim)
    plotlist <- list()
    for (covar in covariates){
        p <- plot_data(plotdt, geom = geom_point, x=x, y=y, color=!!sym(covar),
                        fixed = list(shape=15, size=3))
        p <- p + facet_wrap(~dims, ncol = dimcols, scales = 'free')
        p <- p + xlab(NULL) + ylab(NULL) + ggtitle(covar)
        p <- p + theme(legend.position = 'bottom', legend.title=element_blank())
        plotlist %<>% c(list(p))
    }
    pp <- gridExtra::arrangeGrob(grobs = plotlist, ncol = varcols)
    if (plot) grid::grid.draw(pp)
    invisible(pp)
}

prep_covariates <- function(object, method='pca', ndim=6){

    plotdt <- cbind(sdata(object)[FALSE,], x= character(0), y = character(0))
    projdt <- data.table(sdata(get(method)(object, ndim=ndim)))
    alldims <- names(projdt) %>% extract(stri_detect_fixed(., method)) %>%
                stri_replace_first_fixed(method, '') %>% as.numeric()
    ndim <- min(c(max(alldims), ndim))
    npairs <- ndim %/% 2
    for (idim in seq_len(npairs)){
        dim1 <- idim*2-1
        dim2 <- idim*2
        xvar <- paste0(method, dim1)
        yvar <- paste0(method, dim2)
        tmpdt <- data.table::copy(projdt)
        setnames(tmpdt, c(xvar, yvar), c('x', 'y'))
        tmpdt %<>% extract(, stri_detect_fixed(
                                names(.), method, negate = TRUE), with = FALSE)
        tmpdt$dims <- paste0(dim1, ':', dim2)
        plotdt %<>% rbind(tmpdt)
    }
    plotdt$dims %<>% factor(unique(.))
    plotdt
}


#' Plot samples
#'
#' Plots sample densities and scores
#' @param object    SummarizedExperiment
#' @param x         svar mapped to biplot x (sym, default pca1)
#' @param y         svar mapped to biplot y *sym, default pca2)
#' @param color     svar mapped to biplot color and density fill
#' @param nloadings n loadings added to plot
#' @param ...       passed to plot_data
#' @param plot      whether to print plot
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' plot_samples(object)
#'
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot = FALSE)
#' plot_samples(object)
#' @export
plot_samples <- function(
    object, x=pca1, y=pca2, color=subgroup, nloadings=0, ..., plot=TRUE
){
    p1  <-  biplot(object, x=!!enquo(x), y=!!enquo(y), color=!!enquo(color),
                nloadings = nloadings, ...) +
            theme(legend.position='top')
    p2  <-  plot_sample_densities(object, fill = !!enquo(color), ...) +
            theme(legend.position='none')
    p3 <- plot_detects(object) + theme(legend.position='none')
    p4 <- gglegend(p1)
    p1 <- p1 + theme(legend.position='none')
    pp <- grid.arrange( p4,
                        arrangeGrob(grobs = list(p1, p2, p3), ncol = 3),
                        ncol=1, heights = c(2,8))
    if (plot) grid.draw(pp)
    invisible(pp)
}




