#=============================================================================
#
#                  flip_sign_if_all_exprs_are_negative.
#                  evenify_upwards
#
#=============================================================================

#' Flip sign if all expr values are negative
#' @param x matrix
#' @param verbose TRUE (default) or FALSE
#' @return updated matrix
#' @noRd
flip_sign_if_all_exprs_are_negative <- function(x, verbose=TRUE){
    idx <- !is.na(x)
    if (all(sign(x[idx])==-1)){
        if (verbose) cmessage(
            '\t\tAll values negative: flip signs to prevent singularities.')
        x %<>% multiply_by(-1)
    }
    x
}



#' Is even/odd?
#' @param x integer
#' @return TRUE or FALSE
#' @examples
#' is_even(13)
#' is_even(12)
#' is_odd(13)
#' is_odd(12)
#' @noRd
is_even <- function(x)   (x %% 2) == 0
is_odd  <- function(x)    !is_even(x)


#' Has even/odd length?
#' @param x vector
#' @return logical
#' @examples
#' has_even_length(1:2)
#' has_odd_length(1:2)
#' has_even_length(1:3)
#' has_odd_length(1:3)
#' @noRd
has_even_length <- function(x)   is_even(length(x))
has_odd_length  <- function(x)   is_odd(length(x))


#' Evenify upwards
#' @param x integer
#' @return integer
#' @examples
#' evenify_upwards(3)
#' @noRd
evenify_upwards <- function(x)   if (is_odd(x)) x+1 else x



#=============================================================================
#
#                               merge_sdata
#                               merge_fdata
#
#=============================================================================

merge_fdata <- function(object, df, by = 'feature_id'){
    df %<>% as.data.frame() # convert matrix
    if (!'feature_id' %in% names(df))  df$feature_id <- rownames(df)
    duplicate_cols <- setdiff(intersect(fvars(object), names(df)), 'feature_id')
    fdata(object)[duplicate_cols] <- NULL
    fdata(object) %<>% merge(df, by = by, all.x = TRUE, sort = FALSE)
    fnames(object) <- fdata(object)$feature_id # merging drops them!
    object
}


#'@examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% merge_sdata( data.frame(sample_id = object$sample_id,
#'                                     number = seq_along(object$sample_id)))
#' sdata(object)
#' @noRd
merge_sdata <- function(object, df, by = 'sample_id'){
    df %<>% as.data.frame() # convert matrix to df
    if (!'sample_id' %in% names(df))  df$sample_id <- rownames(df)
    duplicate_cols <- setdiff(intersect(svars(object), names(df)), 'sample_id')
    sdata(object)[duplicate_cols] <- NULL
    sdata(object) %<>% merge(df, by = by, all.x = TRUE, sort = FALSE)
    snames(object) <- object$sample_id # merging drops them!
    if ('subgroup'  %in% svars(object)) object$subgroup  %<>% as.character()
    if ('replicate' %in% svars(object)) object$replicate %<>% as.character()
    object
}


#============================================================================
#
#                    pca, sma, lda, pls, spls, ropls
#
#============================================================================

#' Add PCA, SMA, LDA, or PLS
#'
#' Perform a dimension reduction.
#' Add sample scores, feature loadings, and dimension variances to object.
#'
#' @param object  SummarizedExperiment/Matrix
#' @param ndim    number
#' @param minvar  number
#' @param verbose TRUE (default) or FALSE
#' @return        SummarizedExperiment
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#'
#' pca(exprs(object))
#' sdata(pca(object))
#'
#' pls(object)
#' pls(exprs(object), object$subgroup)
#'
#' lda(object)
#' lda(exprs(object), object$subgroup)
#'
#' sma(object)
#' sma(exprs(object))
#'
#' pca(object, ndim=3)
#' pca(object, ndim=Inf, minvar=5)
#' @author Aditya Bhagwat, Laure Cougnaud (LDA)
#' @export
setGeneric('pca',
            function(object, ...)  standardGeneric("pca"), signature = "object")

#' @rdname pca
#' @export
setGeneric('pls',
           function(object, ...)  standardGeneric("pls"),  signature = "object")

#' @rdname pca
#' @export
setGeneric('lda',
           function(object, ...)  standardGeneric("lda"),  signature = "object")

#' @rdname pca
#' @export
setGeneric('sma',
           function(object, ...)  standardGeneric("sma"),  signature = "object")

#' @rdname pca
#' @export
setMethod("pca", signature("matrix"),
function(object, group = NULL, ndim=2, verbose=TRUE){
# Assert
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_less_than_or_equal_to(ndim, ncol(object))
# Prepare
    object %<>% inf_to_na(verbose=verbose)
    object %<>% nan_to_na(verbose=verbose)
    row_means <- rowMeans(object, na.rm=TRUE)
    col_means <- colWeightedMeans(object, abs(row_means), na.rm = TRUE)
    global_mean <- mean(col_means)
    object %<>% apply(1, '-', col_means)  %>%  # Center columns
                apply(1, '-', row_means)  %>%  # Center rows
                add(global_mean)          %>%  # Add doubly subtracted
                divide_by(sd(., na.rm=TRUE))   # Normalize
# Pca
    pca_res  <- pcaMethods::pca(t(object), nPcs = ndim, scale = 'none',
                                center = FALSE, method = 'nipals')
    samples   <- pca_res@scores
    features  <- pca_res@loadings
    variances <- round(100*pca_res@R2)
    colnames(samples)  <- sprintf('pca%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('pca%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('pca%d', seq_len(length(variances)))
    list(samples=samples, features=features, variances=variances)
})


#' @rdname pca
#' @export
setMethod("pca", signature("SummarizedExperiment"),
function(object, group = NULL, ndim=2, minvar=0, verbose=TRUE){
    # Assert
        assert_is_valid_sumexp(object)
        assert_is_a_number(minvar)
        assert_all_are_in_range(minvar, 0, 100)
        . <- NULL
        if (verbose)  message('\tAdd PCA')
    # Add pca
        object %<>% rm_missing_in_all_samples(verbose = verbose)
        out <- pca(exprs(object), ndim=ndim, verbose=verbose)
        object %<>% merge_sdata(out$samples)
        object %<>% merge_fdata(out$features)
        metadata(object)$pca <- out$variances
    # Filter for minvar
        object %<>% .filter_minvar('pca', minvar)
    # Return
        object
})


#' @rdname pca
#' @export
setMethod("sma", signature("matrix"),
function(object, group = NULL, ndim=NULL, verbose=TRUE){
# Assert
    if (!requireNamespace('mpm', quietly = TRUE)){
        message("First Biocinstaller::install('mpm'). Then re-run.")
        return(object)
    }
# Prepare
    object %<>% minusinf_to_na(verbose = verbose)   # else SVD singular
    object %<>% flip_sign_if_all_exprs_are_negative(verbose = verbose)
    object %<>% extract(rowAlls(!is.na(.)), )
# Transform
    dt <- data.table(feature = rownames(object)) %>% cbind(object)
    mpm_tmp <- mpm::mpm(
        dt, logtrans = FALSE, closure = 'none', center = 'double',
        normal = 'global', row.weight = 'mean', col.weight = 'constant')
    ncomponents <- length(mpm_tmp$contrib)
    mpm_out <- mpm::plot.mpm(mpm_tmp, do.plot=FALSE, dim = seq_len(ncomponents))
# Extract
    samples   <- mpm_out$Columns
    features  <- mpm_out$Rows
    variances <- round(100*mpm_tmp$contrib[seq_len(ncomponents)])
    names(samples)   <- sprintf('sma%d', seq_len(ncol(samples)))
    names(features)  <- sprintf('sma%d', seq_len(ncol(features)))
    names(variances) <- sprintf('sma%d', seq_len(length(variances)))
# Return
    list(samples = samples, features = features, variances = variances)
})


#' @rdname pca
#' @export
setMethod("sma", signature("SummarizedExperiment"),
function(object, ndim=2, minvar=0, verbose=TRUE){
# Assert
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Transform/Restrict/Merge
    out <- sma(exprs(object), verbose=verbose)
    if (is.infinite(ndim)) ndim <- ncol(samples)
    out$samples   %<>% extract(, seq_len(ndim), drop = FALSE)
    out$features  %<>% extract(, seq_len(ndim), drop = FALSE)
    out$variances %<>% extract(  seq_len(ndim))
    object %<>% merge_sdata(out$samples)
    object %<>% merge_fdata(out$features)
    metadata(object)$sma <- out$variances
# Filter for minvar
    object %<>% .filter_minvar('sma', minvar)
# Return
    object
})

#' @rdname pca
#' @export
setMethod("lda", signature("matrix"),
function(object, group, ndim=NULL, verbose=TRUE){
# Prepare
    object %<>% minusinf_to_na(verbose = verbose)         # SVD singular
    object %<>% flip_sign_if_all_exprs_are_negative(verbose = verbose)
    object %<>% extract(rowAlls(!is.na(.)), )
# Transform
    exprs_t  <- t(object)
    lda_out  <- suppressWarnings(MASS::lda(exprs_t, grouping = group))
    features <- lda_out$scaling
    if (ncol(features)==1) features %<>% cbind(LD2 = 0)
    exprs_t %<>% scale(center = colMeans(lda_out$means), scale = FALSE)
    samples  <- exprs_t %*% features
    variances <- round((lda_out$svd^2)/sum(lda_out$svd^2)*100)
    if (length(variances)==1) variances <- c(LD1 = variances, LD2 = 0)
# Rename
    colnames(samples)  <- sprintf('lda%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('lda%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('lda%d', seq_len(length(variances)))
# Return
    list(samples=samples, features=features, variances=variances)
})


#' @rdname pca
#' @export
setMethod("lda", signature("SummarizedExperiment"),
function(object, ndim=2, minvar=0, verbose=TRUE){
# Assert
    assert_is_valid_sumexp(object)
    nsubgroup <- length(subgroup_levels(object))
    if (is.infinite(ndim))  ndim <- nsubgroup - 1
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, nsubgroup-1)
    if (ndim > (nsubgroup-1)) stop(
        sprintf('LDA requires ndim (%d) <= nsubgroup-1 (%d)',ndim, nsubgroup-1))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Transform
    groups <- sdata(object)$subgroup
    out <- lda(exprs(object), groups=groups, verbose=verbose)
    out$samples   %<>% extract(, seq_len(ndim), drop = FALSE)
    out$features  %<>% extract(, seq_len(ndim), drop = FALSE)
    out$variances %<>% extract(  seq_len(ndim))
# Merge
    object %<>% merge_sdata(out$samples)
    object %<>% merge_fdata(out$features)
    metadata(object)$sma <- out$variances
# Filter for minvar
    object %<>% .filter_minvar('lda', minvar)
# Return
    object
})


#' @rdname pca
#' @export
setMethod("pls", signature("matrix"),
function(object, group, ndim=2, verbose=FALSE){
    pls_out <- mixOmics::plsda(t(object), group, ncomp = ndim)
    samples   <- pls_out$variates$X
    features  <- pls_out$loadings$X
    variances <- round(100*pls_out$explained_variance$X)
    colnames(samples)  <- sprintf('pls%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('pls%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('pls%d', seq_len(length(variances)))
    list(samples=samples, features=features, variances=variances)
})


#' @rdname pca
#' @export
setMethod("pls", signature("SummarizedExperiment"),
function(object, ndim=2, minvar=0, verbose=FALSE){
# Assert
    if (!requireNamespace('mixOmics', quietly = TRUE)){
        stop("BiocManager::install('mixOmics'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Transform
    groups <- subgroup_values(object)
    out <- pls(exprs(object), groups=groups, ndim=ndim, verbose=verbose)
# Add/Filter/Return
    object %<>% merge_sdata(out$samples)
    object %<>% merge_fdata(out$features)
    metadata(object)$pls <- out$variances
    object %<>% .filter_minvar('pls', minvar)
    object
})


#' @param object  SummarizedExperiment
#' @param ndim    number
#' @return        SummarizedExperiment
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' spls(object)
#' @noRd
spls <- function(object, ndim = 2, minvar = 0){
# Assert
    if (!requireNamespace('mixOmics', quietly = TRUE)){
        stop("BiocManager::install('mixOmics'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Transform
    x <- t(exprs(object))
    y <- subgroup_values(object)
    pls_out <- mixOmics::splsda( x, y, ncomp = ndim)
    samples   <- pls_out$variates$X
    features  <- pls_out$loadings$X
    variances <- round(100*pls_out$explained_variance$X)
    colnames(samples)  <- sprintf('pls%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('pls%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('pls%d', seq_len(length(variances)))
# Add
    object %<>% merge_sdata(samples)
    object %<>% merge_sdata(features)
    metadata(object)$spls <- variances
# Filter for minvar
    object %<>% .filter_minvar('spls', minvar)
# Return
    object
}


#' @param object  SummarizedExperiment
#' @param ndim    number
#' @return        SummarizedExperiment
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' opls(object)
#' @noRd
opls <- function(object, ndim = 2, minvar = 0){
# Assert
    if (!requireNamespace('ropls', quietly = TRUE)){
        message("BiocManager::install('ropls'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Transform
    x <- t(exprs(object))
    y <- subgroup_values(object)
    pls_out <- ropls::opls(x, y, predI = ndim, permI = 0, fig.pdfC = FALSE)
    samples   <- pls_out@scoreMN
    features  <- pls_out@loadingMN
    variances <- round(pls_out@modelDF$R2X*100)
    colnames(samples)  <- sprintf('pls%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('pls%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('pls%d', seq_len(length(variances)))
# Add
    object %<>% merge_sdata(samples)
    object %<>% merge_fdata(features)
    metadata(object)$opls <- variances
# Filter for minvar
    object %<>% .filter_minvar('opls', minvar)
# Return
    object
}


