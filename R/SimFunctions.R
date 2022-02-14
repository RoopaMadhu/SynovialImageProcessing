## From notebooks
##     SampleSeq.ipynb


## Some functions to extract matrices from this list of lists 
probs_from_params <- function(params) {
    .types <- params %>% map('props') %>% map(names) %>% reduce(union)
    mat <- params %>% map('props') %>% map(function(.x) {
        .x <- unlist(.x)
        res <- matrix(0, nrow=1, ncol=length(.types))
        colnames(res) <- .types
        res[, names(.x)] <- prop.table(.x)
        return(res)
    }) %>% 
        reduce(rbind)
    rownames(mat) <- names(params)
    mat
}


ncells_from_params <- function(params) {
    res <- params %>% map('ncells') %>% map(unlist) %>% map(matrix) %>% reduce(rbind)
    rownames(res) <- names(params)
    return(res)
}


## Names are confusing
##    spots: ST spots
##    region: type of spatial region
##    cluster: type of single cell cluster 
## OUTPUT: count matrix of spots by clusters
sample_composition <- function(spots_by_regions, sim_params) {
    ## Extract useful matrices from sim_params json
    regions_by_clusters <- probs_from_params(sim_params)[colnames(spots_by_regions), , drop = FALSE] 
    type_sizes_mat <- ncells_from_params(sim_params)[colnames(spots_by_regions), , drop = FALSE]

    spot_names <- rownames(spots_by_regions)
    cluster_names <- colnames(regions_by_clusters)
#     region_names <- rownames(regions_by_ncells) ## Needed? 

    ## First, get the number of cells total per spot
    ncells_per_spot <- spots_by_regions %*% type_sizes_mat[colnames(spots_by_regions), , drop = FALSE] %>% map_int(rpois, n=1)
    names(ncells_per_spot) <- spot_names

    ## Then, break it down by proportion of cells in that region
    spots_by_clusters <- map2(ncells_per_spot, data.frame(t(spots_by_regions %*% regions_by_clusters)), function(.n, .probs) {
        rmultinom(1, .n, .probs)
    }) %>% 
        purrr::reduce(Matrix::cbind2) %>% 
        t()

    colnames(spots_by_clusters) <- cluster_names
    rownames(spots_by_clusters) <- spot_names
    return(spots_by_clusters)
    
}

## return_type lets you either return a matrix of 
##     "genes": genes x spots
##     "cells": cells x spots 
sample_cells <- function(spots_by_clusters, sc_meta_data, sc_counts, return_type=c('genes', 'cells')[1], cluster_column='Cluster0') {
    
    if (return_type == 'genes' & missing(sc_counts)) {
        stop('To return a count matrix, please provide sc_counts matrix')
    }
    
    ncells <- nrow(sc_meta_data)
    nspots <- nrow(spots_by_clusters)
    sc_meta_data$ROWNUM <- 1:nrow(sc_meta_data) ## TODO: check for conflict with colname
    sc_meta_data$CLUSTER <- sc_meta_data[[cluster_column]]
    
    ## Sample the cells first, then choose from that pool
    .sampled_cells <- map2(colnames(spots_by_clusters), colSums(spots_by_clusters), function(.cluster, .n) {
        sc_meta_data %>%
            subset(CLUSTER == .cluster) %>% 
            dplyr::sample_n(.n, replace = TRUE) %>% 
            with(ROWNUM)
    })
    names(.sampled_cells) <- colnames(spots_by_clusters)
    
    .counters <- rep(0, length(.sampled_cells))
    names(.counters) <- names(.sampled_cells)

    .i <- c()
    .p <- c(0) ## faster to initialize this as fixed-size vector? 

    system.time({
    # profvis::profvis({
        ## Function doesn't actually return anything
        .tmp <- apply(spots_by_clusters, 1, simplify = FALSE, function(.x) {
            if (any(.x > 0)) {
                .x <- .x[.x > 0]
                .i_new <- map2(names(.x), .x, function(.cluster, .n) {
                    ## sample those n cells
                    .res <- .sampled_cells[[.cluster]][(1+.counters[[.cluster]]):(.n+.counters[[.cluster]])]
                    ## remove them from the pool
                    .counters[.cluster] <<- .counters[.cluster] + .n
                    return(.res)
                }) %>% reduce(c)
                ## NOTE: need global assignment b/c inside function
                .i <<- c(.i, .i_new) 
                .p <<- c(.p, length(.i_new) + tail(.p, 1))
            } else {
                .p <<- c(.p, 0 + tail(.p, 1)) ## add an empty column (length = 0)
            }
        })
    })

    cells_by_spots <- Matrix::sparseMatrix(i=.i, p=.p, x=rep(1, length(.i)), dims = c(ncells, nspots))
    colnames(cells_by_spots) <- rownames(spots_by_clusters)
    if (return_type == 'cells') {
        res <- cells_by_spots
    } else if (return_type == 'genes') {
        res <- sc_counts[, 1:nrow(cells_by_spots)] %*% cells_by_spots
    } else {
        stop(glue('return_type \"{return_type}\" not defined'))
    }
    
    return(res)
}

## Epsilon ranges from 0 (no noise at all) to 1 (signal = noise)
sample_noise_global <- function(genes_by_spots, epsilon = 0.05) {
    nspots <- ncol(genes_by_spots)
    ngenes <- nrow(genes_by_spots)
    noise_mean_numi <- log(median(colSums(genes_by_spots)) * epsilon)
    props_bg <- prop.table(Matrix::rowSums(genes_by_spots))

    ## For each spot, decide on number of reads to sample 
    nreads <- round(rlnorm(nspots, meanlog = noise_mean_numi, sdlog = .2))

    .i <- c()
    .p <- c(0)
    .x <- c()

    ## If we could pre-allocated memory for .i and .x, this would be much faster! 
    noise_vecs <- map(nreads, function(.nreads) {
        .res <- rmultinom(.nreads, n=1, prob=props_bg)
        .i_new <- which(.res != 0)
        if (length(.i_new) > 0) {
            .i <<- c(.i, .i_new)
            .x <<- c(.x, .res[.i_new])
            .p <<- c(.p, length(.i_new) + tail(.p, 1))            
        } else {
            ## add empty column (in case of very low noise)
            .p <<- c(.p, tail(.p, 1)) 
        }
    })

    if (tail(.p, 1) == 0) {
        ## Edge case: return empty matrix 
        noise_mat <- Matrix::Matrix(data = 0, nrow = ngenes, ncol = nspots)
    } else {
        noise_mat <- Matrix::sparseMatrix(i = .i, p = .p, x = .x, dims = c(ngenes, nspots))  
    }
    rownames(noise_mat) <- rownames(genes_by_spots)
    colnames(noise_mat) <- colnames(genes_by_spots)
    return(noise_mat)
}



