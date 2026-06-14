#' Calculate Parallel Distance Matrix via Fork-Based mclapply
#'
#' @param mat A matrix where rows are genomic bins and columns are cells.
#' @param method Character string: "emd", "euclidean", "pearson", or "spearman".
#' @param num_cores Integer specifying CPU cores. Defaults to total available cores minus 1.
#' @return A formal R 'dist' class object.
#'
#' @export

cal_dist <- function(mat, method="emd", num_cores=NULL){
    method <- tolower(method)
    n_cols <- ncol(mat)
    col_names <- colnames(mat)

    # Determine optimal core count
    if (is.null(num_cores)) {
      num_cores <- max(1, parallel::detectCores() - 2)
    }

    # Generate unique pairwise column combinations
    pairs <- combn(n_cols, 2)
    n_pairs <- ncol(pairs)

    # Execute pairwise comparisons using mclapply (fork-based)
    # Instead of indices, we pass the column index of the 'pairs' matrix

    dist_list <- parallel::mclapply(1:n_pairs, function(idx) {
      i <- pairs[1, idx]
      j <- pairs[2, idx]
      vec_i <- mat[, i]
      vec_j <- mat[, j]

      if (method == "emd") {
        return(transport::wasserstein1d(vec_i, vec_j))

      } else if (method == "euclidean") {
        return(sqrt(sum((vec_i - vec_j)^2)))

      } else if (method %in% c("pearson", "spearman")) {
        r <- cor(vec_i, vec_j, method = method, use = "pairwise.complete.obs")
        return(1 - r)
      }
    }, mc.cores = num_cores)

    # Convert list output back to a vector

    dist_values <- unlist(dist_list)

    # Reconstruct the symmetric square matrix
    dist_mat <- matrix(0, nrow = n_cols, ncol = n_cols)
    rownames(dist_mat) <- col_names
    colnames(dist_mat) <- col_names

    for (idx in 1:n_pairs) {
      i <- pairs[1, idx]
      j <- pairs[2, idx]
      dist_mat[i, j] <- dist_values[idx]
      dist_mat[j, i] <- dist_values[idx] # Mirror lower triangle
    }

    return(as.dist(dist_mat))
}





