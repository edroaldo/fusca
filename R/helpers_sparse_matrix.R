#' Apply for sparse matrix.
#'
#' Sparse matrix version of the apply function. Avoid errors when working with
#' large matrix. Apply the function on each element or row.
#'
#' @param X matrix; data.
#' @param FUN function.
#'
#' @return res; the result depends on the function. If the function aggregates,
#' then it will do it by rows and return a vector. If the function does not
#' aggreate, it will return a list with each row.
# apply1_sp <- function(X, FUN) {
#   res <- numeric(nrow(X))
#   # The dgTMatrix is the sparse matrix in the triplet form: x are the non zero
#   # values, i are the row indexes, and j are the column indexes.
#   X2 <- as(X, "dgTMatrix")
#   # Apply to rows (i).
#   tmp <- tapply(X2@x, X2@i, FUN)
#   res[as.integer(names(tmp)) + 1] <- tmp
#   res
# }
apply1_sp <- function(X, FUN) {
  res <- numeric(nrow(X))
  X2 <- as(X, "dgTMatrix")
  tmp <- tapply(X2@x, X2@i, FUN)
  res[as.integer(names(tmp)) + 1] <- tmp
  res
}

#' Apply for sparse matrix.
#'
#' Sparse matrix version of the apply function. Avoid errors when working with
#' large matrix. Apply the function on each element or column. It is the same as
#' the apply1_sp, but for columns. The matrix is transposed.
#'
#' @param X matrix; data.
#' @param FUN function.
#'
#' @return res; the result depends on the function. If the function aggregates,
#' then it will do it by rows and return a vector. If the function does not
#' aggreate, it will return a list with each row.
apply2_sp <- function(X, FUN) {
  res <- numeric(nrow(X))
  # The dgTMatrix is the sparse matrix in the triplet form: x are the non zero
  # values, i are the row indexes, and j are the column indexes.
  X2 <- t(as(X, "dgTMatrix"))
  # Apply to rows (i).
  tmp <- tapply(X2@x, X2@i, FUN)
  res[as.integer(names(tmp)) + 1] <- tmp
  res
}

#' Correlation function for sparse matrix.
#'
#' Code found in \href{https://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r}{Running cor() (or any variant) over a sparse matrix in R}.
#'
#' @param x matrix; data.
sparse.cor <- function(x){
  n <- nrow(x)
  covmat <- (Matrix::crossprod(x)-2*(Matrix::colMeans(x) %o% Matrix::colSums(x))
             +n*Matrix::colMeans(x) %o% Matrix::colMeans(x))/(n-1)
  # Standard deviations of columns.
  sdvec <- sqrt(Matrix::diag(covmat))
  # Correlation matrix.
  covmat/crossprod(Matrix::t(sdvec))
}

#' Correlation function for sparse matrix.
#'
#' Alternative version of sparse.cor1.
#'
#' @param x matrix; data.
sparse.cor2 <- function(x){
  n <- nrow(x)
  covmat <- (Matrix::crossprod(x)-2*(Matrix::colMeans(x) %o% Matrix::colSums(x))
             +n*Matrix::colMeans(x) %o% Matrix::colMeans(x))/(n-1)
  sdvec <- sqrt(Matrix::diag(covmat)) # standard deviations of columns
  covmat/Matrix::crossprod(Matrix::t(sdvec)) # correlation matrix
}


#' Scale sparse matrix dividing it in blocks.
#'
#' An alternative version of scale for sparse matrix. Avoid errors for large
#' matrixes.
#'
#' @param matrix any matix; the data tha will be scaled.
#' @param blocksize numeric; the size of the blocks used in each scaling.
#'
#' @return any matrix; the scaled matrix.
scale.sparse <- function(matrix, blocksize=500){
  chunks <- split(colnames(matrix),
                  ceiling(seq_along(colnames(matrix))/blocksize))
  tmp.data <- lapply(chunks, function(x){
    scale(matrix[,x])
  })
  scaled.matrix <- as(do.call(cbind, tmp.data), 'dgCMatrix')
  return(scaled.matrix)
}
