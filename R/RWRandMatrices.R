## R Code for the Random Walk with Restart Package (RandomWalkRestartMH).

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Functions to compute the matrices and perform the Random Walks.
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## Roxy Documentation comments
#' Computes the adjacency matrix of a multiplex network
#'
#' \code{compute.adjacency.matrix} is a function to compute the adjacency
#' matrix of a multiplex network provided as a \code{Multiplex} object.
#'
#' @usage compute.adjacency.matrix(x,delta = 0.5)
#'
#' @details The parameter \code{delta} sets the probability to change between
#' layers at the next step. If delta = 0, the particle will always remain
#' in the same layer after a non-restart iteration. On the other hand, if
#' delta = 1, the particle will always change between layers, therefore
#' not following the specific edges of each layer.
#'
#' @param x A \code{Multiplex} object describing a multiplex network generated
#' by the function \code{create.multiplex}.
#' @param delta A numeric value between 0 and 1. It sets the probability
#' of performing inter-layer versus intra-layer transitions. It is set by
#' default to 0.5. See more details below.
#'
#' @return A square sparse adjacency matrix created with the \code{Matrix}
#' package.
#'
#' @seealso \code{\link{create.multiplex},\link{normalize.multiplex.adjacency},
#' \link{compute.transition.matrix}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#' compute.adjacency.matrix(multiObject)
#'
#'@import igraph
#'@import Matrix
#'@import parallel
#'@import doParallel
#'@import foreach
#'@import iterators
#'@importFrom methods as
#'@export
compute.adjacency.matrix <- function(x, delta = 0.5) {
    if (!isMultiplex(x) && !isMultiplexHet(x)) {
        stop("Not a Multiplex or Multiplex Heterogeneous object")
    }
    if (delta > 1 || delta <= 0) {
        stop("Delta should be between 0 and 1")
    }

    N <- x$Number_of_Nodes_Multiplex
    L <- x$Number_of_Layers

    Layers_Names <- names(x)[seq(L)]

    counter <- 0
    Layers_List <- lapply(x[Layers_Names],function(x){

        counter <<- counter + 1
        if (is_weighted(x)) {
            Adjacency_Layer <-  as_adjacency_matrix(x,sparse = TRUE, 
            attr = "weight")
        } else {
            Adjacency_Layer <-  as_adjacency_matrix(x,sparse = TRUE)
        }

        Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),
                                        order(colnames(Adjacency_Layer))]
        colnames(Adjacency_Layer) <-
            paste0(colnames(Adjacency_Layer), "_", counter)
        rownames(Adjacency_Layer) <-
            paste0(rownames(Adjacency_Layer), "_", counter)

        Adjacency_Layer <- normalize.multiplex.adjacency(Adjacency_Layer)
        Adjacency_Layer
    })

    ## We impose delta=0 in the monoplex case.
    if (L==1){
      SupraAdjacencyMatrix <- bdiag(unlist(Layers_List))
      SupraAdjacencyMatrix <- as(SupraAdjacencyMatrix, "dgCMatrix")
      return(SupraAdjacencyMatrix)
    }
    
    cl <- makeCluster(detectCores())

    ## IDEM_MATRIX.
    Idem_Matrix <- Matrix::Diagonal(N, x = 1)

    MyColNames <- unlist(lapply(Layers_List, function (x) unlist(colnames(x))))
    MyRowNames <- unlist(lapply(Layers_List, function (x) unlist(rownames(x))))
    names(MyColNames) <- c()
    names(MyRowNames) <- c()
    SupraAdjacencyMatrix <- (1-delta)*(bdiag(unlist(Layers_List)))
    colnames(SupraAdjacencyMatrix) <-MyColNames
    rownames(SupraAdjacencyMatrix) <-MyRowNames

    offdiag <- (delta/(L-1))*Idem_Matrix

    i <- seq_len(L)
    Position_ini_row <- 1 + (i-1)*N
    Position_end_row <- N + (i-1)*N
    j <- seq_len(L)
    Position_ini_col <- 1 + (j-1)*N
    Position_end_col <- N + (j-1)*N
    
    combinations <- expand.grid(seq_len(L), seq_len(L))
    filtered_combinations <- subset(combinations, Var1 != Var2)

    modified_matrices <- foreach(column = seq_len(L), .combine = 'cbind') %dopar% {
        # get rows that match the 
        column_mask_for_rows <- filtered_combinations$Var2 == column
        matching_rows <- filtered_combinations[column_mask_for_rows, ]

        matching_layers <- matching_rows$Var1
        col_ini =  Position_ini_col[column]
        col_end =  Position_end_col[column]

        modified_matrix <- SupraAdjacencyMatrix[, col_ini:col_end]


        col_vector <- c()
        row_vector <- c()
        for (sub_layer in matching_layers){
            row_ini = Position_ini_row[sub_layer]
            row_end = Position_end_row[sub_layer]
            row_vector <- c(row_vector, row_ini:row_end)
            col_vector <- c(col_vector, col_ini:col_end)
        }

        col_vector <- col_vector - (col_ini -1) 

        modified_vectors <- cbind(row_vector, col_vector)
        modified_matrix[modified_vectors] <- unique(offdiag[offdiag>0]) # this is scalar... ought to be 1U

        return(drop0(modified_matrix))

    }
    SupraAdjacencyMatrix = modified_matrices

    stopCluster(cl)

    SupraAdjacencyMatrix <- as(SupraAdjacencyMatrix, "dgCMatrix")

    return(SupraAdjacencyMatrix)
}

## Roxy Documentation comments
#' Computes column normalization of an adjacency matrix
#'
#' \code{normalize.multiplex.adjacency} is a function to compute the column
#' normalization of a sparse matrix of the package \code{Matrix}.
#'
#' @usage normalize.multiplex.adjacency(x)
#'
#' @param x A \code{Matrix} object describing an adjacency matrix of a network.
#'
#' @return A square sparse column normalized matrix created with the
#' \code{Matrix} package.
#'
#' @seealso \code{\link{compute.adjacency.matrix},
#' \link{Random.Walk.Restart.Multiplex}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#' AdjMatrix <- compute.adjacency.matrix(multiObject)
#' normalize.multiplex.adjacency(AdjMatrix)
#'
#'@import Matrix
#'@export

normalize.multiplex.adjacency <- function(x)
{
  if (!is(x,"dgCMatrix")){
    stop("Not a dgCMatrix object of Matrix package")
  }
  
  column_sums <- Matrix::colSums(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)

  col_names <- colnames(x)
  row_names <- rownames(x)
  
  # Find columns with non-zero sums
  mask <- column_sums != 0
  
  # Create a diagonal matrix with 1 / column_sums for non-zero columns, 1 for zero columns
  inv_col_sums <- numeric(length(column_sums))
  inv_col_sums[mask] <- 1 / column_sums[mask]
  inv_col_sums[!mask] <- 1  # or 0 if you want to zero out those columns
  
  D_inv <- Diagonal(x = inv_col_sums)
  
  # Multiply x by D_inv on the right to normalize columns
  x_norm <- x %*% D_inv
  
  colnames(x_norm) <- col_names
  rownames(x_norm) <- row_names

  x_norm
}

## Roxy Documentation comments
#' Computes row normalization of a matrix
#'
#' \code{row.normalize.matrix} is a function to compute the row
#' normalization of a sparse matrix of the package \code{Matrix}.
#'
#' @usage row.normalize.matrix(x)
#'
#' @param x A \code{Matrix} object describing a matrix.
#'
#' @return A sparse column normalized matrix created with the
#' \code{Matrix} package.
#'
#' @seealso \code{\link{compute.adjacency.matrix},
#' \link{Random.Walk.Restart.Multiplex}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#' AdjMatrix <- compute.adjacency.matrix(multiObject)
#' row.normalize.matrix(AdjMatrix)
#'
#'@import Matrix
#'@export

row.normalize.matrix <- function(x) {
  if (!is(x, "dgCMatrix")){
    stop("Not a dgCMatrix object of Matrix package")
  }

  row_sums <- Matrix::rowSums(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)
  row_sums_masked <- row_sums[ row_sums != 0 ]
  mask_names <- names(row_sums_masked)
  x[mask_names, ] <- x[mask_names, ] / row_sums_masked
  return(x)
}

#' Computes the transition matrix of a multiplex network
#'
#' \code{compute.transition.matrix.homogeneous} is a function to compute the
#' transition matrix of a multiplex network provided as a \code{Multiplex} object.
#'
#' @usage compute.transition.matrix.homogeneous(x,delta = 0.5,transpose = TRUE)
#'
#' @details The parameter \code{delta} sets the probability to change between
#' layers at the next step. If delta = 0, the particle will always remain
#' in the same layer after a non-restart iteration. On the other hand, if
#' delta = 1, the particle will always change between layers, therefore
#' not following the specific edges of each layer.
#'
#' @param x A \code{Multiplex} object describing a multiplex network generated
#' by the function \code{create.multiplex}.
#' @param delta A numeric value between 0 and 1. It sets the probability
#' of performing inter-layer versus intra-layer transitions. It is set by
#' default to 0.5. See more details below.
#' @param tranpose A boolean flag. It is used to return a transposed a matrix
#' that is transposed such that the columns are normalized.
#'
#' @return A square sparse adjacency matrix created with the \code{Matrix}
#' package.
#'
#' @seealso \code{\link{create.multiplex},\link{row.normalize.matrix},
#' \link{compute.transition.matrix.heterogeneous}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#' compute.transition.matrix.homogeneous(multiObject)
#'
#'@import igraph
#'@import Matrix
#'@import parallel
#'@import doParallel
#'@import foreach
#'@import iterators
#'@importFrom methods as
#'@export

compute.transition.matrix.homogeneous <- function(x,
                                                  delta = 0.5,
                                                  transpose = TRUE) {
  if (!isMultiplex(x)) {
    stop("Not a Multiplex object")
  }
  if (delta > 1 || delta <= 0) {
    stop("Delta should be between 0 and 1")
  }

  N <- x$Number_of_Nodes_Multiplex
  L <- x$Number_of_Layers

  Layers_Names <- names(x)[seq(L)]

  counter <- 0
  Layers_List <- lapply(x[Layers_Names], function(x) {

    counter <<- counter + 1
    if (is_weighted(x)) {
      Adjacency_Layer <-  as_adjacency_matrix(x, sparse = TRUE, attr = "weight")
    } else {
      Adjacency_Layer <-  as_adjacency_matrix(x, sparse = TRUE)
    }

    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),
                                       order(colnames(Adjacency_Layer))]
    colnames(Adjacency_Layer) <- paste0(colnames(Adjacency_Layer), "_", counter)
    rownames(Adjacency_Layer) <- paste0(rownames(Adjacency_Layer), "_", counter)

    # Adjacency_Layer <- row.normalize.matrix(Adjacency_Layer)
    Adjacency_Layer <- normalize.multiplex.adjacency(Adjacency_Layer) * (1-delta)

    Adjacency_Layer
  })

  # If L==1, we skip the intra-layer computation
  if (L == 1) {
    TransMatrix <- bdiag(unlist(Layers_List))

    # if (transpose) {
    #   TransMatrix <- t(TransMatrix)
    # }

    TransMatrix <- as(TransMatrix, "dgCMatrix")
    return(TransMatrix)
  }

  build_block_matrix_dgCMatrix <- function(diag_matrices, off_diag_matrix) {
    L <- length(diag_matrices)
    
    # Check off_diag_matrix is dgCMatrix
    if (!inherits(off_diag_matrix, "dgCMatrix")) {
      stop("off_diag_matrix must be a dgCMatrix")
    }
    
    # Check all diag_matrices are dgCMatrix
    if (!all(sapply(diag_matrices, function(m) inherits(m, "dgCMatrix")))) {
      stop("All diagonal matrices must be dgCMatrix")
    }
    
    # Get dimensions of diagonal blocks
    row_dims <- sapply(diag_matrices, nrow)
    col_dims <- sapply(diag_matrices, ncol)
    
    off_rows <- nrow(off_diag_matrix)
    off_cols <- ncol(off_diag_matrix)
    
    if (!(off_rows == row_dims[1] && off_cols == col_dims[1])) {
      stop("off_diag_matrix dimensions must match diagonal block dimensions")
    }
    
    # Compute block offsets
    row_offsets <- c(0, cumsum(row_dims))
    col_offsets <- c(0, cumsum(col_dims))
        
    # Compute total nnz = sum of all diagonal blocks nnz + number of off-diagonal blocks * nnz(off_diag_matrix)
    nnz_diag <- sum(sapply(diag_matrices, function(m) length(m@x)))
    nnz_off <- length(off_diag_matrix@x)
    total_off_blocks <- L * L - L
    total_nnz <- nnz_diag + total_off_blocks * nnz_off
    
    # Preallocate vectors
    i_all <- integer(total_nnz)
    j_all <- integer(total_nnz)
    x_all <- numeric(total_nnz)
    
    pos <- 1  # position in vectors to write next entries
    
    # Helper function to copy dgCMatrix triplets into output vectors with offset
    copy_block <- function(mat, row_off, col_off) {
      # mat@p gives column pointers, mat@i row indices (0-based), mat@x values
      nnz_block <- length(mat@x)
      # Preallocate local vectors
      i_block <- integer(nnz_block)
      j_block <- integer(nnz_block)
      x_block <- numeric(nnz_block)
      
      k <- 1
      ncol_mat <- length(mat@p) - 1
      for (col in seq_len(ncol_mat)) {
        start <- mat@p[col] + 1         # +1 for R 1-based indexing
        end <- mat@p[col + 1]
        if (end >= start) {
          idx_range <- start:end
          n_entries <- length(idx_range)
          
          # Row indices: mat@i is 0-based, add 1 for R indexing, then add row offset
          i_block[k:(k + n_entries - 1)] <- mat@i[idx_range] + 1 + row_off
          # Column indices: current column index + col offset
          j_block[k:(k + n_entries - 1)] <- col + col_off
          # Values
          x_block[k:(k + n_entries - 1)] <- mat@x[idx_range]
          
          k <- k + n_entries
        }
      }
      list(i = i_block, j = j_block, x = x_block)
    }
    
    # Loop over blocks to fill preallocated vectors
    for (i in seq_len(L)) {
      print(paste0("Combining blocks in layer", i))
      for (j in seq_len(L)) {


        if (i == j) {
          # Diagonal block
          block_trip <- copy_block(diag_matrices[[i]], row_offsets[i], col_offsets[j])
        } else {
          # Off diagonal block
          block_trip <- copy_block(off_diag_matrix, row_offsets[i], col_offsets[j])
        }
        
        nnz_block <- length(block_trip$x)
        idx_range <- pos:(pos + nnz_block - 1)
        i_all[idx_range] <- block_trip$i
        j_all[idx_range] <- block_trip$j
        x_all[idx_range] <- block_trip$x
        pos <- pos + nnz_block
      }
    }
    
    total_rows <- sum(row_dims)
    total_cols <- sum(col_dims)
    
    sparseMatrix(i = i_all, j = j_all, x = x_all, dims = c(total_rows, total_cols))
  }

  MyColNames <- unname(unlist(lapply(Layers_List, function(x) unlist(colnames(x)))))
  MyRowNames <- unname(unlist(lapply(Layers_List, function(x) unlist(rownames(x)))))

  ## IDEM_MATRIX.
  offdiag <- Matrix::Diagonal(N, x = (delta / (L - 1)) )
  offdiag <- as(offdiag, "dgCMatrix")

  TransMatrix <- build_block_matrix_dgCMatrix(Layers_List, offdiag)
  print(paste0("Constructed non-normalized matrix"))

  # Column normalize to account for nodes with zero edges in a layer
  TransMatrix <- normalize.multiplex.adjacency(TransMatrix)

  colnames(TransMatrix) <- MyColNames
  rownames(TransMatrix) <- MyRowNames

  # if (transpose) {
  #   TransMatrix <- t(TransMatrix)
  # }

  TransMatrix <- as(TransMatrix, "dgCMatrix")

  return(TransMatrix)
}

#' Computes the transition matrix of a heterogeneous multiplex network
#'
#' \code{compute.transition.matrix.heterogeneous} is a function to compute the
#' transition matrix of a multiplex network provided as a \code{Multiplex} object.
#'
#' @usage compute.transition.matrix.heterogeneous(x,delta = 0.5,transpose = TRUE)
#'
#' @details The parameter \code{delta} sets the probability to change between
#' layers at the next step. If delta = 0, the particle will always remain
#' in the same layer after a non-restart iteration. On the other hand, if
#' delta = 1, the particle will always change between layers, therefore
#' not following the specific edges of each layer.
#'
#' @param x A \code{Multiplex} object describing a multiplex network generated
#' by the function \code{create.multiplex}.
#' @param delta A numeric value between 0 and 1. It sets the probability
#' of performing inter-layer versus intra-layer transitions. It is set by
#' default to 0.5. See more details below.
#' @param tranpose A boolean flag. It is used to return a transposed a matrix
#' that is transposed such that the columns are normalized.
#'
#' @return A square sparse adjacency matrix created with the \code{Matrix}
#' package.
#'
#' @seealso \code{\link{create.multiplex},\link{row.normalize.matrix},
#' \link{compute.transition.matrix.homogeneous}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#' compute.transition.matrix.heterogeneous(multiObject)
#'
#'@import igraph
#'@import Matrix
#'@import parallel
#'@import doParallel
#'@import foreach
#'@import iterators
#'@importFrom methods as
#'@export

compute.transition.matrix.heterogeneous <- function(x,
                                                    delta1 = 0.5,
                                                    delta2 = 0.5,
                                                    lambda = 0.5) {
  if (!isMultiplexHet(x)) {
    stop("Not a Multiplex Heterogeneous object")
  }

  if (delta1 > 1 || delta1 <= 0) {
    stop("Delta1 should be between 0 and 1")
  }

  if (delta2 > 1 || delta2 <= 0) {
    stop("Delta2 should be between 0 and 1")
  }

  if (lambda > 1 || lambda <= 0) {
    stop("Lambda should be between 0 and 1")
  }

  # Compute transition matrix for homogeneous networks
  Transition_Multiplex_Network1 <-
    compute.transition.matrix.homogeneous(x$Multiplex1, delta1, FALSE)
  Transition_Multiplex_Network2 <-
    compute.transition.matrix.homogeneous(x$Multiplex2, delta2, FALSE)

  # Normalize Supra Bipartite Matrix
  Transition_Multiplex1_Multiplex2 <-
    row.normalize.matrix(x$BipartiteNetwork)
  Transition_Multiplex2_Multiplex1 <-
    row.normalize.matrix(t(x$BipartiteNetwork))

  # Scale each matrix
  Transition_Multiplex_Network1 <- (1 - lambda) * Transition_Multiplex_Network1
  Transition_Multiplex_Network2 <- (1 - lambda) * Transition_Multiplex_Network2
  Transition_Multiplex1_Multiplex2 <- lambda * Transition_Multiplex1_Multiplex2
  Transition_Multiplex2_Multiplex1 <- lambda * Transition_Multiplex2_Multiplex1

  Transition_Multiplex_Heterogeneous_Matrix_1 <-
    cbind(Transition_Multiplex_Network1, Transition_Multiplex1_Multiplex2)
  Transition_Multiplex_Heterogeneous_Matrix_2 <-
    cbind(Transition_Multiplex2_Multiplex1, Transition_Multiplex_Network2)
  Transition_Multiplex_Heterogeneous_Matrix <-
    rbind(Transition_Multiplex_Heterogeneous_Matrix_1,
          Transition_Multiplex_Heterogeneous_Matrix_2)

  Transition_Multiplex_Heterogeneous_Matrix <-
    row.normalize.matrix(Transition_Multiplex_Heterogeneous_Matrix)
  Transition_Multiplex_Heterogeneous_Matrix <-
    t(Transition_Multiplex_Heterogeneous_Matrix)
  
  Transition_Multiplex_Heterogeneous_Matrix <- as(Transition_Multiplex_Heterogeneous_Matrix, "dgCMatrix")

  return(Transition_Multiplex_Heterogeneous_Matrix)
}

## Roxy Documentation comments
#' Performs Random Walk with Restart on a Multiplex Network
#'
#' \code{Random.Walk.Restart.Multiplex} is a function to perform a Random Walk
#' with Restart on a Multiplex network (on a \code{Multiplex} object). See
#' more details about the algorithm below.
#'
#' @usage Random.Walk.Restart.Multiplex(...)
#'
#' @details Random Walk with Restart simulates an imaginary particle that
#' starts on a seed(s) node(s) and follows randomly the edges of a network. At
#' each step, there is a restart probability, r, meaning that the particle comes
#' back to the seed(s). The extension to multiplex networks allows the particle
#' to explore different monoplex networks (layers). At each step, the particle
#' can also jump to the same node in a different layer.
#'
#'
#' \itemize{
#' \item \code{Seeds}: A vector containing the name of the different seed
#' node(s). It's mandatory to provide at least one seed. The seed(s) node(s)
#' should belong to any of the layers. The length of this vector should be
#' smaller than the total number of nodes in the multiplex network.
#' \item \code{r}: A numeric value representing the restart probability on the
#' seeds for the random walker. It must be between 0 and 1. It is set by default
#' to 0.7, which is the most common value in this kind of approaches. It means
#' that, at each step, the walker has a 70\% of probability of coming back to
#' one of the seeds.
#' \item \code{tau}: A numeric vector containing the probability of restarting
#' in the nodes of the different layers of the multiplex. In the example below,
#' we define the node 1 as the seed node. However, we can find this node in both
#' layers. Therefore, the walker can restart in any of these seed nodes. It is
#' a way to give different relevance (weight) to the different layers.
#' }
#'
#'
#' @param x An object of the \code{Matrix} package describing a column
#' normalized adjacency matrix of a multiplex network.
#' @param MultiplexObject A \code{Multiplex} object generated by the function
#' \code{create.multiplex} representing a multiplex network.
#' @param Seeds A vector containing the names of the seeds for the Random
#' Walk algorithm. See more details below.
#' @param r A numeric value between 0 and 1. It sets the probability of
#' restarting to a seed node after each step. See more details below.
#' @param tau A vector containing the probability of restart on the seeds
#' of the different layers (layers weights). It must have the same length than
#' the number of layers of the multiplex network. The sum of its components
#' divided by the number of layers must be 1. See more details below.
#' @param weights An optional numeric vector containing  weights for all seeds.
#' It must have the same length as the \code{Seeds} vector.
#' These should all be positive values > 0. If not set then all seeds are 
#' treated equally by the RWR. If weights are set then the walker is biased to 
#' nodes with higher weight. See more details below.
#' @param MeanType The user can choose one of the following options: 
#' c("Geometric","Arithmetic","Sum"). These options represent the different way
#' to combine the RWR score for the same node in different layers. By default 
#' and recommended Geometric (Geometric Mean.). Arithmetic is the arithmetic 
#' mean and sum just sum all the scores for the same node across the different
#' layers. 
#' @param DispResults The user can choose one of the following options: 
#' c("TopScores","Alphabetic"). These options represent the way the RWR results
#' would be presented. By default, and recommended, the nodes would be ordered
#' by score. This option is also required to properly run the 
#' \code{create.multiplexNetwork.topResults} and 
#' \code{create.multiplexHetNetwork.topResults} functions
#' @param ... Further arguments passed to \code{Random.Walk.Restart.Multiplex}
#'
#' @return A \code{RWRM_Results} object. It contains a sorted ranking of all
#' the nodes of the multiplex network, except the seeds, along with their score.
#' In addition, it contains in a different field the nodes used as seeds.
#'
#' @seealso \code{\link{create.multiplex},
#' \link{compute.adjacency.matrix}, \link{normalize.multiplex.adjacency},
#' \link{isRWRM_Results}, \link{Random.Walk.Restart.MultiplexHet}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#' AdjMatrix <- compute.adjacency.matrix(multiObject)
#' AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
#' SeedNodes <- c(1)
#' Random.Walk.Restart.Multiplex(AdjMatrixNorm,multiObject,SeedNodes)
#'
#'@import igraph
#'@import Matrix
#'@rdname Random.Walk.Restart.Multiplex
#'@export
Random.Walk.Restart.Multiplex <- function(...) {
    UseMethod("Random.Walk.Restart.Multiplex")
}

#'@rdname Random.Walk.Restart.Multiplex
#'@export
Random.Walk.Restart.Multiplex.default <- function(x, MultiplexObject, Seeds, 
    r=0.7,tau,weights=1,MeanType="Geometric", DispResults="TopScores", ...){
        
    ### We control the different values.
    if (!is(x,"dgCMatrix")){
        stop("Not a dgCMatrix object of Matrix package")
    }
        
    if (!isMultiplex(MultiplexObject)) {
        stop("Not a Multiplex object")
    }
        
    L <- MultiplexObject$Number_of_Layers
    N <- MultiplexObject$Number_of_Nodes
        
    Seeds <- as.character(Seeds)
    if (length(Seeds) < 1 | length(Seeds) >= N){
        stop("The length of the vector containing the seed nodes is not 
           correct")
    } else {
        if (!all(Seeds %in% MultiplexObject$Pool_of_Nodes)){
            stop("Some of the seeds are not nodes of the network")
        }
    }
        
    if (r >= 1 || r <= 0) {
        stop("Restart partameter should be between 0 and 1")
    }
        
    if(missing(tau)){
        tau <- rep(1,L)/L
    } else {
        tau <- as.numeric(tau)
        if (sum(tau)/L != 1) {
            stop("The sum of the components of tau divided by the number of 
             layers should be 1")
        }
    }
        
    if(!(MeanType %in% c("Geometric","Arithmetic","Sum"))){
        stop("The type mean should be Geometric, Arithmetic or Sum")
    }
        
    if(!(DispResults %in% c("TopScores","Alphabetic"))){
        stop("The way to display RWRM results should be TopScores or
       Alphabetic")
    }
        
    ## We define the threshold and the number maximum of iterations for
    ## the random walker.
    Threeshold <- 1e-10
    NetworkSize <- ncol(x)
        
    ## We initialize the variables to control the flux in the RW algo.
    residue <- 1
    iter <- 1
        
    ## We compute the scores for the different seeds
    ## Seeds are given initial weightings to bias the walker.
    Seeds_Score <- get.seed.scoresMultiplex(Seeds,L,tau, weights)
    
    ## We define the prox_vector(The vector we will move after the first RWR
    ## iteration. We start from The seed. We have to take in account
    ## that the walker with restart in some of the Seed nodes, depending on
    ## the score we gave in that file).
    prox_vector <- matrix(0,nrow = NetworkSize,ncol=1)
    prox_vector[which(colnames(x) %in% Seeds_Score[,1])] <- (Seeds_Score[,2])
    prox_vector  <- prox_vector/sum(prox_vector)
    
    restart_vector <-  prox_vector
    
    debug=F
    dfdebug <- data.frame()  ### remove
    while(residue >= Threeshold){
        old_prox_vector <- prox_vector
        prox_vector <- (1-r)*(x %*% prox_vector) + r*restart_vector
        
        residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
        if(debug) {
            print(iter)
            dfdebug <- rbind(dfdebug, data.frame(node=rownames(prox_vector), 
                                                       score=prox_vector[,1],
                                                       isSeed = colnames(x) %in% Seeds_Score[,1],
                                                       iter=iter))
        }
        iter <- iter + 1;
    }
    if(debug){
        print("converged!") 
        write.table(dfdebug,"~/devwork/_SCRATCH/RWRdebug.txt", sep="\t", col.names = T, row.names = F, quote=F)
    }
    NodeNames <- character(length = N)
    Score = numeric(length = N)
        
    rank_global <- data.frame(NodeNames = NodeNames, Score = Score)
    rank_global$NodeNames <- gsub("_1$", "", row.names(prox_vector)[seq_len(N)])
        
    if (MeanType=="Geometric"){
        rank_global$Score <- geometric.mean(as.vector(prox_vector[,1]),L,N)    
    } else {
        if (MeanType=="Arithmetic") {
            rank_global$Score <- regular.mean(as.vector(prox_vector[,1]),L,N)    
        } else {
            rank_global$Score <- sumValues(as.vector(prox_vector[,1]),L,N)    
        }
    }
        
    if (DispResults=="TopScores"){
        ## We sort the nodes according to their score.
        Global_results <- 
            rank_global[with(rank_global, order(-Score, NodeNames)), ]
            
        ### We remove the seed nodes from the Ranking and we write the results.
        Global_results <- 
            Global_results[which(!Global_results$NodeNames %in% Seeds),]
    } else {
        Global_results <- rank_global    
    }
        
    rownames(Global_results) <- c()
    
    RWRM_ranking <- list(RWRM_Results = Global_results,Seed_Nodes = Seeds)
        
    class(RWRM_ranking) <- "RWRM_Results"
    return(RWRM_ranking)
}



#' @method print RWRM_Results
#' @export
print.RWRM_Results <- function(x,...)
{
    cat("Top 10 ranked Nodes:\n")
    print(head(x$RWRM_Results,10))
    cat("\nSeed Nodes used:\n")
    print(x$Seed_Nodes)
}

## Roxy Documentation comments
#' Computes the transition matrix of a multiplex and heterogeneous network
#'
#' \code{compute.transition.matrix} is a function to compute the transition
#' matrix of a multiplex heterogeneous network provided as a \code{MultiplexHet}
#' object.
#'
#' @usage compute.transition.matrix(x,lambda = 0.5, delta1=0.5, delta2=0.5)
#'
#' @details We clarify the role of the different parameters in this point:
#' \itemize{
#' \item \code{lambda}: For a given node, if a bipartite association exists,
#' the particle can either jump between networks or stay in the current
#' graph with a probability given by this parameter. The closer lambda is to
#' one, the higher is the probability of jumping between networks following
#' bipartite interactions.
#' \item \code{delta}: This parameter sets the probability to change between
#' layers at the next step. If delta = 0, the particle will always remain
#' in the same layer after a non-restart iteration. On the other hand, if
#' delta = 1, the particle will always change between layers, therefore
#' not following the specific edges of each layer.
#' }
#'
#' @param x A \code{MultiplexHet} object describing a multiplex and
#' heterogeneous network generated by the function
#' \code{create.multiplexHet}.
#'
#' @param lambda A numeric value between 0 and 1. It sets the probability of
#' jumping within a network or change to the other network of the heterogeneous
#' system. It is set by default to 0.5. See more details below.
#' @param delta1 A numeric value between 0 and 1. It sets the probability
#' of performing inter-layer versus intra-layer transitions in the first 
#' multiplex. It is set by default to 0.5. See more details below.
#' @param delta2 A numeric value between 0 and 1. It sets the probability
#' of performing inter-layer versus intra-layer transitions in the second 
#' multiplex. It is set by default to 0.5. See more details below.
#' 
#' @return A square sparse transition matrix created with the \code{Matrix}
#' package. It is the transition matrix for the Random Walk with Restart on
#' Multiplex and Heterogeneous networks algorithm.
#'
#' @seealso \code{\link{create.multiplexHet},
#' \link{compute.adjacency.matrix}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject_1 <- create.multiplex(list(m1=m1,m2=m2))
#' h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
#' bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))
#' multiObject_2 <- create.multiplex(list(h1=h1))
#' multiHetObject <- create.multiplexHet(multiObject_1, multiObject_2,
#'     bipartite_relations)
#' compute.transition.matrix(multiHetObject)
#'
#'@import igraph
#'@import Matrix
#'@export
compute.transition.matrix <- function(x,lambda = 0.5, delta1=0.5,delta2=0.5)
{
    if (!isMultiplexHet(x)) {
        stop("Not a Multiplex Heterogeneous object")
    }
    
    if (delta1 > 1 || delta1 <= 0) {
        stop("Delta1 should be between 0 and 1")
    }
    
    if (delta2 > 1 || delta2 <= 0) {
        stop("Delta2 should be between 0 and 1")
    }
    
    if (lambda > 1 || lambda <= 0) {
        stop("Lambda should be between 0 and 1")
    }
    
    Number_Nodes_1 <- x$Multiplex1$Number_of_Nodes_Multiplex
    Number_Layers_1 <- x$Multiplex1$Number_of_Layers
    Number_Nodes_2 <- x$Multiplex2$Number_of_Nodes_Multiplex
    Number_Layers_2 <- x$Multiplex2$Number_of_Layers
    SupraBipartiteMatrix <- x$BipartiteNetwork
    
    message("Computing adjacency matrix of the first Multiplex network...")
    AdjMatrix_Multiplex1 <- compute.adjacency.matrix(x$Multiplex1,delta1)
    AdjMatrix_Multiplex1 <- normalize.multiplex.adjacency(AdjMatrix_Multiplex1)
    
    message("Computing adjacency matrix of the second Multiplex network...")
    ## We have to sort the adjacency matrix
    AdjMatrix_Multiplex2 <- compute.adjacency.matrix(x$Multiplex2,delta2)
    AdjMatrix_Multiplex2 <- normalize.multiplex.adjacency(AdjMatrix_Multiplex2)
    
    ## Transition Matrix for the inter-subnetworks links
    message("Computing inter-subnetworks transitions...")
    Transition_Multiplex1_Multiplex2 <- 
        get.transition.multiplex1.multiplex2(Number_Nodes_1,Number_Layers_1,
            Number_Nodes_2,Number_Layers_2,SupraBipartiteMatrix,lambda)
    
    Transition_Multiplex2_Multiplex1 <- 
        get.transition.multiplex2.multiplex1(Number_Nodes_1,Number_Layers_1,
            Number_Nodes_2, Number_Layers_2, SupraBipartiteMatrix,lambda)
    
    ## Transition Matrix for the intra-subnetworks links
    message("Computing intra-subnetworks transitions...")
    Transition_Multiplex_Network1 <- 
        get.transition.multiplex(Number_Nodes_1, Number_Layers_1, lambda, 
            AdjMatrix_Multiplex1,SupraBipartiteMatrix)
    Transition_Multiplex_Network2 <- 
        get.transition.multiplex(Number_Nodes_2,Number_Layers_2, lambda,
            t(AdjMatrix_Multiplex2),t(SupraBipartiteMatrix))
    
    ## We generate the global transiction matrix and we return it.
    message("Combining inter e intra layer probabilities into the global 
          Transition Matix")
    Transition_Multiplex_Heterogeneous_Matrix_1 <-
        cbind(Transition_Multiplex_Network1, Transition_Multiplex1_Multiplex2)
    Transition_Multiplex_Heterogeneous_Matrix_2 <-
        cbind(Transition_Multiplex2_Multiplex1, Transition_Multiplex_Network2)
    Transition_Multiplex_Heterogeneous_Matrix <-
        rbind(Transition_Multiplex_Heterogeneous_Matrix_1,
              Transition_Multiplex_Heterogeneous_Matrix_2)
    
    return(t(Transition_Multiplex_Heterogeneous_Matrix))
}

## Roxy Documentation comments
#' Performs Random Walk with Restart on a Multiplex and Heterogeneous
#' Network
#'
#' \code{Random.Walk.Restart.MultiplexHet} is a function to perform a Random
#' Walk with Restart on a Multiplex and Heterogeneous network (on a
#' \code{MultiplexHet} object). See more details about the algorithm below.
#'
#' @usage Random.Walk.Restart.MultiplexHet(...)
#'
#' @details Random Walk with Restart simulates an imaginary particle which
#' starts on a seed(s) node(s) and follows randomly the edges of a network. At
#' each step, there is a restart probability, r, meaning that the particle comes
#' back to the seed(s). The extension to multiplex networks allows the particle
#' to explore different monoplex networks (layers). At each step, the particle
#' can also jump to the same node in a different layer. The extension to
#' heterogeneous networks allows the particle to jump between nodes of different
#' nature thanks to bipartite relationships between them. We can combine both,
#' the multiplex and heterogeneous extension, by allowing the particle to jump
#' from a node in every layer of the multiplex network to the other network, and
#' the other way around.
#'
#' \itemize{
#' \item \code{Multiplex_Seed_Nodes}: A vector containing the name of the
#' different seed node(s) of the multiplex network. It's mandatory to provide at
#' least one seed (taking in account both types of seeds) The seed(s) node(s)
#' should belong to any of the layers of the multiplex network. The length of
#' this vector should be smaller than the total number of nodes in the multiplex
#' network.
#' \item \code{SecondNet_Seed_Nodes}: A vector containing the name of the
#' different seed node(s) of the second network. It's mandatory to provide at
#' least one seed (taking in account both types of seeds) The seed(s) node(s)
#' should belong to the second network. The length of this vector should be
#' smaller than the total number of nodes in the second network.
#' \item \code{r}: A numeric value representing the restart probability on the
#' seeds for the random walker. It must be between 0 and 1. It is set by default
#' to 0.7, which is the most common value in this kind of approaches. It means
#' that, at each step, the walker has a 70\% of probability of coming back to
#' one of the seeds.
#' \item \code{tau}: A numeric vector containing the probability of restarting
#' in the nodes of the different layers of the multiplex. In the example below,
#' we define the node 1 as the seed node. However, we can find this node in both
#' layers. Therefore, the walker can restart in any of these seed nodes. It is
#' a way to give different relevance (weight) to the different layers.
#' \item \code{eta}: A numeric value between 0 and 1 controlling the
#' probability of restarting in the nodes of each network. In the example below,
#' we define the node 1 as a multiplex seed node and "E" as a second network
#' seed node. Therefore, the walker can restart either in the seed 1 or in the
#' seed "E" with different probabilities (it is a way to give more relevance
#' to the different components of the heterogeneous system). If eta < 0.5
#' the particle will be more likely to restart in one of the multiplex seeds.
#' }
#'
#' @param x An object of the \code{Matrix} package describing the possible
#' transitions in a multiplex and heterogeneous network.
#' @param MultiplexHet_Object A \code{MultiplexHet} object generated by the
#' function \code{create.multiplexHet} representing a multiplex
#' and heterogeneous network.
#' @param Multiplex1_Seeds A vector containing the names of the seeds of
#' the first multiplex network for the Random Walk algorithm. See more details 
#' below.
#' @param Multiplex2_Seeds A vector containing the names of the seeds of
#' the second multiplex network for the Random Walk algorithm. See more details 
#' below.' 
#' @param r A numeric value between 0 and 1. It sets the probability of
#' restarting to a seed node after each step. See more details below.
#' @param tau1 A vector containing the probability of restart on the seeds
#' of the different multiplex layers (layers weights) for the first multiplex.
#' It must have the same length than the number of layers of the multiplex. 
#' network. The sum of its components divided by the number of layers must be 1.
#' See more details below.
#' @param tau2 A vector containing the probability of restart on the seeds
#' of the different multiplex layers (layers weights) for the second multiplex.
#' It must have the same length than the number of layers of the multiplex.
#' network. The sum of its components divided by the number of layers must be 1.
#' See more details below.
#' @param eta A numeric value between 0 and 1. It controls the probability of
#' restarting in each network of the heterogeneous system (Multiplex or second
#' network). See more details below.
#' @param MeanType The user can choose one of the following options: 
#' c("Geometric","Arithmetic","Sum"). These options represent the different way
#' to combine the RWR score for the same node in different layers. By default 
#' and recommended Geometric (Geometric Mean.). Arithmetic is the arithmetic 
#' mean and sum just sum all the scores for the same node across the different
#' layers. 
#' @param DispResults The user can choose one of the following options: 
#' c("TopScores","Alphabetic"). These options represent the way the RWR results
#' would be presented. By default, and recommended, the nodes would be ordered
#' by score. This option is also required to properly run the 

#' @param ... Further arguments passed to
#' \code{Random.Walk.Restart.MultiplexHet}
#'
#' @return A \code{RWRMH_Results} object. It contains three sorted rankings: 
#' i) The first one contains the global results, i.e. the nodes of both 
#' multiplex networks along with their score; 
#' ii) The second one contains the nodes of the first multiplex network,
#' except the seeds, along with their score.
#' iii) The last one contains the nodes of the second multiplex network, 
#' excepting the seeds, along with their score
#' In addition, it contains one more field describing the nodes used as seeds.
#'
#' @seealso \code{\link{create.multiplexHet},
#' \link{compute.transition.matrix}, \link{Random.Walk.Restart.Multiplex}
#' \link{isRWRMH_Results}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject_1 <- create.multiplex(list(m1=m1,m2=m2))
#' h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
#' bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))
#' multiObject_2 <- create.multiplex(list(h1=h1))
#' multiHetObject <- create.multiplexHet(multiObject_1, multiObject_2,
#'     bipartite_relations)
#' MultiHetTranMatrix <- compute.transition.matrix(multiHetObject)
#' Multiplex1_Seeds <- c(1)
#' Multiplex2_Seeds <- c("E")
#' Random.Walk.Restart.MultiplexHet(MultiHetTranMatrix,
#'     multiHetObject,Multiplex1_Seeds,Multiplex2_Seeds)
#'
#'@import igraph
#'@import Matrix
#'@rdname Random.Walk.Restart.MultiplexHet
#'@export
Random.Walk.Restart.MultiplexHet <- function(...){
    UseMethod("Random.Walk.Restart.MultiplexHet")
}

#'@rdname Random.Walk.Restart.MultiplexHet
#'@export
Random.Walk.Restart.MultiplexHet.default <- function(x, MultiplexHet_Object, 
    Multiplex1_Seeds, Multiplex2_Seeds, r=0.7,tau1,tau2,eta=0.5,
    MeanType="Geometric", DispResults="TopScores",...){
        
    ## We control the different values.
    if (!"dgCMatrix" %in% class(x)){
      stop("Not a dgCMatrix object of Matrix package")
    }
    
    if (!isMultiplexHet(MultiplexHet_Object)) {
      stop("Not a Multiplex Heterogeneous object")
    }
        
    NumberLayers1 <- MultiplexHet_Object$Multiplex1$Number_of_Layers
    NumberNodes1 <- MultiplexHet_Object$Multiplex1$Number_of_Nodes_Multiplex
    NumberLayers2 <- MultiplexHet_Object$Multiplex2$Number_of_Layers
    NumberNodes2 <- MultiplexHet_Object$Multiplex2$Number_of_Nodes_Multiplex
        
    All_nodes_Multiplex1 <- MultiplexHet_Object$Multiplex1$Pool_of_Nodes
    All_nodes_Multiplex2 <- MultiplexHet_Object$Multiplex2$Pool_of_Nodes
        
    MultiplexSeeds1 <- as.character(Multiplex1_Seeds)
    MultiplexSeeds2 <- as.character(Multiplex2_Seeds)
        
    if (length(MultiplexSeeds1) < 1 & length(MultiplexSeeds2) < 1){
      stop("You did not provided any seeds")
    } 
    if (length(MultiplexSeeds1) >= NumberNodes1 ||
        length(MultiplexSeeds2) >= NumberNodes2) {
      stop("The length of some of the vectors containing the seed nodes 
            is not correct")
    }
    if (!all(MultiplexSeeds1 %in% All_nodes_Multiplex1)){
      stop("Some of the input seeds are not nodes of the first input network")
    } 
    if (!all(All_nodes_Multiplex2 %in% All_nodes_Multiplex2)){
      stop("Some of the inputs seeds are not nodes of the second input network")
    } 
    if (r >= 1 || r <= 0) {
      stop("Restart partameter should be between 0 and 1")
    }
    if (eta >= 1 || eta <= 0) {
      stop("Eta partameter should be between 0 and 1")
    }   
    if (missing(tau1)){
      tau1 <- rep(1,NumberLayers1)/NumberLayers1
    } else {
      tau1 <- as.numeric(tau1)
      if (sum(tau1)/NumberLayers1 != 1) {
        stop("The sum of the components of tau1 divided by the number of layers should be 1")
      }
    }

    if (missing(tau2)){
      tau2 <- rep(1,NumberLayers2)/NumberLayers2
    } else {
      tau2 <- as.numeric(tau2)
      if (sum(tau2)/NumberLayers2 != 1) {
        stop("The sum of the components of tau2 divided by the number of layers should be 1")
      }
    }

    if (!(MeanType %in% c("Geometric","Arithmetic","Sum"))) {
      stop("The type mean should be Geometric, Arithmetic or Sum")
    }

    if (!(DispResults %in% c("TopScores","Alphabetic"))) {
      stop("The way to display RWRM results should be TopScores or Alphabetic")
    }

    ## We define the threshold and the number maximum of iterations
    ## for the random walker.
    Threeshold <- 1e-10
    NetworkSize <- ncol(x)

    ## We initialize the variables to control the flux in the RW algo.
    residue <- 1
    iter <- 1

    ## We compute the scores for the different seeds.
    Seeds_Score <- get.seed.scores.multHet(Multiplex1_Seeds, Multiplex2_Seeds,eta,
        NumberLayers1,NumberLayers2,tau1,tau2)

    ## We define the prox_vector(The vector we will move after the first
    ## RWR iteration. We start from The seed. We have to take in account
    ## that the walker with restart in some of the Seed genes,
    ## depending on the score we gave in that file).
    prox_vector <- matrix(0,nrow = NetworkSize,ncol=1)

    prox_vector[which(colnames(x) %in% Seeds_Score[,1])] <- (Seeds_Score[,2])

    prox_vector  <- prox_vector/sum(prox_vector)
    restart_vector <-  prox_vector

    while (residue >= Threeshold) {
      old_prox_vector <- prox_vector
      prox_vector <- (1-r)*(x %*% prox_vector) + r*restart_vector
      residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
      iter <- iter + 1;
    }

    IndexSep <- NumberNodes1*NumberLayers1
    prox_vector_1 <- prox_vector[1:IndexSep,]
    prox_vector_2 <- prox_vector[(IndexSep+1):nrow(prox_vector),]

    NodeNames1 <- gsub("_1$", "",names(prox_vector_1)[seq_len(NumberNodes1)])
    NodeNames2 <- gsub("_1$", "",names(prox_vector_2)[seq_len(NumberNodes2)])

    if (MeanType=="Geometric") {
      rank_global1 <- geometric.mean(prox_vector_1,NumberLayers1,NumberNodes1)    
      rank_global2 <- geometric.mean(prox_vector_2,NumberLayers2,NumberNodes2)   
    } else if (MeanType=="Arithmetic") {
      rank_global1 <- regular.mean(prox_vector_1,NumberLayers1,NumberNodes1)    
      rank_global2 <- regular.mean(prox_vector_2,NumberLayers2,NumberNodes2)     
    } else {
      rank_global1 <- sumValues(prox_vector_1,NumberLayers1,NumberNodes1)    
      rank_global2 <- sumValues(prox_vector_2,NumberLayers2,NumberNodes2) 
    }

    Global_results <- data.frame(NodeNames = c(NodeNames1,NodeNames2), 
        Score = c(rank_global1,rank_global2))

    Multiplex1_results <-data.frame(NodeNames = NodeNames1,Score = rank_global1)
    Multiplex2_results <-data.frame(NodeNames = NodeNames2,Score = rank_global2)
    Seeds <- c(Multiplex1_Seeds, Multiplex2_Seeds)

    if (DispResults=="TopScores") {
      ## We sort the nodes according to their score.
      Global_results <- 
          Global_results[with(Global_results, order(-Score, NodeNames)), ]
      Multiplex1_results <-
          Multiplex1_results[with(Multiplex1_results, order(-Score, NodeNames)),] 
      Multiplex2_results <-
          Multiplex2_results[with(Multiplex2_results, order(-Score, NodeNames)),] 

      ### We remove the seed from the individual networks 
      Multiplex1_results <- 
          Multiplex1_results[which(!Multiplex1_results$NodeNames %in% Seeds),]
      Multiplex2_results <- 
          Multiplex2_results[which(!Multiplex2_results$NodeNames %in% Seeds),]      
    } else {
      Global_results <- Global_results
      Multiplex1_results <- Multiplex1_results
      Multiplex2_results <- Multiplex2_results
    }

    rownames(Global_results) <- c()

    RWRMH_ranking <- list(RWRMH_GlobalResults = Global_results, 
        RWRMH_Multiplex1 = Multiplex1_results, 
        RWRMH_Multiplex2 = Multiplex2_results,
        Seed_Nodes = Seeds)

    class(RWRMH_ranking) <- "RWRMH_Results"
    return(RWRMH_ranking)
}

#' @method print RWRMH_Results
#' @export
print.RWRMH_Results <- function(x,...)
{
    cat("Top 10 ranked global nodes:\n")
    print(head(x$RWRMH_GlobalResults,10))
    cat("\nTop 10 ranked nodes from the first Multiplex:\n")
    print(head(x$RWRMH_Multiplex1,10))
    cat("\nTop 10 ranked nodes from the second Multiplex:\n")
    print(head(x$RWRMH_Multiplex2,10))
    cat("\nSeeds used:\n")
    print(x$Seed_Nodes)
}
