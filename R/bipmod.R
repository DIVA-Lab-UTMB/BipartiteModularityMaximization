#' Partition bipartite network into non-overlapping biclusters, by optimizing bipartite modularity.
#'
#'
#' This function partitions a bipartite network into non-overlapping biclusters by optimizing bipartite modularity defined in Barber (2007) using the bipartite version of the algorithm described in Treviño (2015).
#'
#'
#' The function takes as input a bipartite network represented as an incidence matrix (using a matrix or a data frame) with non-negative values (the row sums and column sums must be positive, to ensure there are no disconnected nodes). The function partitions the rows and columns into non-overlapping submatrices (biclusters), and outputs the membership of rows and columns to a partition, and modularity (Q) representing the quality of the partitioning.
#'
#' @param incid_mat Incidence matrix of a bipartite network.
#' @param ITER A positive integer representing the number of iterations used to maximizing modularity, (default=10).
#'
#' @return MODULARITY Modularity value (Q).
#' @return ASSIGN Integer labels representing partition of rows followed by columns in same order as incidence matrix.
#' @export
#'
#' @examples
#' data(example_data)
#' bipmod(example_data)
#' @references Barber, M. J. (2007). Modularity and community detection in bipartite networks. Physical Review E, 76(6), 066102. <doi:10.1103/PhysRevE.76.066102>
#' @references Trevino, S., Nyberg, A., Del Genio, C. I., & Bassler, K. E. (2015). Fast and accurate determination of modularity and its effect size. Journal of Statistical Mechanics: Theory and Experiment, 2015(2), P02003. <doi:10.1088/1742-5468/2015/02/P02003>
bipmod=function(incid_mat,ITER=10){
  return(CoClust(nrow(incid_mat),ncol(incid_mat),as.double(t(incid_mat)),ITER=ITER))
}

.onUnload <- function (libpath) {
  library.dynam.unload("BipartiteModularityMaximization", libpath)
}
