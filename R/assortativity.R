##- Heatmaps fns
Var1=Var2=value=NULL

##---- mat_to_heatmap ----
#' Heatmap from matrix in base R
#' @param m matrix
#' @param xlab label for columns
#' @param ylab label for rows
#' @param values display values in cells
#' @param ... For title, color, etc. plot options
#' @seealso \code{\link{level_gg}}
#' @examples
#' m <- matrix(1:9, nrow=3, byrow = TRUE)
#' colnames(m) <- paste0("C", 1:3)
#' rownames(m) <- paste0("R", 1:3)
#' m
#' mat_to_heatmap(m)
#' @export
mat_to_heatmap <- function(m, xlab = "", ylab = "", values = TRUE, ...){
  graphics::image(1:ncol(m), 1:nrow(m), t(m), axes = FALSE, xlab = xlab, ylab = ylab, ...)
  graphics::axis(1, 1:ncol(m), colnames(m))
  graphics::axis(2, 1:nrow(m), rownames(m))
  if(values){
    for (x in 1:ncol(m))
      for (y in 1:nrow(m))
        graphics::text(x, y, m[y,x])
  }
}

##---- edgelist_to_mixingmatrix ----
#' Weighted edgelist to mixing matrix
#'
#' Add weights for each combination of \code{from} and \code{to} pairs in a normalised matrix. Used for mixing matrix of probability that there is one transmission between donor_j (column) to recipient_i (row) ( total sums to 1 ). Or donor = row and recip = column if \code{transpose = FALSE}.
#' @param df dataframe in form [from, to, weight]
#' @param transpose If TRUE, \code{from} in column and \code{to} in row
#' @examples
#' ( df <- data.frame(from = c('a','a','b','b'), to = c('a','b','b','b') , ip = c(.2,.8,.5,.5)) )
#' ( m <- edgelist_to_mixingmatrix(df) ) # mixing matrix
#' m / rowSums(m)[row(m)] # If transpose = TRUE,
#' # P(donor_j|recipient_i) or
#' # (probability that recipient row has donor column).
#' # Every row sums to 1
#' @export
# edgelist_to_mixingmatrix <- function(df, transpose = TRUE){
#   ud <- as.vector(unique(df[,1][!is.na(df[,1])]))
#   ur <- as.vector(unique(df[,2][!is.na(df[,2])]))
#   ud <- ud[order(ud)]
#   ur <- ur[order(ur)]
#   m <- matrix(0, nrow = length(ud), ncol = length(ur))
#   rownames(m) <- ud
#   colnames(m) <- ur
#   for (i in 1:nrow(df)){
#     ii <- which(ud == df[i, 1])
#     jj <- which(ur == df[i, 2])
#     w <- df[i, 3]
#     m[ii, jj] <- m[ii, jj] + w
#   }
#   if(transpose) m <- t(m)
#   return( m / sum(m) )
# }
##- igraph version
# edgelist_to_mixingmatrix <- function(df, transpose = TRUE){
#   g <- igraph::graph.data.frame(df[,1:2])
#   igraph::E(g)$weight <- df[,3]
#   m2 <- igraph::as_adjacency_matrix(g, attr="weight")
#   m3 <- as.matrix(m2[sort(rownames(m2)), sort(colnames(m2))])
#   m <- m3 / sum(m3)
#   if(transpose) m <- t(m)
#   m
# }
##- tapply() version, quicker
edgelist_to_mixingmatrix <- function(df, transpose = TRUE){
  m <- tapply(df[, 3], list(df[,1], df[,2]), sum)
  m[is.na(m)] <- 0L
  if(transpose) m <- t(m)
  m / sum(m)
}

##---- mat2assortmat ----
#' Assortativity matrix
#'
#' Excess assortment between categories relative to random allocation
#' @param mat mixing matrix
#' @export
mat2assortmat <- function(mat){
rs <- rowSums(mat)
cs <- colSums(mat)
s <- sum(mat)
M <- matrix(0, nrow=nrow(mat), ncol = ncol(mat))
A <- matrix(0, nrow=nrow(mat), ncol = ncol(mat))
colnames(A) <- colnames(mat)
rownames(A) <- rownames(mat)
k <- nrow(mat)
for (i in 1:k) for (j in 1:k){
  M[i,j] <- rs[i] * cs[j] / s
  A[i,j] <- ( mat[i,j] - M[i,j] )/ M[i,j]
}
A
}

##---- newman_r ----
#' Newman assortativity coefficient from original paper
#' @references Newman, M.E.J., 2003. Mixing patterns in networks. Phys. Rev. E 67, 026126.
#' \url{https://doi.org/10.1103/PhysRevE.67.026126}
#' @param mat mixing matrix
#' @examples
#' ##- Table 1
#' e <- matrix(data = c(0.258,0.016,0.035,0.013,
#'   0.012,0.157,0.058,0.019,
#'   0.013,0.023,0.306,0.035,
#'   0.005,0.007,0.024,0.016),
#'   nrow = 4, ncol = 4, byrow = TRUE)
#' rownames(e) <- c("MenBlack", "MenHispanic", "MenWhite", "MenOther")
#' colnames(e) <- c("WomenBlack", "WomenHispanic", "WomenWhite", "WomenOther")
#' e
#' newman_r(e)
#'
#' ##- Bipartite graph
#' g <- igraph::graph_from_incidence_matrix(e, weighted = TRUE)
#' igraph::V(g)$color <- igraph::V(g)$type + 1
#' igraph::plot.igraph(g, edge.width = igraph::E(g)$weight*20)
#'
#' ##- example from Newmann p.12
#' M <- matrix(data = c(50,50,0,50,50,0,0,0,2)/202, nrow = 3, byrow = TRUE)
#' M
#' ##- expected value under random allocation:
#' # $E = \sum_i{m_{ij}} \otimes \sum_ j{m_{ij}} / \sum_{ij}{m_{ij}}$
#' E <- rowSums(M) %o% colSums(M)/sum(M)
#' E
#' ##- assortativity matrix
#' A <- (M - E)/E
#' A
#' identical(mat2assortmat(M), A) # TRUE
#' ##- assortativity coef: $r = (\Tr(M) - \Tr(E)) / (1 - \Tr(E))$
#' r <- (sum(diag(M)) - sum(diag(E)))/ (1 - sum(diag(E)))
#' r
#' identical(r, newman_r(M)) # 0.029
#' @export
newman_r <- function(mat){
  mat <- mat / sum(mat)
  sum_aibi <- sum(rowSums(mat)*colSums(mat))
  r <- (sum(diag(mat)) - sum_aibi) / (1 - sum_aibi)
  return(r)
}

##---- plots_grid ----
#' Don't know what this does
#'
#' yet
#' @inheritParams level_gg
#' @param ... extra parms for lattice::levelplot
#' @importFrom grDevices heat.colors
#' @export
plots_grid <- function(mat, contrast = 1, transpose = TRUE, ... ){
  if(transpose) mat <- t(mat)
  p <-  vector("list", length(mat))
  for (i in 1:length(mat)){
    M <- mat[[i]]
    if(contrast!=1) {
      M <- (abs(M)^(1/contrast)) * sign(M)
    }
    p[[i]] <- lattice::levelplot( M, col.regions = heat.colors, main = names(mat)[i], ...)
  }
  return(p)
}

##---- level_gg ----
#' Heatmap from matrix with ggplot
#'
#' Forced to square rendering
#' @seealso \code{\link{mat_to_heatmap}}
#' @param mat matrix
#' @param contrast change heatmap contrast (independent of plotted values)
#' @param title main title
#' @param legend logical for legend (if contrast transformation, legend shows untransformed values)
#' @param transpose transpose matrix
#' @param labx label x
#' @param laby label y
#' @param labz label of values, passed to legend
#' @param rotate_y rotate y-axis labels
#' @param abbrev use abbreviate or identity
#' @param values write values (untransformed by contrast) in cells
#' @param fmt format values, i.e. round, signif, etc. [default: identity]
#' @param textsize size of values
#' @param col1 color low values [default: red]
#' @param col2 color high values [default: yellow]
#' @examples
#' m <- matrix(1:9, nrow=3, byrow = TRUE)
#' colnames(m) <- paste0("Column", 1:3)
#' rownames(m) <- paste0("Row long", 1:3)
#' m
#' level_gg(m)
#' level_gg(m, legend = FALSE, values = TRUE)
#' @export
level_gg <- function(mat, contrast = 1, title = NULL,
                     legend = FALSE,
                     transpose = TRUE,
                     labx = "donor", laby = "recipient", labz = "ip",
                     rotate_y = TRUE,
                     abbrev = identity,
                     values = FALSE,
                     fmt = identity,
                     textsize = NULL,
                     col1 = "red", col2 = "yellow"){
  trf <- function(M) (abs(M)^(1/contrast)) * sign(M) #log10 #identity
  untrf <- function(y) signif(sign(y)*exp(contrast*log(abs(y))), 2) # untransformed for legend (reciprocate)
  if(transpose) mat <- t(mat)
  df0 <- reshape2::melt(mat)
  longData <- reshape2::melt(trf(mat))
  g <- ggplot(longData, aes(x = factor(Var1), y = factor(Var2)))
  g0 <- g + geom_raster(aes(fill = value)) +
    scale_fill_gradient(name = labz, labels = untrf,
                        low = col1, high = col2) +
    labs(x = labx, y= laby, title = title) +
    scale_x_discrete(expand = c(0, 0), labels = abbrev) +
    scale_y_discrete(expand = c(0, 0), labels = abbrev) +
    theme_bw()  +
    theme(aspect.ratio=1) # force square output
  if(rotate_y){g0 <- g0 +  theme(axis.text.y = element_text(angle = 90, hjust = 0.5))}
  if(values){
    if(is.null(textsize)){
      g0 <- g0 + geom_text(aes(label = fmt(df0$value))) # auto size
    } else {
      g0 <- g0 + geom_text(aes(label = fmt(df0$value)), size = textsize)
    }
  }
  if(legend){
    return(g0)
  } else {
    g1 <- g0 + theme(legend.position="none")
    return(g1)
  }
}

##---- rorder ----
#' reorder matrix factor level
#' @param m matrix
#' @param o character vector of rownames/colnames
#' @param normalize logical to normalize if matrix subsetted
#' @importFrom utils type.convert
#' @export
rorder <- function(m, o, normalize = FALSE){
  if(is.integer(type.convert(rownames(m)))){
    rownames(m) <- colnames(m) <- o
    om <- m
  } else {
    om <- m[o[o %in% rownames(m)], o[o %in% rownames(m)]]
  }
  if(normalize) om <- om / sum(om)
  return(om)
}

##---- interpret_mat_gg ----
#' Series of 4 heatmaps and assortativity coefficient.
#'
#' \itemize{
#' \item mixing: is the probability that there is one transmission between donor of type j to a recipient of type i (cumulative infector probabilities, sums to 1).
#' \item std_recip: is the probability for a recipient of type i to have a donor of type j (Each row sums to 1).
#' \item std_donor: is the probability for a donor of type j to transmit to a recipient of type i (Each col sums to 1).
#' \item assortment: is the excess assortment between types relative to random allocation. Summarized by Newman's assortativity coefficient: r.
#' }
#' @param l_all_mixmat List of mixing matrices computed by \code{\link{edgelist_to_mixingmatrix}}, by variable and by group
#' @param var First list level of \code{l_all_mixmat}
#' @param var_order Optional named list (by variable) of row/col names to reorder matrices
#' @param group Second list level of \code{l_all_mixmat}
#' @param ... option to \code{\link{level_gg}}, i.e. contrast = 2, transpose = FALSE, tl = "something", abbrev_y = abbreviate
#' @export
interpret_mat_gg <- function(l_all_mixmat, var, var_order = NULL, group, ...){
  l_mixmat <- l_all_mixmat[[var]]
  if(!is.null(var_order)) { # re-order matrix
    l_mixmat <- lapply(l_mixmat, function(mat)
      rorder(m = mat, o = var_order[[var]], normalize = TRUE))
  }
  l_transmat_r <- lapply(l_mixmat, function(m) m / rowSums(m)[row(m)])
  l_transmat_c <- lapply(l_mixmat, function(m) m / colSums(m)[col(m)])
  l_assortmat <- lapply(l_mixmat, mat2assortmat)
  ll <- list(mixing = l_mixmat[[group]], std_recip = l_transmat_r[[group]], std_donor = l_transmat_c[[group]], assortment = l_assortmat[[group]])
  pp <- lapply(1:length(ll), function(i) level_gg(mat = ll[[i]], title = names(ll)[i], legend = FALSE,...))
  print(paste(var, ":", "Newman r = ", round(newman_r(ll[[1]]), 3) ))
  do.call(gridExtra::grid.arrange, c(pp, ncol = ifelse(length(pp)==1,1,2)))
}

##---- mat_to_igraph ----
#' Bidirectional weighted graph from matrix
#'
#' With \code{igraph}. See net0_sims_slv.R
#' @param m adjacency matrix
#' @param layout default: layout.fruchterman.reingold
#' @param labels weights as labels
#' @param ... passed to plot.igraph.
#' @examples
#' k <- 3
#' m <- matrix(runif(k*k,0,1), nrow = k, byrow = TRUE)
#' colnames(m) <- rownames(m) <- LETTERS[1:k]
#' m
#' mat_to_igraph(m)
#' @export
mat_to_igraph <- function(m, layout = NULL, labels = FALSE, ...){
  #source("R/custom_circle_igraph.R")

  {
    # http://lists.gnu.org/archive/html/igraph-help/2013-03/msg00030.html
    # library(igraph)

    mycircle <- function(coords, v=NULL, params) {
      vertex.color <- params("vertex", "color")
      if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
      }
      vertex.size  <- 1/200 * params("vertex", "size")
      if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
      }
      vertex.frame.color <- params("vertex", "frame.color")
      if (length(vertex.frame.color) != 1 && !is.null(v)) {
        vertex.frame.color <- vertex.frame.color[v]
      }
      vertex.frame.width <- params("vertex", "frame.width")
      if (length(vertex.frame.width) != 1 && !is.null(v)) {
        vertex.frame.width <- vertex.frame.width[v]
      }

      mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
             vertex.size, vertex.frame.width,
             FUN=function(x, y, bg, fg, size, lwd) {
               graphics::symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                                 circles=size, add=TRUE, inches=FALSE)
             })
    }

    igraph::add.vertex.shape("fcircle", clip=igraph::igraph.shape.noclip,
                             plot=mycircle, parameters=list(vertex.frame.color=1,
                                                            vertex.frame.width=1))

    # plot(graph.ring(10), vertex.shape="fcircle", vertex.frame.color=rainbow(10),
    #     vertex.frame.width=1:10/2)
  }

  g <- igraph::graph_from_adjacency_matrix(m, mode = "directed", weighted = TRUE)
  if(is.null(layout)) layout <- igraph::layout_with_fr(g) # layout.fruchterman.reingold
  # l <- layout_with_kk(g) # Kamada Kawai
  # l <- layout_with_dh(g)
  # l <- layout_nicely(g)

  ## edge width from 0.5 to 5
  ( w <- igraph::E(g)$weight / min( igraph::E(g)$weight ) )
  ( k = log(10) / log(max(w)) )
  igraph::E(g)$width <- abs(w)^k/2 * sign(w)
  igraph::E(g)$arrow.width <- max(igraph::E(g)$width)
  igraph::E(g)$arrow.size <- 0.5
  igraph::E(g)$curved <- 0.3
  igraph::E(g)$color = "grey"
  ## vertices
  igraph::V(g)$color="white"
  igraph::V(g)$label.color="black"
  igraph::V(g)$label= igraph::V(g)$name # c("Aaaa", "B", "Cccccccc")
  igraph::V(g)$label.family="sans"
  igraph::V(g)$label.font=2 #bold
  igraph::V(g)$shape="fcircle"
  igraph::V(g)$frame.width=2
  igraph::V(g)$frame.color="black"

  if(labels) igraph::E(g)$label <- round(igraph::E(g)$weight * 100) # no way to control label placement
  igraph::plot.igraph(g, layout = layout, ...)
  ## http://igraph.org/r/doc/plot.common.html
}
