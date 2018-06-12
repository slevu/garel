#require(Rcpp)
#require(RcppArmadillo)
#require(ape)
#sourceCpp( 'dist2root.cpp' )

#' Root to tip distance
#' @param tre phylo tree
#' @useDynLib garel
#' @importFrom Rcpp sourceCpp
#' @export
dist2root <- function( tre )
{
  with( tre,{
    n <- length(tip.label)
    edge.length <- c( edge.length, 0 ) #root
    edge <- rbind( edge, c( -1, n+1) ) # root
    node2edge.index <- sapply( 1:(n + Nnode ), function(u)  which(edge[,2]==u) )
    node2parent <- sapply( node2edge.index, function(i) edge[i,1]  )
    node2edge.length <- edge.length[node2edge.index]
    setNames( d2root_cpp(n, node2parent, node2edge.length)[,1], tip.label )
  })
}


# example usage
if (F )
{
  n <- 10
  tre <- rtree( n )
  d0 <- dist2root( tre )
  # compare to ape func's
  d1 <- dist.nodes(tre)[n+1,1:n]
  plot( d0, d1 )
  abline( b= 1, a = 0)

  # larger tree
  n <- 1e3
  tre <- rtree( n )
  st.0 <- system.time( { d0 <- dist2root( tre ) } )
  # compare to ape func's
  st.1 <- system.time( { d1 <- dist.nodes(tre)[n+1,1:n] } )
  X11(); plot( d0, d1 )
  abline( b= 1, a = 0)
  print( st.0 )
  print( st.1)

  # and a BIG tree
  n <- 2e4
  tre <- rtree( n )
  st.0 <- system.time( { d0 <- dist2root( tre ) } ) # 22 sec
  # compare to ape func's
  #~ DOES NOT RUN, Error: cannot allocate vector of size 11.9 Gb
  #~ st.1 <- system.time( { d1 <- dist.nodes(tre)[n+1,1:n] } )
  print( st.0 )

}
