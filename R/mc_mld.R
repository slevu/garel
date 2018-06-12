##---- single_mc_mld ----
#' Monte-Carlo simulation of likely donor to each recipient
#'
#' Counts how many times each donor is the most likely donor to a recipient.
#' For each recipient i, one donor or "outside source" is sampled with probabilities = infector probs for each link  + (1 - indegree)
#' @param x Edge list dataframe of infector probabilities with columns (donor, recip, ip)
#' @return Vector of MLD count named by donor ID
#' @examples
#' d <- as.character(1:3)
#' r <- as.character(1:3)
#' W <- expand.grid('donor' = d, 'recip' = r, stringsAsFactors = FALSE)
#' W$ip <- runif(dim(W)[1])/length(d) # arbitrary ip
#' W # edge list
#' tapply(W$ip, W$recip, sum) # in-degree for each recipient
#' Nsim=10
#' replicate(Nsim, single_mc_mld(x = W))
#' @export
single_mc_mld <- function(x = W){ # x = edge list c(donor, recip, ip)
  nd <- length(unique(x$donor))
  ud <- setNames(rep(0, nd + 1),
                 c(unique(x$donor), 'out'))# list of unique donors
  ur <- unique(x$recip) # unique recip
  ##- list of df named by recip
  # system.time( df.r <- split(x, x$recip) ) # too slow
  ## data.table version
  DT <- data.table::as.data.table(x)
  data.table::setkey(DT, recip)
  df.r <- DT[, list(list(.SD)), by=recip]$V1
  data.table::setattr(df.r, 'names', unique(DT$recip)) # sets names by reference
  ##- sample within donors' names + "outside source", with infector probs + (1 - indegree)
  mld <- sapply(df.r, function(r){
    sample(c(r$donor, 'out'), size = 1,
           prob = c(r$ip, max(0, 1 - sum(r$ip)) ) )
  }) ## result: for each recipient either the mld or out
  n.mld <- setNames(as.vector(table(mld)), names(table(mld))) # count
  ## add count of being mld to zero vector of all donors
  v <- c(ud, n.mld)
  new.ud <- tapply(v, names(v), sum)
  stopifnot(sum(new.ud) == length(ur)) # one donor found by recip
  return( setNames(as.vector(new.ud), names(new.ud)) ) ## vector form of 1d array
}
