##---- norm_test ----
#' Test for funr
#' @param n blabla
#' @importFrom stats rnorm
#' @export
norm_test <- function(n){
  x = rnorm(n)
  print(x)
}

##---- check_st ----
#' Check subtype from LANL sequences labels
#' @param x DNAbin object or table of BLAST results
#' @param k number of subtype to return
#' @return proportion of \code{k} most frequent
# #' @export
check_st <-  function(x, k = 6){
  if(class(x) == "DNAbin") {
    subtype <- sapply(rownames(x), function(z) strsplit(z, '\\.')[[1]][1])
  } else {
    subtype <- sapply(x$sacc, function(z) strsplit(z, '\\.')[[1]][1])
  }
  tb <- table(subtype)
  round(prop.table(tb[order(tb, decreasing = TRUE)]),2)[1:k]
}

##---- list.depth ----
list.depth <- function(this, thisdepth = 0) {
  # http://stackoverflow.com/a/13433689/1270695
  if(!is.list(this)) {
    return(thisdepth)
  } else {
    return(max(unlist(lapply(this, list.depth, thisdepth = thisdepth+1))))
  }
}
# l1 <- list(a = 1, b =2)
# l2 <- list(l1, l1)
# list.depth(l1);list.depth(l2)
# l1[['b']]

##---- date.from.sacc.year ----
##- extract dates from accession number
date.from.sacc.year <- function(x){
  ddmmyy <- paste('01', '01', unlist(strsplit(x, '.', fixed = T))[3], sep = '-')
  fmt <- ifelse(nchar(ddmmyy) == 8, "%d-%m-%y", "%d-%m-%Y" )
  da <- as.Date( ddmmyy , format = fmt)
  ## midpoint year (since 01-01-0000)
  as.numeric(format(da, '%Y')) + 0.5
}

##---- get_lanl_dates ----
#' Extract and format dates from sequences labels
#'
#' For dating trees
#' @param lstinfile List of data files (in package) containing LANL sequences
#' @param outfile Path to named vector of dates (numeric in year; for LANL, midpoint year of accession)
#' @return Path to output
#' @importFrom utils data
#' @export
get_lanl_dates <- function(lstinfile = c("db_hiv1_A_noUK", "db_hiv1_B_noUK", "db_hiv1_CRF02AG_noUK", "db_hiv1_C_noUK"),
                           outfile = "data/lanl_dates.rds"){
  if(file.exists(outfile)){
    print(paste("file already there:", outfile))
  } else {
    ##- get data from package
    d <- data(list = lstinfile, package = 'tenbrit', envir = environment())
    lstdates <- lapply(d, function(x){
      lbl <- rownames(get(x))
      dates <- sapply(lbl, date.from.sacc.year)
      # print(head(dates))
      return(dates)
    })
    alldates <- do.call(c, lstdates)
    ##- add outgroup sequences labels
    data("ref_hiv_msa", package = "tenbrit", envir = environment())
    nm.outgroup <- unlist(lapply(ref_hiv_msa, function(x) rownames(x)), use.names = FALSE)
    dates.outgroup <- setNames(sapply(nm.outgroup, date.from.sacc.year), nm.outgroup)
    alldates <- c(alldates, dates.outgroup)
    saveRDS(alldates, file = outfile)
    print(paste("Save dates in", outfile))
  }
  return(outfile)
}

##---- get_uk_dates ----
#' Extract and format dates from sequences labels
#'
#' For dating trees
#' @param outfile Path to named vector of dates (numeric in year since 01-01-000)
#' @return Path to output
#' @export
##- dates for uk sequences in LSD format
##- get dates in calendar year
##---|-----------------|-----|
##   0                1970  2012

get_uk_dates <- function(outfile = "data/uk_dates.rds"){
  if(file.exists(outfile)){
    print(paste("file already there:", outfile))
  } else {
    #load(infile) # infile = "data/resi.rda"
    data("resi", package = 'tenbrit', envir = environment())
    #- origin
    origin <- as.Date('0000-01-01')
    #- import and time difference in year unit
    dates_uk <- setNames( as.numeric(resi$dateres - origin)/365.25 , resi$testindex )
    #head(dates_uk)
    saveRDS(dates_uk, file = outfile)
    print(paste("Save dates in", outfile))
  }
  return(outfile)
}

##---- name.list.group ----
#' name list of files
#' @param lst List of paths containing group info
#' @param regroup Regroup elements in a list of vectors for each unique group
#' @return Named list
#' @export
name.list.group <- function(lst, regroup = FALSE){
  .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", unlist(lst) )
  names(lst) <- regmatches(unlist(lst), .m)
  if (regroup){
    rlst <- lapply(setNames(unique(names(lst)), unique(names(lst)) ), function(l) unlist(unname(lst[names(lst) == l])) )
    return(rlst)
  } else {
    return(lst)
  }
}

##---- get_subtype ----
get_subtype <- function(infile){
  .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", infile)
  ST <- unlist(strsplit(regmatches(infile, .m), "_"))[2] # subtype
  return(ST)
}

##---- get_outgroup_for_lsd ----
## with number of labels first
## Usage
## test <- get_outgroup_for_lsd()
## file.show(test)
get_outgroup_for_lsd <- function(infile, subtype = NULL,
                                 listmsa = NULL, outfile = NULL,
                                 verbose = TRUE){
  if(is.null(subtype)){
    ST <- get_subtype(infile)
    if (verbose) print(paste("subtype of object is", ST))
    OGST <- outgroup_choice(ST)
    if (verbose) print(paste("subtype of outgroup is", OGST))
  } else {
    OGST <- subtype
  }
  if(is.null(outfile)){
    outfile <- paste0("data/", "LSD_outgroup_", OGST,".txt")
    if (file.exists(outfile)){
      print(paste(outfile, 'already exists'))
      return(outfile)
    }
  }
  if(is.null(listmsa)) {
    data("ref_hiv_msa", package = "tenbrit", envir = environment())
    listmsa <- ref_hiv_msa
  }
  lbl <- rownames(listmsa[[OGST]])
  m <- matrix(c(length(lbl),lbl), byrow = FALSE)
  if (verbose) print(paste("Write labels of outgroup in", outfile))
  write.table(m, file = outfile, quote = FALSE, row.names = FALSE, col.names = FALSE)
  return(outfile)
  }

##---- unfactor ----
unfactorDataFrame <- function( x ) {
  x <- data.frame( lapply(x, as.character), stringsAsFactors = FALSE)
  x <- data.frame( lapply(x, utils::type.convert, as.is = TRUE),
                   stringsAsFactors = FALSE)
  return(x)
}

##---- load_msa ----
#' Load MSA in various formats
#'
#' load \code{seq} object from file in .R, .fas/fasta, .rds. Remove HXB2 if any
#' @param pathmsa Path to file
#' @param verbose print content of seq
load_msa <- function(pathmsa, verbose = FALSE){
  if( grepl(".*\\.R$", basename(pathmsa)) ){
    load(pathmsa) # as seq
    if(class(seq)!='DNAbin') stop('need DNAbin')
  } else if ( grepl(".*\\.(fas|fasta)$", basename(pathmsa)) ){
    seq <- ape::read.dna(pathmsa, format = 'fasta')
  } else if (grepl(".*\\.rds$", basename(pathmsa)) ){
    l <- readRDS(pathmsa)
    w <- which(sapply(l, class)=='DNAbin')
    seq <- l[[w]]
    if(class(seq)!='DNAbin') stop('need DNAbin')
    seq <- seq[rownames(seq)!='HXB2',]
  } else {
    stop('Need filename as .R, .fas/fasta or .rds')
  }
  if (verbose) print(seq)
  return(seq)
}

##---- count_bstrees ----
#' Count files existing
#' @param path Path to trees
#' @param pattern For list.files regex
#' @param n_expected vector of numbers of bootstrap expected
#' @return list of table(group, n, file), number of files by group, number missing
#' @importFrom stats setNames
#' @export
count_bstrees <- function(path = "data/BSTREES", pattern = "ExaML_result.*finaltree\\.\\d{3}$", n_expected = 1:100){
  v <- list.files(path = path, pattern = pattern, full.names = TRUE)
  .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", v)
  gr <- regmatches(v, .m)
  nbs <- sub(".*(\\d{3}).*", "\\1", v)
  here <- data.frame(group = gr, n = as.numeric(nbs), file = v, stringsAsFactors = FALSE)
  ##- n by group
  n_by_group <- table(here$group)
  ##- number missing
  n_missing <- lapply(setNames(unique(here$group), unique(here$group)), function(x){
    setdiff(n_expected, here[here$group == x,]$n)
  })
  return(list(tab = here, n_by_group = n_by_group, n_missing = n_missing))
}

##---- get_notnaive ----
#' Get vector of \code{testindex} of treated patients
get_notnaive <- function(){
  data(list = "df", package = 'tenbrit', envir = environment())
  not_naive <- as.character(df$testindex[df$status != 'Naive'])
  return(not_naive)
}

##---- revert list of lists  -----
##- see: http://goo.gl/sFMWzD
revert_list <- function(ls) {
  # get sub-elements in same order
  x <- lapply(ls, `[`, names(ls[[1]]))
  # stack and reslice
  apply(do.call(rbind, x), 2, as.list)
}


