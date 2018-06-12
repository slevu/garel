# dates = NULL; outgroup = NULL; verbose = TRUE; ci = '-f 10'
##---- lsd_fn ----
#' @title Generic LSD usage
#' @description Version 0.3beta
#' @param tr Phylo tree
#' @param dates Dates of tips, by default extracted from accession label [default=NULL]
#' @param outgroup Labels of outgroup tips
#' @param ci Confidence interval calculation (\code{ci = '-f 10'} or \code{ci = NULL})
#' @param root.constraint In the form \code{root.constraint = c("mrca(a,b,c,d)", "b(2000,2001)")}
#' @param verbose Write details
#' @export
lsd_fn <- function(tr, dates = NULL, outgroup = NULL, root.constraint = NULL, verbose = TRUE, ci = '-f 10'){ # ci = NULL
  # root.constraint = c("mrca(a,b,c,d)", "b(2000,2001)")
  INPUTDAT <-  tempfile(pattern = "LSD_dates_", fileext = ".txt")
  INPUTNWK <- tempfile(pattern = "LSD_tree_", fileext = ".nwk")
  if(is.null(dates)){
    dates <- sapply(tr$tip.label, date.from.sacc.year)
  }
  if(!is.null(outgroup)){
    INPUTOG <-  tempfile(pattern = "LSD_outgroup_", fileext = ".txt")
    m <- matrix(c(length(outgroup), outgroup), byrow = FALSE)
    if (verbose) print(paste("Write labels of outgroup in", INPUTOG))
    write.table(m, file = INPUTOG, quote = FALSE, row.names = FALSE, col.names = FALSE)
    arg_out <- paste("-g", INPUTOG)
    dates <- dates[!names(dates) %in% outgroup]
  } else arg_out <- "-r a"

  write.table( rbind(cbind( names(dates), dates), root.constraint), file = INPUTDAT, col.names=c(as.character(length(dates) + (length(root.constraint)!=0)), '') , quote = F, row.names=F)

  if (verbose) print(paste(length(dates), "dates as LSD input"))
  ape::write.tree(tr, file = INPUTNWK)
  #file.show(INPUTDAT); file.show(INPUTNWK)
  if( any(grep("Darwin", Sys.info())) ){ # Mac
    lsd.cmd <- '~/Documents/softwares/lsd-0.3beta/src/lsd'
  } else { # HPC
    lsd.cmd <- '/home/slevu/lsd-0.3beta/src/lsd'
  }
  arg0 <- paste("-i", INPUTNWK, "-d", INPUTDAT)
  lsd.args <- paste(arg0, arg_out, "-c -s 1200 -v 2", ci)
  cmd <- path.expand(lsd.cmd)
  parms <- path.expand(lsd.args)
  if (verbose) print(paste(cmd, parms))
  system2(cmd, args = parms, stderr = verbose, stdout = verbose)
  ##- parse results
  tmp <- basename(INPUTNWK)
  results <- list.files(dirname(INPUTNWK), full.names = TRUE)
  to <- results[grepl(tmp, results) & grepl('result|ingroup', results)]
  ##- output rates and dated tree
  metrics <- get.lsd.metrics(txt = to[ !(grepl('newick|nexus|ingroup', to)) ], ci = ci)
  if (verbose) print(metrics)
  ##- import dated tree
  dated_tree <- ape::read.tree( to[ grepl('date', to) & grepl('newick', to) ] )
  ##- output
  out <- list(metrics = metrics, dated_tree = dated_tree )
  return(out)
}

##---- find_tmrca_outliers ----
#' @title Check which seqs generate outliers in tmrca
#' @description Test UK sequences one-by-one with group of lanl seq closest to outgroup
#' @param pathmsa msa file in dnabin (*.R) or fasta (*.fas/fasta) format
#' @export
find_tmrca_outliers <- function(pathmsa = 'data/ML_TREE0/MSM_CRF02AG_nodrm_ref_og.R'){
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
  ##- find subtype of MSA, subtype of outgroup, and blast DB
  ST <- get_subtype(pathmsa)
  OGST <- outgroup_choice(ST)
  DB <- find_db(pathmsa)
  ##- import uk+lanl+og MSA
  data("ref_hiv_msa", package = "tenbrit", envir = environment())
  og <- rownames(ref_hiv_msa[[OGST]])
  lanl_lbl <- rownames(seq)[grepl('[A-Z]', rownames(seq)) & !(rownames(seq) %in% og)]
  uk_lbl <- rownames(seq)[!grepl('[A-Z]', rownames(seq))]
  ##- select lanl closest to og = hits from BLAST outgroup against subtype DB included in MSA
  a <- blast_dnabin(seq = ref_hiv_msa[[OGST]], path.database = DB, target_seqs = 5)
  include <- a$unique_hits[a$unique_hits %in% lanl_lbl]
  print(paste(length(include), 'global sequences used in tree'))
  ##- put UK date in sacc type new labels, for date extraction of lsd_fn
  dates_uk <- readRDS(get_uk_dates())
  da <- round(dates_uk[names(dates_uk) %in% uk_lbl],0)
  new_uk_lbl <- sapply(1:length(uk_lbl), function(i) paste('_', uk_lbl[i], da[i], sep = '.'))
  rownames(seq)[rownames(seq) %in% uk_lbl] <- new_uk_lbl
  ##- iterate neighborjoining tree
  ##- make tmrca estimates
  fast_tmrca <- function(lbl, plots = TRUE){
    # lbl = new_uk_lbl[1]; plots = TRUE
    suk <- lbl # NULL
    sla <- include
    msa <- seq[c(sla, suk, og),]
    d <- ape::dist.dna(msa, 'TN93')
    if(any(is.na(as.matrix(d)))){
      ## misaligned outlier at distance stage
      return(c(suk, rep(-1, 8)))
    }
    tr <- ape::nj(d)
    tr <- ape::root(tr, og, resolve.root = TRUE)
    if(plots) plot_tree_col_og(tr, tiplab = TRUE)
    # rc <- c(paste0('mrca(', paste(include, collapse = ','), ')'), "b(1900,2000)")
    o <- lsd_fn(tr = tr,  outgroup = og, root.constraint = NULL, verbose = FALSE)
    metrics <- o[[1]]
    dated_tree <- o[[2]]
    tip_checked <- setNames(ifelse(length(suk) == 0, 0L, suk), 'tip_checked')
    at_root <- setNames(dated_tree$tip.label[ape::Ntip(dated_tree)], 'at_root')
    if(plots) plot_tree_col_og(dated_tree, tmrca =  metrics['tMRCA'], tiplab = TRUE)
    return(c(tip_checked, at_root, metrics))
  }
  print('start trees and lsd')
  st <- system.time(
    mrca_results <- lapply(1:length(new_uk_lbl), function(i){
      if(i %% 100==0) {
        print(paste0("iteration: ", i))
      }
      fast_tmrca(lbl = new_uk_lbl[i], plots = FALSE)
    })
  )
  print(st)
  .ch <- unfactorDataFrame(as.data.frame(do.call(rbind, mrca_results)))
  tmrcas <- setNames(.ch$tMRCA, unname(sapply(.ch[,1], function(x) unlist(strsplit(x,'\\.'))[2])))
  # head(tmrcas[order(tmrcas)]); summary(tmrcas)
  ##- outliers
  if(FALSE) plot(tmrcas[tmrcas > 0])
  gru <- grubbs.flag(tmrcas)
  if(FALSE) plot(gru$X [gru$X > 0], col = gru[gru$X > 0,]$Outlier + 1)
  out.tmrca <- rownames(gru[gru$Outlier,])
  print(dput(out.tmrca))
  return(out.tmrca)
}

# MSAS <- list.files(path = 'data/ML_TREE0', pattern = '_nodrm_ref_og\\.(R|fas)$', full.names = TRUE)
# outliers.tmrca <- find_tmrca_outliers(pathmsa = MSAS[1])

