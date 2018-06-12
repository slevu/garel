# require(ape)
# # require(seqinr)
qacc=resi=J=W=recip=df=patientindex=rolldate=dateres=datecd4=datevl=NULL ## get rid of check notes
ref_hiv_msa=NULL

##---- check.seqs.input ----
#' @title Check if input is sequences (DNAbin in compressed R file)
#'
#' @description Check if input is \code{rds} or \code{rda}, and contains a DNAbin object
#' @param path.seqs Path to R file
#' @return DNAbin object
#' @import ape
# #' @export
check.seqs.input <- function(path.seqs = "msa.rds"){
  ##- read seqs
  if (grepl("^.*\\.(rda|rdata)$", path.seqs, ignore.case = TRUE)) {
    tmp <- load(path.seqs)
    eval(parse(text = paste("seqs <-", tmp)))
  }
  else if (grepl("^.*\\.rds$", path.seqs, ignore.case = TRUE)) {
    seqs <- readRDS(path.seqs)
  }
  else stop("need a rds or rda as input")
  if(class(seqs) != "DNAbin") stop("need dnabin object")
  return(seqs)
}

##---- dnabinToMatrix ----
#' Convert DNAbin to matrix form
#'
#' Padding with gaps to match max(length)
#' @param l Object in DNAbin format
#' @param gaps format of gaps in MSA (default = '-')
#' @return DNAbin in matrix format
#' @export
dnabinToMatrix <- function(l, gaps = '-'){
  if(class(l) != "DNAbin") stop('need DNAbin object')
  #- check same length of seqs
  if ( length(unique(sapply(l, length))) != 1 ) {
    m <- max(sapply(l, length)) #- max length
    dna <- ape::as.matrix.DNAbin(ape::as.DNAbin(lapply(l, function(x){
      #- padding with gaps
      c(ape::as.character.DNAbin(x), rep(gaps, m - length(x)) )
    } ) ) )
  } else { #- if same length
    dna <- ape::as.matrix.DNAbin(l)
  }
  return(dna)
}

##---- dashToN ----
#' Substitute characters for gaps in DNAbin object
#'
#' By default substitute dashes to 'N' either for list or matrix forms of DNAbin
#' @param s Object in DNAbin format
#' @param a old gap character (default = '-')
#' @param b new gap character (default = 'N')
#' @return DNAbin object
#' @export
dashToN <- function(s, a = '-', b = 'n'){
  if(class(s) != "DNAbin") stop('need DNAbin object')
  if(is.null(dim(s))){ # if dnabin as list
    r <- ape::as.DNAbin(
      lapply(ape::as.character.DNAbin(s), function(x) gsub(a, b, x) )
    )
  } else{ # if dnabin as matrix
    r <- ape::as.DNAbin( gsub(a, b, s, perl = TRUE) )
  }
  return(r)
}

##---- mafft_it ----
#' Align new sequences to existing alignment
#'
#' With MAFFT v7.
#' @param seq new sequences in DNAbin format (default = hxb2)
#' @param profile existing alignment (profile in MAFFT jargon) in DNAbin format
#' @param args arguments for MAFFT v7. Use \code{--thread 16} in HPC. See \url{https://mafft.cbrc.jp/alignment/server/add.html} and \url{https://mafft.cbrc.jp/alignment/software/addsequences.html}
#' @param verbose print stderr
#' @return DNAbin object
#' @export
mafft_it <- function(
  seq = ape::read.dna(system.file("ext", "hxb2.txt", package = "garel", mustWork = TRUE), format = "fasta"),
  profile,
  args = '--auto --thread -1 --keeplength --add',
  verbose = TRUE)
  {
  ##- check mafft cmd
  if( nchar(Sys.which("qsub")) ) { # if on HPC
    cmd <- "module load mafft/7.271 && mafft"
  } else {
    if( !grepl('mafft', system2('which', 'mafft', stdout = TRUE) )) stop(' Need mafft cmd installed')
    cmd <- "mafft"
  }
  ## temp files names
  nm.profile <- tempfile(pattern = "profile", fileext = '.fas')
  nm.seq <- tempfile(pattern = "seq", fileext = '.fas')
  nm.both <- tempfile(pattern = "both", fileext = '.fas')
  ## convert to fasta
  ape::write.dna(profile, nm.profile, format = 'fasta')
  ape::write.dna(seq, nm.seq, format = 'fasta')
  ## run mafft
  arg <- paste(args, nm.seq, nm.profile, '>', nm.both)
  system(paste(cmd, arg), ignore.stderr = !verbose)
  ## return dnabin
  both <- ape::read.dna(nm.both, format = 'fasta' )
  return(both)
}

##---- splitForDRM ----
#' Split MSA in subsets
#'
#' Split MSA in subsets, each aligned with HXB2
#' @param path.seqs Path to RDA or RDS file containing large alignment in DNAbin format
#' @param path.out Path of folder (based on \code{path.seqs}) to save the splitted MSAs
#' @param path.hxb2 Path to HXB2 txt file (default = location in \code{garel} package)
#' @param n.by.split Number of sequences by split (default = 100)
#' @param verbose blabla
#' @return List of splits and save splits as RDS in \code{path.out}
#' @export
splitForDRM <- function(path.seqs, path.out = ".", path.hxb2 = system.file("ext", "hxb2.txt", package = "garel", mustWork = TRUE), n.by.split = 100, verbose = TRUE) #, debut = 2253)
  {
  if(verbose) print(path.seqs)
  if (path.out == ".") {
    path.out <- paste0(sub(".rds|.rda", "", path.seqs), "_splits")
    dir.create(path.out, recursive = TRUE)
  } else {
    if (!file.exists(path.out)) stop("need an existing output folder")
  }
  seqs <- check.seqs.input(path.seqs)

  ##- force seqs to matrix
  seqsm <- dnabinToMatrix(seqs)
  max.length.seqs <- dim(seqsm)[2]

  ##- extract HXB2 fragment
  ##- if not in it, load it and align it
  if( any(grep("K03455", rownames(seqsm))) ){
    hxb2 <- seqsm[grep("K03455", rownames(seqsm)), ]
    rownames(hxb2) <- "HXB2"
    if(verbose) print('found K03455 as HXB2 reference')

  } else {

    hxb2_full <- ape::read.dna(file = path.hxb2, format = "fasta") #[, debut:(debut + max.length.seqs -1)] # fragment no longer needed if mafft keep length of profile
    rownames(hxb2_full) <- "HXB2"
    ## Align hxb2 sequence with seqsm profile of n.by.split/total sequences = one random split
    ##- get 100 longer seqs
    seqsm_nogaps <- ape::del.gaps(seqsm)
    lg <- sapply(seqsm_nogaps, length) # length when no gaps
    pro <- seqsm[order(lg, decreasing = TRUE)[1:n.by.split], ]

    ##- align hxb2 with 100 longest seqs profile
    # sp <- sample(dim(seqsm)[1], n.by.split) # random
    st <- system.time(
    rsplit <- mafft_it(seq = hxb2_full, profile = pro, verbose = FALSE) # seqsm[sp,])
    )
    if(verbose) print(paste('Align HXB2 in', round(st[3]), 's'))
    ## separate aligned hxb2
    hxb2 <- rsplit[grep("HXB2", rownames(rsplit)), ]
  }

  ##- merge aligned hxb2 sequence with
  ##- splits of full seqs
  a <- rownames(seqsm[grep("K03455", rownames(seqsm), invert = TRUE), ])
  sp <- split(a, ceiling(seq_along(a) / n.by.split  ))
  seqsm_split <- lapply(sp, function(x) rbind(hxb2, seqsm[x,]) )

  ##- save each split separately
  name.splits <- lapply(1:length(seqsm_split), function(i){
    nm <- paste0(path.out, "/", "split", i, "_", sub(".rds|.rda", "", basename(path.seqs), ignore.case = TRUE), ".rds" )
    saveRDS(seqsm_split[[i]], file = nm)
    return(nm)
  })
  if(verbose) print(paste0("writing ", length(name.splits), " splits in ", path.out, "/"))
  return(name.splits)
}

##---- rds2fasta ----
#' Convert DNAbin in fasta
#'
#' Convert in place \code{*.rds} to \code{*.fas}
#' @param infile Path to RDS file containing DNAbin object
#' @param verbose Indicates path of output file
#' @return Path to fasta file
#' @export
rds2fasta <- function(infile = "test.rds", verbose = TRUE){
  if(!grepl("^.*\\.(rds|RDS)$", infile)) stop("need a rds file as infile")
  ##- name fasta file
  outfile <- sub('.rds', '.fas', infile)
  ##- read dnabin
  a <- readRDS(infile)
  if(class(a)!="DNAbin") stop("need a DNAbin object as infile")
  ##- write fasta
  if(verbose) print(paste('writing', outfile))
  ape::write.dna(a, file = outfile, format = 'fasta')
  return(outfile)
}

##---- bigphylo_drm ----
#' Apply \code{big.phylo::seq.rm.drugresistance}
#'
#' Convert in place \code{*.rds} to \code{nodrm_*.rds}
#' @param infile Path of RDS file containing DNAbin object
#' @param verbose gives computation time
#' @param wait sleep time option to allow qsub dependency
#' @return Path of RDS file containing sequences without DRM
#' @export
bigphylo_drm <- function(infile = "test.rds", wait = 0, verbose = TRUE){
  if(!grepl("^.*\\.(rds|RDS)$", infile)) stop("need a rds file as infile")
  ##- name output file
  outname <- sub("\\..*", "", basename(infile)) # without extension
  outfile <- paste0(dirname(infile), "/", outname, "_nodrm.rds")
  ##- read dnabin
  a <- readRDS(infile)
  if(class(a)!="DNAbin") stop("need a DNAbin object as infile")
  ##-- run
  st <- system.time(
    split_nodrm <- big.phylo::seq.rm.drugresistance(a, outfile = NA)
  )
  if(verbose){
    print(st)
    print(paste0("writing ", outfile))
  }
  saveRDS(split_nodrm, file = outfile)
  Sys.sleep(wait)
  return(outfile)
}

##---- merge_splits ----
#' Merge split sequences
#'
#' Merge multiple sequences as \code{*.rds} in one MSA
#' @param nsplits Optional number of splits to wait for before merging
#' @param indir Path of RDS files containing split sequences
#' @param ptrn Regex pattern for selecting which files are splits [default = \code{"nodrm.*\\.rds$"}]
#' @param outname Name of output [default = \code{merged_seqs_nodrm}]
#' @param verbose Write stuff
#' @return Path of merged sequences RDS file
#' @export
merge_splits <- function(nsplits = 0, indir = "data/MSM_CRF02AG_splits", ptrn = ".*split.*_nodrm.rds", outname = "merged_seqs_nodrm", verbose = TRUE){
  ##- read splits with "nodrm"
  ##- each contains HXB2
  print(nsplits)
  if(nsplits == 0){
    nm.splits <- list.files(indir, pattern = ptrn, full.names = TRUE)
  } else {
    if (verbose) print(paste("Waiting for", nsplits, "nodrm files"))
    if (verbose) print(Sys.time())
    while ( length(list.files(indir, pattern = ptrn)) < nsplits ) {
      Sys.sleep(5)
    }
    if (verbose) print(Sys.time())
    nm.splits <- list.files(indir, pattern = ptrn, full.names = TRUE)
  }
  if(any(sapply(nm.splits, function(x) !grepl("^.*\\.(rds|RDS)$", x)))) stop("need rds files")
  if (verbose) print(paste(length(nm.splits), "splits"))
  ##- all splits in a list
  l <- lapply(nm.splits, function(x){
    a <- readRDS(x)[[2]] # $nodr.seq
    if(class(a)!="DNAbin") stop("need a DNAbin object as infile")
    b <- dashToN(a) # for makeblastdb
    return( b[rownames(b) != "HXB2", ]) # exclude all HXB2 but keep "K03455"
  })
  ##- rbind all splits
  c <- do.call(rbind, l)
  ##- save
  outfile <- paste0(indir, '/', outname, "_nodrm.rds")
  if (verbose) print(paste('save', outfile))
  saveRDS(c, file = outfile)
  return(outfile)
  # write.dna(c, file = sub(".rds", "_nodrm.fasta", nm.in), format = "fasta" )
}

##---- make_blast_db ----
#' Make blast DB
#'
#' Need an installation of ncbi-blast. Turns a RDS file containing MSA into a BLAST database in fasta format
#' @param pathfile Path of RDS files containing sequences
#' @param pathout Folder to store BLAST databases
#' @export
make_blast_db <- function(pathfile, pathout = "blastDB"){
  if( !nchar(Sys.which('makeblastdb')) ) stop('need ncbi-blast')
  if (pathout!="") {
    dir.create(pathout, showWarnings = FALSE)
    newpathfile <- paste(pathout, basename(pathfile), sep = '/')
    file.copy(pathfile, newpathfile , overwrite = TRUE)
    pathfile <- newpathfile
  }
  FAS <- rds2fasta(pathfile)
  arg <- paste("makeblastdb -in", FAS, "-dbtype nucl -parse_seqids")
  #  makeblastdb -in $FAS -dbtype nucl -parse_seqids
  system(arg)
  print(paste("Write DB in", FAS))
}

##---- blast_dnabin ----
#' BLAST sequences against database
#'
#' @description Runs \code{blastn} of \code{seq}(query sequences) against a database. Returns tables of
#' \describe{
#' \item{all_hits}{all matches per query sequence}
#' \item{best_hits}{best match per query sequence}
#' \item{unique_hits}{unique database sequences found as best match}
#' }
#' @param seq Input seqs in dnabin format
#' @param path.database Path to database in fasta
#' @param metric criterium for best hit. Among \code{c('pident', 'bitscore', 'evalue')} [default = \code{bitscore}]
#' @param target_seqs Maximum number of aligned sequences to keep, per query [default = 5]
#' @param subtype Print distribution of subtypes of hits
#' @param verbose explain stuff
#' @return Output a list of tables
#' @export
#' @import data.table
blast_dnabin <- function(seq, path.database, metric = 'bitscore', target_seqs = 5, subtype = FALSE, verbose = TRUE){
  if( !grepl('blastn', system2('which', 'blastn', stdout = TRUE) )) stop(' Need blastn cmd installed')
  ##- choose optimization for best_hits
  if (metric %in% c('pident', 'bitscore')) {
    FUN <-  which.max
    } else if (metric == 'evalue') {
      FUN <-  which.min
    } else { stop('Need metric in pident, evalue or bitscore') }
  ##- write fasta
  path.query <- tempfile(pattern = 'query', fileext = '.fas')
  ## del.gaps of seqs !
  ape::write.dna(ape::del.gaps(seq), file = path.query, format = 'fasta')
  ##- out file
  blastout <- tempfile(pattern = 'blastout', fileext = '.fas')
  ##- database
  ##- cmd
  cmd <- 'blastn'
  args <- paste('-db', path.database, '-query',  path.query, '-max_target_seqs', target_seqs, '-outfmt "7 qacc sacc pident evalue bitscore"') # -max_hsps 5 only for blast 2.5.0+
  if (verbose) print(paste(cmd, args))
  ##- blast
  st.blast <- system.time(
    system2(cmd, args, stdout = blastout )
  )
  if (verbose) print(st.blast)
  ##- parse result
  # columns  <-  c('qacc', 'sacc', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
  columns  <-  c('qacc', 'sacc', 'pident', 'evalue', 'bitscore')
  all_hits <- read.table(blastout, stringsAsFactors = FALSE)
  colnames(all_hits) <- columns
  ##- unique blast hits (based on max % identity)
  best_hits <- data.table::setDT(all_hits)[data.table::setDT(all_hits)[, .I[FUN(get(metric))], by = qacc]$V1]
  ## see: https://stackoverflow.com/a/28191051 and https://stackoverflow.com/a/16325932
  unique_all <- unique(all_hits$sacc)
  unique_best <- unique(best_hits$sacc)
  w <- length(unique_all)
  x <- length(unique_best)
  y <- length(unique(all_hits$qacc))
  z <- length(all_hits$sacc)
  if (verbose) print(paste0(w, ' (', round(w/z*100, 1), '%)',  ' unique hits and ', x, ' (', round(x/y*100, 1), '%)', ' unique best hits for ', y , ' query sequences'))
  ##- check subtype of hits
  if(subtype) {if (verbose) print(check_st(best_hits, 10))}
  return(list(all_hits = all_hits, unique_hits = unique_all, best_hits = unique_best))
}
##--- end ---

##---- blast_splits ----
#' BLAST splits against database
#'
#' Remove hxb2, apply blastn, output best hits
#' @param infile Path to query sequence (split) in RDS format
#' @param database Path to database in fasta
#' @param target_seqs Maximum number of aligned sequences to keep, per query [default = 5]
#' @param verbose explain stuff
#' @return Path to output list of tables
#' @export
blast_splits <- function(infile = "test.rds", database = "blastDB/hiv1_FLT_noUK_nodrm.fasta", target_seqs = 5, verbose = TRUE){
  ##- name output
  OUT <- sub('.rds', '_blast.out', infile)
  ##- read dnabin
  a <- dnabinToMatrix(readRDS(infile)[[2]]) # $nodr.seq
  a <- a[rownames(a) != 'HXB2',] # remove hxb2
  if (verbose) print(paste("query has", nrow(a), "sequences"))
  ##- list of all_hits, best_hits and unique sacc
  lst <- blast_dnabin(seq = a, path.database = database, target_seqs = target_seqs)
  hits <- lst[['all_hits']] # lst[['best_hits']]
  ##- write output (only best_hits)
  write.table(hits, file = OUT, row.names = FALSE)
  if (verbose) print(paste('writing', OUT))
  return(OUT)
}

##---- merge_blasts ----
#' Merge BLAST splits against database
#'
#' Merge split BLAST output with best hits
#' @param nblasts Optional number of splits to wait for before merging
#' @param indir Folder containing splits of BLAST output
#' @param outfile Path to output file based on original MSA name
#' @param verbose explain stuff
#' @return Table of unique hits
#' @export
merge_blasts <- function(nblasts = 0, indir = "data/MSM_CRF02AG_splits", outfile = "data/MSM_CRF02AG_splits/MSM_CRF02AG_nodrm.rds", verbose = TRUE){
  if(nblasts == 0){
    nm.blasts <- list.files(indir, pattern='blast.out$', full.names = TRUE )
  } else {
    print('waiting for all blast.out files')
    print(Sys.time())
    while ( length(list.files(indir, pattern='blast.out$')) < nblasts ) {
      Sys.sleep(5)
    }
    print(Sys.time())
    nm.blasts <- list.files(indir, pattern='blast.out$', full.names = TRUE )
  }
  ##- read blast splits
  if(verbose) print(paste("Read", length(nm.blasts), "blast results"))
  ##- merge from list
  tt <- lapply(nm.blasts, function(x) read.table(x, header = TRUE, stringsAsFactors = FALSE))
  a <- do.call(rbind, tt)
  a <- a[a[,1] != 'HXB2', ] # remove HXB2s
  ##- name output
  OUTUNI <- sub('.rds', '_unique_blast_hits.rds', outfile)
  #OUTBEST <- sub('.rds', '_best_blast_hits.rds', outfile)
  ##- parse result
  ##- unique blast hits (based on max bitscore)
  # best_hits <- a
  unique_hits <- unique(a$sacc)
  if (verbose) print(paste("Save vector of", length(unique_hits), "hits in", OUTUNI))
  ##- save unique, all and best hits
  saveRDS(unique_hits, file = OUTUNI)
  #saveRDS(best_hits, file = OUTBEST)
}

##---- find_db ----
#' Find BLAST database
#'
#' @description Find subtype on infile name and then database
#' @param infile Path to MSA
#' @param dbdir Folder of BLAST databases
#' @param verbose explain stuff
#' @return Path to database
find_db <- function(infile = "path.msa", dbdir = "blastDB", verbose = TRUE){
  ST1 <-  regmatches(basename(infile), regexpr(pattern = "_CRF02AG.|_C.|_A.|_B.", basename(infile)) )
  ST <- sub("\\.", "_", ST1)
  DB <- list.files(dbdir, paste0(".*db.*", ST, ".*\\.fas$"), full.names = TRUE)
  if (verbose) print(paste("DB =", DB))
  return(DB)
}

##---- prep_blast_align ----
#' Process MSA and its corresponding BLAST match file for alignment
#'
#' create two fasta files
#' @param infile Path to MSA
#' @param subtype Describe subtype distribution of hits
#' @param verbose explain stuff
#' @return List of path to fasta files
# #' @export
prep_blast_align <- function(infile = "data/MSM_CRF02AG_splits/MSM_CRF02AG_nodrm.rds", subtype = FALSE, verbose = TRUE){
  ##- find hits and database from infile name
  HITS <- sub("nodrm.rds", "nodrm_unique_blast_hits.rds", infile)
  if (!file.exists(HITS)) stop(paste('need corresponding BLAST results:', HITS))
  DB <- find_db(infile)
  ##- load DB
  lanl <- ape::read.dna(file = DB, format = 'fasta')
  ##- results from blast
  .seqs <- check.seqs.input(infile)
  .hits <- readRDS(HITS)
  ##- uk seqs (removing HXB2)
  uk <- ape::del.colgapsonly(dashToN(.seqs[rownames(.seqs) != 'HXB2',], 'n', '-'))
  ##- lanl reference seqs
  ref <- ape::del.colgapsonly( dashToN(lanl[rownames(lanl) %in% unique(.hits),] , 'n', '-') )
  if (verbose) print(paste("Align", dim(ref)[1], "LANL sequences and", dim(uk)[1], "sequences"))
  ##- check subtype of hits
  if(subtype) print(check_st(ref))
  ##---- align
  REFSEQS <- paste0(dirname(infile), '/', 'ref_', sub('_nodrm.rds', '', basename(infile)), '.fas')
  UKSEQS <- paste0(dirname(infile), '/', 'uk_', sub('_nodrm.rds', '', basename(infile)), '.fas')
  UKREF <- sub("nodrm.rds", "nodrm_ref.fas", infile)
  ##- save as fasta for clustal, mafft ...
  write.dna(ref, REFSEQS, format = "fasta")
  write.dna(uk, UKSEQS, format = "fasta")
  print(paste("Save:", REFSEQS, "and", UKSEQS))
  return(list(fas1 = REFSEQS, fas2 = UKSEQS, out = UKREF))
}

##---- align_mafft ----
#' Align with \code{add} from MAFFT
#'
#' Add \code{fas1} to \code{fas2}. Need \code{mafft/7.271}
#' @param fas1 Path to MSA in fasta format (small)
#' @param fas2 Path to MSA in fasta format (big)
#' @param out Path to output in fasta format
#' @param verbose explain stuff
#' @return Path to output fasta file
#' @export
align_mafft <- function(fas1 = "data/MSM_CRF02AG_splits/ref_MSM_CRF02AG.fas",
                        fas2 = "data/MSM_CRF02AG_splits/uk_MSM_CRF02AG.fas",
                        out = "data/MSM_CRF02AG_splits/MSM_CRF02AG_nodrm_ref.fas",
                        verbose = TRUE){
  ##- check mafft cmd
  if( nchar(Sys.which("qsub")) ) { # if on HPC
    cmd <- "module load mafft/7.271 && mafft"
  } else {
    if( !grepl('mafft', system2('which', 'mafft', stdout = TRUE) )) stop(' Need mafft cmd installed')
    cmd <- "mafft"
  }
  arg <- paste("--auto --thread 16 --keeplength --add", fas1, fas2, ">", out)
  system(paste(cmd, arg))
  return(out)
}

##---- prep_and_align_blast ----
#' Process MSA and align
#'
#' combine \code{prep_blast_align} and \code{align_mafft}
#' @param infile Path to MSA
#' @return Path to output
#' @export
prep_and_align_blast <- function(infile = "data/MSM_CRF02AG_splits/MSM_CRF02AG_nodrm.rds"){
  lst_fas <- prep_blast_align(infile)
  out <- align_mafft(lst_fas[['fas1']], lst_fas[['fas2']], lst_fas[['out']])
  return(out)
}

#
#
# ###-------------------------------#
# ###--- Start ucsd_hivclust ---###
# #-------------------------------#
# ##----------------------------##
# ###--- Apply UCSD
# ###--- HIV clustering
# ###--- to a dataframe of
# ###--- pairwise distances
# ##----------------------------##
#
# ###- input: path to RDS file of edge list
# ###- output: cluster assignements
# ###--- loops through distance thresholds ( = vector of threshold based on quantiles)
#
# ucsd_hivclust <- function(path.el, thr = NA, k = NA, out = "", quant = c(5e-5, 1e-4, 5e-4, 1e-3, 1e-2, 1e-1) ){
#
#   ## read RDS file
#   el <- readRDS(file = path.el)
#
#   ## get rid of distance > k
#   if(!is.na(k)) {
#     el <- el[el$distance < k,]
#     print(paste("restricting to distances <", k , sep = ''))
#   }
#
#   ## get var.name for output path (Regex !)
#   var.name <- paste("d", substr(path.el,
#                                 regexpr("\\/[^\\/]*$", path.el)[[1]][1] +1,
#                                 regexpr("\\_el.rds", path.el)[[1]][1] - 1 ),
#                     sep = '')
#
#   ## Either thr given or based on quantiles
#   if(is.na(thr[1])){
#
#     # choose threshold based on quantiles (first 4 for now)
#     qt <- quantile(el$distance,
#                    probs = quant )
#     thr <- signif(qt[1:4]*2, 1)/2 ## ensure 0.010, 0.015, etc.
#     print("thresholds based on quantiles")
#
#     ## get rid of distance > larger quantile
#     #      k <- round(qt[length(qt)], 2)
#     #      subel <- el[el$distance < k,]
#
#   } else {
#     thr <-  thr
#     qt <- "Not used"
#   }
#
#   ## write csv without rownames and get input path
#   inputCSV <- tempfile(pattern = "input", fileext = ".csv")
#   write.csv(el, file = inputCSV, row.names = FALSE )
#
#   ## full path needed
#   exec <- '~/Documents/softwares/hivclustering/scripts/hivnetworkcsv'
#
#   ####---- loop threshold (first 4 qt = 0.05, 0.1, 1, 10)
#
#   ## empty results
#   cmd <- vector( mode= "character" )
#   warn <- list()
#
#   for ( t in thr ){
#     print(paste("threshold =", t))
#     ## output
#     outputCSV <- paste(out, var.name, "_ucsd_hivclust_output_",
#                        t, ".csv", sep = '')
#
#     ## parms
#     parms <- paste("-t", t, "-f plain")
#
#     ## command
#     cmd_hivclustering <- paste(exec, "-i", inputCSV, "-c", outputCSV, parms )
#
#     print(cmd_hivclustering)
#
#     ## too long, the fitting of degree
#     ## distribution takes most time
#     ## => just issue command
#     ## to run on terminal
#     if ( !file.exists(outputCSV) ){
#       stderr <-  system(
#         paste(cmd_hivclustering, "2>&1"),
#         intern = TRUE)
#       print(paste("Writing", outputCSV ))
#
#       # save commands and 'stderr' warnings
#       cmd <- c(cmd, cmd_hivclustering)
#       warn <- c(warn, c(t, stderr))
#     } else print(paste("output already exists:", outputCSV))
#   }
#
#   ###--- bin table in one list
#   cl <- list()
#   for(i in 1:length(thr)){
#     ## add table i
#     CSV <- paste(out, var.name,
#                  "_ucsd_hivclust_output_", thr[i],
#                  ".csv", sep = '')
#     if(file.exists(CSV)){
#       cl[[i]] <- read.csv(CSV)
#       # name threshold
#       names(cl)[i] <- thr[i]
#     }
#   }
#
#   return(list(qt = qt, cmd = cmd, warn = warn, cl = cl))
# }
# #-----------------------------#
# ###--- End ucsd_hivclust ---###
# #-----------------------------#
#
# #--------------------------------#
# ###--- Start TreeToEdgeList ---###
# #--------------------------------#
# ##-----------------------------------------##
# ####---- function converting timed tree
# ####---- in edge list of pairwise distances
# ##-----------------------------------------##
#
# ###--- input: timed tree, substitution rate
# ###--- output: returns path for dataframe (ID1, ID2, distance)
#
# ###--- Scales the time based distance by consensus mutation rate to obtain 'standard' substitution per site distances
# # tree.name <- "simtree" tree.name <- "uktree"
# TreeToEdgeList <- function(t, name.output = NA, rate = 1, seqlength = NA,
#                            output = "data/", fig = "figure/" ,
#                            stats = TRUE, plot = TRUE ){
#
#   ###--- object names for output path
#   if (!is.na(name.output)){
#     tree.name <- name.output
#   } else tree.name <- substitute(t)
#   dist.mat.name <- paste(tree.name, "_dist", sep = '')
#   dist.mat.path <- paste(output, dist.mat.name,".rds", sep = '')
#   plot.name <- paste(fig, dist.mat.name, ".png", sep ='')
#   out_edge_list <- paste(output, tree.name, "_el.rds", sep = '')
#
#   ## if file exists, don't run and return the path
#   if (file.exists(out_edge_list) ){
#     print(paste('Edge list already exists here:', out_edge_list))
#     return(out_edge_list)
#   } else {
#
#     ###--- for time-based tree,
#     ##- estimate subst along branches
#     if (rate != 1){
#       print(paste('distances scaled by rate =', rate, 'subst/site/day'))
#       t$edge.length <- rpois(length(t$edge.length),
#                              t$edge.length * rate * seqlength ) / seqlength ## (random number of subst)/ number sites
#     }
#
#     ###--- get distances between tips
#     if (file.exists(dist.mat.path)){
#       d <- readRDS( dist.mat.path )
#     } else {
#       d <- as.dist(cophenetic.phylo(t))
#       print("Matrix of cophenetic distances was saved")
#       saveRDS(d, file = dist.mat.path )
#     }
#
#     ###--- stats
#     if(stats == TRUE){ print(summary(d)) }
#     if(plot == TRUE){
#       print(paste("png plot saved:", plot.name))
#       png(filename = plot.name, type="cairo",
#           units="in", width=5, height=4,
#           pointsize=12, res=96)
#       hist(d, breaks = 50, xlim = c(0, 1),
#            xlab = "distance", ylab = "frequency",
#            main = '', #tree.name # Tree's distances in subst/site
#            col = "grey")
#       dev.off()
#     }
#
#     ###--- normalize (if not normalized)
#     if(max(d) > 1) {
#       d <- round(d / max(d), 4)
#     }
#
#     ###--- as matrix
#     m <- as.matrix(d)
#
#     ###--- keep only the lower triangle by
#     ## filling upper with NA
#     m[upper.tri(m, diag=TRUE)] <- NA
#
#     ###--- create edge list if not there
#     require(reshape2)
#     if (file.exists(out_edge_list)){
#       el <- readRDS(out_edge_list)
#     } else {
#       # system.time(
#       el <- melt(m, na.rm = TRUE)
#       # ) # the na.rm removal takes most of time
#       colnames(el) <- c('ID1', 'ID2', 'distance')
#       saveRDS(el, file = out_edge_list)
#       print(paste("Edge list was saved:", out_edge_list))
#     }
#     return(out_edge_list)
#   }
# }
# #------------------------------#
# ###--- End TreeToEdgeList ---###
# #------------------------------#
#
# ##---- find.tree.outliers ----
# find.tree.outliers <- function(name_mltree){
#
#   print(basename(name_mltree))
#   ##- read ML tree
#   tr <- read.tree(name_mltree)
#
#   ##- distinguish outgroup OG
#   og <- tr$tip.label[grep("[A-Z]", tr$tip.label )]
#
#   ##- distances
#   d <- as.dist(cophenetic.phylo(tr))
#   print( summary(d) )
#
#   ##- detect outlier:
#   ## index of greater sum of distances
#   ## select uk seqs with sum(dist) > outgroup sequences !!
#
#   mat <- as.matrix(d, diag = FALSE)
#   su <- rowSums(mat)
#   # head( su[order(su, decreasing = TRUE)], 20 )
#
#   # compare with max sum(dist) for outgroup
#   lim <- max(su[names(su) %in% og])
#   outliers <- names(which(su > lim) )
#   if( length(outliers) == 0){
#     print('no outlier')
#   } else {
#     print(paste(length(outliers), 'outliers'))
#     return (outliers)
#   }
# }
# ##---- end find.tree.outliers ----#
#
# ##--- start rm.tree.outliers ----
# rm.tree.outliers <- function(name_mltree){
#   outliers <- find.tree.outliers(name_mltree)
#   ##- read ML tree
#   tr <- read.tree(name_mltree)
#   ##- tree without outliers
#   return(drop.tip(tr, outliers))
# }
# ##--- end rm.tree.outliers ----#
#
# ##--- start find.tn93.outliers ----
# ## list of outliers UK and outliers LANL OUTGROUP
# find.tn93.outliers <- function(msa){
#   # debug: msa = seq
#
#   ## length(msa)
#   nseq <- dim(msa)[1]
#   print(msa)
#
#   ##- distinguish outgroup OG
#   og <- rownames(msa)[grep("[A-Z]", rownames(msa))]
#
#   ##- distances
#   el <- ucsd_tn93(msa)
#   print( summary(el$Distance) )
#
#   ##--- outliers UK
#   ##- detect outlier: rule 1
#   ## index of greater sum of distances
#   ## select uk seqs with sum(dist) > outgroup sequences
#   su <- EdgelistToSumDist(el)
#   # head( su[order(su, decreasing = TRUE)], 30 )
#
#   # compare with max sum(dist) for outgroup
#   lim <- max(su[names(su) %in% og])
#   outliers <- names(which(su > lim) )
#
#   ##- detect outlier: rule 2
#   ## add rownames not in su because thr cannot exceed 1
#   notin <- setdiff(rownames(msa), names(su)) # or rownames(msa)[!(rownames(msa) %in% names(su))]
#   outliers <- c(outliers, notin)
#
#   ##- detect outlier: rule 3
#   ## tip name still in edge list but number of connections truncated by distance > 1
#   ## outlier if number connection >= that of outgroup
#   ##- Add connections from both from and to
#   ls <- list(tapply(el$Distance, el$ID1, length), tapply(el$Distance, el$ID2, length))
#   ##- unique from/to or ID1/ID2
#   colnames <- unique(unlist(lapply(ls, names)))
#   ##- list of same size/names vectors
#   mm <- matrix( unlist( lapply(ls, function(x) { x[colnames] }) ),
#                 nrow = length(ls), byrow = TRUE)
#   ##- number of connections from and to, for same node
#   links <- setNames( colSums(mm, na.rm = TRUE), colnames )
#   # head( links[order(links)])
#   lim.links <- min(links[names(links) %in% og])
#   outliers3 <- names(which(links < lim.links & !(names(links) %in% og) ) )
#   outliers <- unique(c(outliers, outliers3)) # final outliers UK vector
#
#   ##---- outliers outgroup
#   lim.og <- max(su[!(names(su) %in% og) & !(names(su) %in% outliers)])
#   outliers.og <- names(su[(names(su) %in% og) & su > lim.og])
#
#   if( length(outliers) == 0){
#     print('no outlier UK')
#   } else {
#     print(paste(length(outliers), 'outliers UK'))
#   }
#   print(paste(length(outliers.og), 'outliers OUTGROUP'))
#   return (list(uk = outliers, outgroup = outliers.og) )
# }
# ##--- end find.tn93.outliers ----#
#
# ##--- start rm.tn93.outliers ----
# rm.tn93.outliers <- function(name_mltree, name_msa){
#   load(name_msa) # seq
#   outliers <- find.tn93.outliers(seq)
#   print(outliers)
#   total.outliers <- do.call(c, outliers)
#   ##- read ML tree
#   tr <- read.tree(name_mltree)
#   ##- tree without outliers
#   return(drop.tip(tr, total.outliers))
# }
# ##--- end rm.tn93.outliers ----#
#
# MatrixToEdgelist <- function(m){
#   # for symmetric distance matrix
#   # debug m <- d_ape
#   require(reshape2)
#   if(class(m) != 'matrix') stop ('need a matrix')
#   ##- get upper triangle
#   m[ lower.tri(m, diag = FALSE) ] <- 0
#   ##- melt
#   el <-  melt(m)
#   ##- get rid of zeros
#   el <-  el[el$value != 0,]
#   rownames(el) <- NULL
#   return(el)
#   # names(el) <- c('ID1', 'ID2', 'Distance')
# }
#
# ucsd_tn93 <- function(msa, thr = 1){
#   # msa = rnd_seq
#   require(reshape2)
#   require(data.table)
#   require(ape)
#   FAS <- tempfile(fileext = '.fas')
#   OUT <- tempfile(fileext = '.csv')
#   write.dna(msa, file = FAS, format = 'fasta')
#   exec <- "tn93"
#   if(is.null(thr)){
#     parms <- "-q"
#   } else {
#     parms <- paste("-q -t", thr) # -q to track progress
#   }
#   args <-  paste(parms, "-o", OUT, FAS)
#   print(paste(exec, args))
#
#   stm <- system.time(
#     system2(exec, args)
#   )
#   print('Internal time')
#   print(stm)
#   dist <- fread(OUT, colClasses=list(character=1:2)) # faster than ( read.csv(OUT) )
#   return(as.data.frame(dist))
# }
#
# ape_tn93 <- function(msa, thr = 1){
#   # debug msa = rnd_seq
#   require(ape)
#   # require(reshape2)
#   d_ape <- dist.dna(msa, model = 'TN93', pairwise.deletion = TRUE, as.matrix = TRUE)
#   el <- MatrixToEdgelist(d_ape)
#   names(el) <- c('ID1', 'ID2', 'Distance')
#   return(el)
# }
#
# EdgelistToSumDist <- function(el){
#   ##- Add distances from both from and to
#   ls <- list(tapply(el$Distance, el$ID1, sum), tapply(el$Distance, el$ID2, sum))
#   ##- unique from/to or ID1/ID2
#   colnames <- unique(unlist(lapply(ls, names)))
#   ##- list of same size/names vectors
#   mm <- matrix( unlist( lapply(ls, function(x) { x[colnames] }) ),
#                 nrow = length(ls), byrow = TRUE)
#   ##- sum of distances from + to for same node
#   ds <- setNames( colSums(mm, na.rm = TRUE), colnames )
#   return(ds)
# }
#

##---- plot.tree.col.og ----
#' Plot tree with external sequences colored
#' @param tr phylo tree
#' @param tmrca date of root to add time axis
#' @param title title of plot
#' @param tiplab logical to show tip labels
#' @export
plot_tree_col_og <- function(tr, tmrca = NULL, title = NULL, tiplab = FALSE){
  tr2 <- ape::ladderize(tr)
  l <- tr2$tip.label
  wh <- ape::which.edge( tr2, grep("[A-Z].*", l ) ) ## assuming names outgroup contain letters
  ecolo <- rep('black', ape::Nedge(tr2))
  ecolo[wh] <- 'red'
  tcolo <- rep('black', ape::Ntip(tr2))
  tcolo[grep("[A-Z].*", l )] <- 'red'
  ape::plot.phylo(tr2, show.tip.label = tiplab, tip.color = tcolo, cex = .6, label.offset = min(tr2$edge.length)/4, edge.color = ecolo, edge.width = 2, main = title)
  if( !(is.null(tmrca)) ) ape::axisPhylo(root.time = tmrca, backward = F)
}

##---- grubbs.flag ----
#' Iterative grubbs outliers searching
#' @param x named vector of values
#' @return data.frame of value with Outlier logical

grubbs.flag <- function(x) {
  # require(outliers)
  if(is.null(names(x))){ x <- setNames(x, 1:length(x)) } # force named vector
  outliers <- NULL
  test <- x
  grubbs.result <- outliers::grubbs.test(test)
  pv <- grubbs.result$p.value # named
  while(pv < 0.05) {
    outliers <- c(outliers, names(pv))
    test <- x[!names(x) %in% outliers]
    grubbs.result <- outliers::grubbs.test(test)
    pv <- grubbs.result$p.value
  }
  return(data.frame(X = x, Outlier=(names(x) %in% outliers)))
}

##---- find.rtt.outliers ----
#' Root-to-tip regression outliers
#'
#' Superseded by \code{find.rtt.outliers2}, with cpp distance calculation
#' @param tree phylo tree
#' @param dates with tip names
#' @param mult multiplier of sd to define outliers
#' @param plots plot
#' @param title title of plot
#' @return output list of outliers from mult*SD and from grubbs
#' @importFrom graphics plot
#' @importFrom stats sd
#' @export
find.rtt.outliers <- function(tree, dates, mult = 2, plots = TRUE, title = NULL){
  ##- find root
  # require(parallel)
  nc <- parallel::detectCores()
  system.time(
    t_rooted <- ape::rtt(tree, dates, objective = 'correlation', ncpu = nc - 1)
  )
  ##- regression between RTT distance and dates
  ds <- ape::dist.nodes(t_rooted)
  rttd <- ds[ape::Ntip(t_rooted) + 1, 1:ape::Ntip(t_rooted)] ## Root-to-tip distances
  # This assumes that the root node is the first node in the numbered list 1:Nnode,
  #  *after* the tip nodes, which are labeled 1:ntips
  names(rttd) <- tree$tip.label

  ##- regression
  y <- rttd[order(names(rttd))]
  x <- dates[order(names(dates))]
  m1 <- stats::lm(y ~ x)  # linear model

  ##- rate and tMRCA
  rate <- round(m1$coef[2], 5)
  tmrca <- round(- m1$coef[1] / m1$coef[2], 2)
  print(paste('rate', rate, 'tMRCA', tmrca))
  ## same as TempEst if both with rms calculation

  ##- residuals
  rs <- stats::residuals(m1)
  # head(rs[order(abs(rs), decreasing = TRUE)])

  ##- outliers based on standard deviation
  ## mult = 2 standard deviations
  sd2 <- mult * sd(rs)
  outs <- ifelse(abs(rs) > sd2, 1, 0)
  sum(outs)
  outliers.sd <- names(x)[as.logical(outs)]
  # outliers.sd
  if (plots) plot(rs, col = outs + 1, main = title)
  ##- outliers based on grubbs algo
  liste <- grubbs.flag(x = rs)
  out.grubbs <- rownames(liste[liste$Outlier,])
  if (plots) plot(liste$X, col = liste$Outlier + 1, main = title)
  return(list(SD=outliers.sd, Grubbs=out.grubbs))
}
##-- end find.rtt.outliers --

##---- get.lsd.metrics -----
####- extract rate amd tMRCA
## depends on LSD version
get.lsd.metrics <- function(txt, ci = NULL){
  string <- readLines(con = txt)
  r <- string[grep('tMRCA', string)]
  if(length(r) == 1){
    nums <- gregexpr('[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?', r) # regex scientific notation: http://www.regular-expressions.info/floatingpoint.html
    words <- gregexpr('([A-Za-z]|\\s|_){2,}', r) # words of at least two letters or space or underscore
    result <- setNames( as.numeric(unlist(regmatches(r, nums))),
                        trimws(unlist(regmatches(r, words))) ) # trim spaces
  } else {
    r.last <-  r[length(r)]
    words <- gregexpr('([A-Za-z]|\\s|_){2,}', r.last)
    w <- trimws(unlist(regmatches(r.last, words)))
    nums <- gregexpr('[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?', r.last)
    n <- as.numeric(unlist(regmatches(r.last, nums)))
    if(is.null(ci)){
      result <- setNames(n, w)
    } else {
      result <- setNames(n, c(w[1], paste0(w[1],'.lo'), paste0(w[1],'.up'), w[2], paste0(w[2],'.lo'), paste0(w[2],'.up'), w[3]))
    }
  }
  return(result)
}


##---- lsd.examl.tree ----
## debug:
# TREE="data/ML_TREE0/ExaML_result.MSM_A_nodrm_ref_og_fltd.finaltree.000"; dates = NA; rooting = TRUE; outgroup = TRUE; ci = '-f 10'; lsd.version = '0.3'; verbose = TRUE

#' Date tree with LSD
#' @param TREE Path to tree
#' @param dates Dates of tips [default = NA]. By default found by \code{dates_of_tree}
#' @param rooting Logical indicating that rooting is done on outgroup of sequences (passed to \code{dates_of_rtree})
#' @param outgroup Outgroup option of LSD (LSD outputs a tree without outgroup)
#' @param ci option for confidence interval
#' @param lsd.version [default = 0.3 (beta)]
#' @param verbose Details
#' @return NULL. Save list of [[1]] metrics and [[2]] dated tree
#' @export
lsd.examl.tree <- function(TREE,
                           dates = NA,
                           rooting = TRUE,
                           outgroup = TRUE,
                           ci = '-f 10', #NULL,
                           lsd.version = '0.3',
                           verbose = TRUE){

  ##- names
  PATHTREE <- dirname(TREE)
  NAMETREE <- basename(TREE)
  INPUTDAT <-  tempfile(pattern = "LSD_dates_", fileext = ".txt")
  INPUTNWK <- tempfile(pattern = "LSD_tree_", fileext = ".nwk")
  OUTTREE <- paste0(NAMETREE, '_dated_tree')
  OUTMETRICS <- paste0(NAMETREE, '_lsd_metrics')
  ##- import tree and find dates
  if(is.na(dates)){
    lst <- dates_of_rtree(TREE, rooting = rooting)
    tr <- lst$tree
    dates <- lst$dates
    MAIN <- lst$group
  } else {
    .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", NAMETREE)
    MAIN <- regmatches(NAMETREE, .m) # name of group
    tr <- ape::read.tree(TREE)
  }
  if(length(dates)!=ape::Ntip(tr)) warning("Not the same number of tree tips and dates")
  ## write dates for LSD   # "The first line must indicate the number of leaves"
  ## don't need dates of outgroup in LSD
  if(outgroup){
    INPUTOG <- get_outgroup_for_lsd(TREE, verbose = verbose)
    names.outgroup <- as.vector(read.table(INPUTOG)[-1,1])
    dates <- dates[!names(dates) %in% names.outgroup]
    }
  write.table( cbind( names(dates), dates), file = INPUTDAT, col.names=c(as.character(length(dates)), '') , quote = F, row.names=F)
  if (verbose) print(paste(length(dates), "dates as LSD input"))
  ## create newick tree
  ape::write.tree(tr, file = INPUTNWK)
  ## command by version
  ## default
  if(!(lsd.version %in% c('0.1', '0.2', '0.3'))) stop('LSD version not found')
  arg0 <- paste("-i", INPUTNWK, "-d", INPUTDAT)
  if( any(grep("Darwin", Sys.info())) ){ # Mac
    lsd.cmd <- setNames(
      c('lsd',
        '~/Documents/softwares/lsd-0.2.2/src/lsd',
        '~/Documents/softwares/lsd-0.3beta/src/lsd'),
      c('0.1', '0.2', '0.3'))
  } else { # HPC
    lsd.cmd <- setNames(
      c('lsd',
        '~/Documents/softwares/lsd-0.2.2/src/lsd',
        '/home/slevu/lsd-0.3beta/src/lsd'),
      c('0.1', '0.2', '0.3'))
  }
  if(outgroup){
    lsd.args <- setNames(
      c(paste(arg0, "-g", INPUTOG, "-c -s 1200 -v"),
        paste(arg0, "-g", INPUTOG, "-c -s 1200 -v"),
        paste(arg0, "-g", INPUTOG, "-c -s 1200 -v 2", ci)),
      c('0.1', '0.2', '0.3'))
  } else {
    lsd.args <- setNames(
      c(paste(arg0, "-r -c -s 1200 -v"),
        paste(arg0, "-r a -c -s 1200 -v"),
        paste(arg0, "-r a -c -s 1200 -v 2", ci)),
      c('0.1', '0.2', '0.3'))
  }
  cmd <- path.expand(lsd.cmd[ lsd.version ])
  parms <- path.expand(lsd.args[ lsd.version ])
  ##- run
  if (verbose) print(paste('LSD', lsd.version))
  if (verbose) print(paste(cmd, parms))
  system2(cmd, args = parms)
  ##- rename and copy results (depends on LSD version)
  if(lsd.version == '0.1'){
    tmp <- sub(pattern = '.nwk', '', basename(INPUTNWK)) # basename of tmp tree to be replaced
  } else {
    tmp <- basename(INPUTNWK) # basename of tmp tree to be replaced
  }
  results <- list.files(tempdir(), full.names = TRUE)
  ##- files to be copied # version independent
  from <- results[grepl(tmp, results) & grepl('result|ingroup', results)]
  ##- create dir for results
  dir.create(paste(PATHTREE, 'lsd', sep = '/'), showWarnings = FALSE)
  if (verbose) print("create `lsd` folder if does not exist")
  ##- strip nameTREE
  # sub('ExaML_result.sub_uk_', '', x = nameTREE)
  to <- gsub(tmp, paste(PATHTREE, 'lsd', NAMETREE, sep ='/'), basename(from) )
  file.copy(from, to, overwrite = TRUE)
  if (verbose) print(paste0('copy result files to ', PATHTREE, '/lsd'))
  # file.remove( list.files(tempdir(), full.names = TRUE) )
  ##- output rates and dated tree
  metrics <- get.lsd.metrics(txt = to[ !(grepl('newick|nexus|ingroup', to)) ], ci = ci) # get txt file, same as `grep('newick|nexus', to, invert = TRUE)`
  if (verbose) print(metrics)
  ##- import dated tree
  dated_tree <- ape::read.tree( to[ grepl('date', to) & grepl('newick', to) ] )
  ##- output
  out <- list(metrics, dated_tree )
  names(out) <- c(OUTMETRICS, OUTTREE)
  nm.outfile <- paste0(TREE, '_lsd.rds')
  print(paste("save", nm.outfile))
  saveRDS(out, file = nm.outfile)
}

##---- find.distant.tips ----
##- check ancestral tips
##- find tips with greater terminal edges distance
find.distant.tips <- function(tr){
  # debug: tr = dated_tree

  ## number of tips
  n.tips <- length(tr$tip.label) # Ntip
  ## select terminal edges
  selec <- tr$edge[,1] <= n.tips | tr$edge[,2] <= n.tips
  ## length of terminal edges
  el <- tr$edge.length[selec]
  ## order of edge length
  ord <- order(el, decreasing = TRUE)
  ## ordered named vector of lengths
  v <- setNames(el[ord], tr$tip.label[ord] )
  return( v )
}

##---- outliers.lsd ----
outliers.lsd <- function(lst = 'lsd_results_no_outliers_G'){
  list.outliers <- vector("list", length(lst))
  # debug: i = 1
  for(i in 1:length(lst)){
    lsd_out <- readRDS( lst[i] )

    ##- extract info
    .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", names(lsd_out)[1])
    MAIN <- regmatches(names(lsd_out)[1], .m) # name of group
    e <- find.distant.tips(lsd_out[[2]])
    k <- which.max( grepl(pattern = "[A-Z]", names(e)) ) # first outgroup
    list.outliers[i] <- list(names(e)[1:(k - 1)])
    names(list.outliers)[i] <- MAIN
  }
  return(list.outliers)
}


##---- find_mrca_outliers ----
#' @title Find tips that should not be ancestral
#' @description Find tips producing nodes ancestral to MRCA from global external sequences
#' @param tree dated tree
# #' @inheritParams lsd.examl.tree
#' @return problematic tip labels
#' @export
find_mrca_outliers <- function(tree){
  ##- filter
  #- find mrca of non-uk, if labeled with LETTERS
  m <- ape::getMRCA(tree, tip = tree$tip.label[grepl("[A-Z]", tree$tip.label)])
  # if m is root, do nothing
  if(m == ape::Ntip(tree)+1){
    print('Root at MRCA of global, nothing to return')
    return(NULL)
  } else {
    #- tips descendants of m
    dm <- phytools::getDescendants(tree, m)[phytools::getDescendants(tree, m) <= ape::Ntip(tree)]
    #- tips remaining
    x <- tree$tip.label[-dm]
    print(paste(length(x), "tips producing node ancestral to global MRCA"))
    return(x)
  }
}


##---- metrics.lsd ----
#' Put LSD metrics from list in a table
#' @param lst List of LSD results
#' @export
metrics.lsd <- function(lst = "lsd_results"){
  l <- length(readRDS( lst[1] )[[1]]) # read first for length results vector
  lsd_metrics <- data.frame(matrix(0, nrow = length(lst), ncol = l + 2 ))
  # i = 1
  for(i in 1:length(lst)){
    lsd_out <- readRDS( lst[i] )
    ##- extract info
    .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", names(lsd_out)[1])
    MAIN <- regmatches(names(lsd_out)[1], .m) # name of group
    lsd_metrics[i, 1] <- MAIN
    .n <- regexpr("[0-9]{3}", names(lsd_out)[1])
    BSN <- as.numeric(regmatches(names(lsd_out)[1], .n)) # number BS
    lsd_metrics[i, 2] <- BSN
    lsd_metrics[i, 3:(l+2)] <- signif(lsd_out[[1]], 4)
  }
  colnames(lsd_metrics) <- c('group', 'bs', names(lsd_out[[1]]))
  return(lsd_metrics)
}

##---- plot_all_lsd_tree ----
#' from list of LSD results, plot all dated trees at once
#' @param lst list of dated trees
#' @export
plot_all_lsd_tree <- function(lst = "lsd_results_og_nooout"){
  graphics::par(mfrow = c(2, ceiling(length(lst)/2) ),
      mar = c(4,3,3,2)+0.1, oma = c(0, 0, 2, 0), bty = 'n') # b,l,t,r
  for(i in 1:length(lst)){
    lsd_out <- readRDS( lst[i] )
    x <- names(lsd_out)[1]
    ##- extract info
    .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", x)
    MAIN <- regmatches(x, .m) # name of group
    N <- sub(".*(\\d{3}).*", "\\1", x) # number bootstrap
    metrics <- lsd_out[[1]]
    dated_tree <- lsd_out[[2]]
    ##- plot
    plot_tree_col_og(dated_tree, tmrca = metrics['tMRCA'], title = paste(MAIN, N))
  }
  # mtext(paste(depvar, "(y) vs co-variates (x) "), outer = TRUE, cex = 1)
  # dev.off()
}

##---- lsd.probs ----
#' @title Find ancient nodes with uk only descendants
#' @description having dates of each node
#' @param lst List of LSD results
#' @export
lsd.probs <- function(lst = "lsd_results_og_nooout"){
  ##- initialize
  list.probs <- vector("list", length(lst))
  for(i in 1:length(lst)){
    lsd_out <- readRDS( lst[i] )
    ##- extract info
    .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", names(lsd_out)[1])
    MAIN <- regmatches(names(lsd_out)[1], .m) # name of group
    # metrics <- lsd_out[[1]]
    dated_tree <- lsd_out[[2]]
    list.probs[i] <- list(find_mrca_outliers(dated_tree))
    names(list.probs)[i] <- MAIN
  }
  return(list.probs)
}

##---- find.rtt.outliers2 ----
#' @title C++ version of Root-to-tip regression outliers
#' @description with dist2root.r
#' @inheritParams find.rtt.outliers
#' @param verbose explain
#' @export
find.rtt.outliers2 <- function(tree, dates, mult = 2, title = NULL, verbose = TRUE, plots = FALSE){
  #source( 'dist2root.R' )
  rttd <- dist2root( tree )
  ##- regression
  y <- rttd[order(names(rttd))]
  x <- dates[order(names(dates))]
  if (plots) plot(x,y, main = "RTT dist vs dates")
  m1 <- stats::lm(y ~ x)  # linear model
  ##- rate and tMRCA
  rate <- round(m1$coef[2], 5)
  tmrca <- round(- m1$coef[1] / m1$coef[2], 2)
  if (verbose) print(paste('Parameters from RTT lm: rate', rate, 'tMRCA', tmrca))
  ##- residuals
  rs <- stats::residuals(m1)
  # head(rs[order(abs(rs), decreasing = TRUE)])
  ##- outliers based on standard deviation
  ## mult = 2 standard deviations
  sd2 <- mult * sd(rs)
  outs <- ifelse(abs(rs) > sd2, 1, 0)
  if (verbose) print(paste("Number outside", mult, "* sd:", sum(outs)))
  out.sd <- names(x)[as.logical(outs)]
  if(plots) {
    if(sum(outs) != 0L){
    plot(rs, col = outs + 1, main = title)
    }
  }
  ##- outliers based on grubbs algo
  liste <- grubbs.flag(x = rs)
  out.grubbs <- rownames(liste[liste$Outlier,])
  if (verbose) print(paste("Number Grubbs outliers:", length(out.grubbs)))
  if(plots){
    if(length(out.grubbs) != 0L){
    plot(liste$X, col = liste$Outlier + 1, main = title)
    }
  }
  lst <- list(SD=out.sd, Grubbs=out.grubbs)
  return(lst)
}

##---- dates_of_tree ----
#' Find dates of tree for UK and LANL sequences
#' @param ALLDATES_LANL Path to dates file from LANL sequences
#' @param ALLDATES_UK Path to dates from UK sequences
#' @param rooting rooting on outgroup
#' @inheritParams lsd.examl.tree
#' @return List of (tree, dates, name of tree, number tips)
dates_of_rtree <- function(TREE,
                           ALLDATES_LANL = get_lanl_dates(),
                           ALLDATES_UK = get_uk_dates(),
                           rooting = FALSE,
                           verbose = TRUE
){
  ##- names
  NAMETREE <- basename(TREE)
  .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", NAMETREE)
  MAIN <- regmatches(NAMETREE, .m) # name of group
  ##- import
  tr <- ape::read.tree(TREE)
  if (verbose) print(paste("Initial", MAIN, "tree has", ape::Ntip(tr), "tips"))
  ##- import dates
  dates_uk <- readRDS(ALLDATES_UK)
  dates_og  <- readRDS(ALLDATES_LANL)
  ## bind dates
  dates <-  c(dates_uk[names(dates_uk) %in% tr$tip.label ],
              dates_og[names(dates_og) %in% tr$tip.label & !(is.na(dates_og))] ) ## without NA
  dates <- dates[unique(names(dates))]
  ##- save tree with non null dates
  tr <- ape::drop.tip(tr, tip = tr$tip.label[ !(tr$tip.label %in% names(dates)) ])
  n <- ape::Ntip(tr)
  if (verbose) print(paste("Tree has", n, "tips with a date"))
  if(rooting){
    INPUTOG <- get_outgroup_for_lsd(TREE, verbose = verbose)
    names.outgroup <- as.vector(read.table(INPUTOG)[-1,1])
    if(!all(names.outgroup %in% tr$tip.label)) warning("Outgroup not found in tip labels")
    tr <-  ape::root(tr, names.outgroup, resolve.root = TRUE)
  } else {
    names.outgroup = NA
  }
  return( list(tree = tr, dates = dates, group = MAIN, names.outgroup = names.outgroup, ntips = n) )
}

##---- rm_rtt_outliers ----
#' Remove RTT outliers from tree
#' @param OUTFILE Path to save list of outliers, in rds
#' @param rooting Passed to dates_of_rtree
#' @inheritParams lsd.examl.tree
#' @return Path to filtered tree
#' @export
rm_rtt_outliers <- function(TREE, outgroup = TRUE,
                            OUTFILE = NULL,
                            rooting = TRUE,
                            verbose = TRUE){
  lst <- dates_of_rtree(TREE, rooting = rooting)
  tr <- lst$tree
  dates <- lst$dates
  MAIN <- lst$group
  if (is.null(OUTFILE)) OUTFILE <- paste0(dirname(TREE), '/', MAIN, "_rtt.outliers.rds")
  NEWTREE <- sub("\\.finaltree", "_noRTTout\\.finaltree", TREE)
    ##- remove outliers
      if(outgroup){ # but don't remove the outgroup
      INPUTOG <- get_outgroup_for_lsd(TREE, verbose = verbose)
      names.outgroup <- as.vector(read.table(INPUTOG)[-1,1])
      tr0 <- ape::drop.tip(tr, tip = names.outgroup)
      dates0 <- dates[names(dates) %in% tr0$tip.label]
    } else {
      tr0 <- tr
      dates0 <- dates
    }
    outliers.rtt <- find.rtt.outliers2(tree = tr0, dates = dates0, title = MAIN)
    if (verbose){
      print(paste('remove', length(outliers.rtt[['Grubbs']]), 'Grubbs outliers from RTT regression:'))
      dput(outliers.rtt[['Grubbs']])
      print(paste("Outliers saved in", OUTFILE))
      saveRDS(outliers.rtt, file = OUTFILE)
    }
    ##- remove on tr outliers found on tr0
    tr <- ape::drop.tip(tr, outliers.rtt[['Grubbs']])
    #dates <- dates[names(dates) %in% tr$tip.label]
    if (verbose){
      print(paste("After RTT regression, tree has", ape::Ntip(tr), "tips"))
      print(paste("Tree saved in", NEWTREE))
    }
    ape::write.tree(tr, NEWTREE)
    return( NEWTREE )
}

##---- rm_mrca_outliers ----
#' Remove mrca outliers from tree
#' @param OUTFILE Path to save list of outliers, in rds
#' @inheritParams lsd.examl.tree
#' @return Path to filtered tree
#' @export
rm_mrca_outliers <- function(TREE, outgroup = TRUE,
                             OUTFILE = NULL, verbose = TRUE){
  ##- names
  NAMETREE <- basename(TREE)
  .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", NAMETREE)
  MAIN <- regmatches(NAMETREE, .m) # name of group
  if (is.null(OUTFILE)) OUTFILE <- paste0(dirname(TREE), '/', MAIN, "_mrca.outliers.rds")
  NEWTREE <- sub("\\.finaltree", "_noMRCAout\\.finaltree", TREE)
  ##- read
  tree <- ape::read.tree(TREE)
  ##- evaluate without outgroup
  if (outgroup) {
    INPUTOG <- get_outgroup_for_lsd(TREE, verbose = verbose)
    names.outgroup <- as.vector(read.table(INPUTOG)[-1,1])
    if(!all(names.outgroup %in% tree$tip.label)) warning("Outgroup not found in tip labels")
    tree <-  ape::root(tree, names.outgroup, resolve.root = TRUE)
    tree <- ape::drop.tip(tree, names.outgroup)
  }
  outliers.mrca <- find_mrca_outliers(tree = tree)
  if (verbose){
    print(paste('remove', length(outliers.mrca), 'outliers from mrca filtering:'))
    dput(outliers.mrca)
    print(paste("Outliers saved in", OUTFILE))
    saveRDS(outliers.mrca, file = OUTFILE)
  }
  ##- remove outliers
  tr <- ape::drop.tip(tree, outliers.mrca)
  if (verbose){
    print(paste("After mrca filtering, tree has", ape::Ntip(tr), "tips"))
    print(paste("Tree saved in", NEWTREE))
  }
  ape::write.tree(tr, NEWTREE)
  return(NEWTREE)
}

##---- filter_tree_outliers ----
#' Remove RTT and MRCA outliers from tree
#' @param TREE Path to tree
#' @param outgroup Logical indicating that rooting is on outgroup of sequences
#' @param OUTFILE Path to save list of outliers, in rds
#' @param rooting passed to dates_of_rtree
#' @param verbose stuff done
#' @return List of (group name, path to filtered tree, ntips original tree, ntips filtered tree, n rtt outliers, n mrca outliers, proportion tips filtered)
#' @export
filter_tree_outliers <- function(TREE,
                             outgroup = TRUE,
                             OUTFILE = NULL,
                             rooting = TRUE,
                             verbose = TRUE){
  ##- get dates
  lst <- dates_of_rtree(TREE, rooting = rooting)
  tr <- lst$tree
  dates <- lst$dates
  MAIN <- lst$group
  names.outgroup <- lst$names.outgroup
  ntips <- lst$ntips
  ##- names
  if (is.null(OUTFILE)) {
    OUTFILE <- paste0(dirname(TREE), '/', MAIN, "_fltd.outliers.rds") }
  NEWTREE <- sub("\\.finaltree", "_fltd\\.finaltree", TREE)
  ##- filter without outgroup
  if(outgroup){ # but don't remove the outgroup
    tr0 <- ape::drop.tip(tr, tip = names.outgroup)
    dates0 <- dates[names(dates) %in% tr0$tip.label]
  } else {
    tr0 <- tr
    dates0 <- dates
  }
  ##- RTT outliers
  outliers.rtt <- find.rtt.outliers2(tree = tr0, dates = dates0, title = MAIN)
  or <- length(outliers.rtt[['Grubbs']])
  if (verbose){
    print(paste('found', or, 'Grubbs outliers from RTT regression:'))
    dput(outliers.rtt[['Grubbs']])
  }
  ##- MRCA outliers
  outliers.mrca <- find_mrca_outliers(tree = tr0)
  om <- length(outliers.mrca)
  if (verbose){
    print(paste('found', om, 'outliers from mrca filtering:'))
    dput(outliers.mrca)
  }
  ##- remove on tr outliers found on tr0
  outliers <- c(outliers.rtt, mrca = list(outliers.mrca), filtered = list(outliers.rtt[['Grubbs']])) # list(unique(c(outliers.rtt[['Grubbs']], outliers.mrca))))
  tr <- ape::drop.tip(tr, outliers[['filtered']] )
  nf <- ape::Ntip(tr)
  if (verbose){
    print(paste("After filtering, tree has", nf, "tips"))
    print(paste("Tree saved in", NEWTREE))
    print(paste("Outliers saved in", OUTFILE))
  }
  saveRDS(outliers, file = OUTFILE)
  ape::write.tree(tr, NEWTREE)
  return(list(group = MAIN,
              filtered_tree = NEWTREE,
              filtered_outliers = outliers,
              counts = setNames(c(ntips, nf, or, om, signif((ntips - nf)/ntips*100, 2)), c("n_init_tree", "n_fltd_tree", "n_rtt_out", "n_mrca_out", "pct_fltd"))
              )
  )
}

##---- filter_tree_outliers2 ----
#' DO NOT USE: Remove RTT and MRCA outliers from tree + branch length outliers
#' @inheritParams filter_tree_outliers
# #' @export
filter_tree_outliers2 <- function(TREE,
                                 outgroup = TRUE,
                                 OUTFILE = NULL,
                                 rooting = TRUE,
                                 verbose = TRUE){
  ##- get dates
  lst <- dates_of_rtree(TREE, rooting = rooting)
  tr <- lst$tree
  dates <- lst$dates
  MAIN <- lst$group
  names.outgroup <- lst$names.outgroup
  ntips <- lst$ntips
  ##- names
  if (is.null(OUTFILE)) {
    OUTFILE <- paste0(dirname(TREE), '/', MAIN, "_fltd.outliers.rds") }
  NEWTREE <- sub("\\.finaltree", "_fltd\\.finaltree", TREE)
  ##- filter without outgroup
  if(outgroup){ # but don't remove the outgroup
    tr0 <- ape::drop.tip(tr, tip = names.outgroup)
    dates0 <- dates[names(dates) %in% tr0$tip.label]
  } else {
    tr0 <- tr
    dates0 <- dates
  }
  ##- RTT outliers
  outliers.rtt <- find.rtt.outliers2(tree = tr0, dates = dates0, title = MAIN)
  or <- length(outliers.rtt[['Grubbs']])
  if (verbose){
    print(paste('found', or, 'Grubbs outliers from RTT regression:'))
    dput(outliers.rtt[['Grubbs']])
  }
  ##- MRCA outliers
  outliers.mrca <- find_mrca_outliers(tree = tr0)
  om <- length(outliers.mrca)
  if (verbose){
    print(paste('found', om, 'outliers from mrca filtering:'))
    dput(outliers.mrca)
  }
  ##- branch length outliers
  e <- find.distant.tips(tr0)
  k <- which.max( grepl(pattern = "[A-Z]", names(e)) ) # first outgroup
  outliers.branch <- ifelse(k == 1, NULL, names(e)[1:(k - 1)])
  ob <- length(outliers.branch)
  if (verbose){
    print(paste('found', ob, 'outliers from branch length filtering:'))
    dput(outliers.branch)
  }
  ##- remove on tr outliers found on tr0
  outliers <- c(outliers.rtt, mrca = list(outliers.mrca), branch = list(outliers.branch), filtered = list(unique(c(outliers.rtt[['Grubbs']], outliers.mrca, outliers.branch))))
  tr <- ape::drop.tip(tr, outliers[['filtered']])
  nf <- ape::Ntip(tr)
  if (verbose){
    print(paste("After filtering, tree has", nf, "tips"))
    print(paste("Tree saved in", NEWTREE))
    print(paste("Outliers saved in", OUTFILE))
  }
  saveRDS(outliers, file = OUTFILE)
  ape::write.tree(tr, NEWTREE)
  return(list(group = MAIN,
              filtered_tree = NEWTREE,
              filtered_outliers = outliers[['filtered']],
              counts = setNames(c(ntips, nf, or, om, ob, signif((ntips - nf)/ntips*100, 2)), c("n_init_tree", "n_fltd_tree", "n_rtt_out", "n_mrca_out", "n_branch_out", "pct_fltd")) )
  )
}

##---- rm_outliers_msa ----
#' Remove outliers from MSA in DNAbin
#' @param pathmsa Path to MSA saved as \code{seq} DNAbin object
#' @param pathoutliers Path to outliers contained in labels of \code{seq}, vector saved in RDS format
#' @param suffix Suffix added to filename
#' @param verbose print details
#' @return save a fasta file
#' @export
rm_outliers_msa <- function(pathmsa, pathoutliers, suffix = "fltd", verbose = TRUE){
  NMOUT <- paste0( dirname(pathmsa), '/', sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(pathmsa)), "_", suffix, ".fas"  )
  seq <- load_msa(pathmsa)
  outs <- readRDS(pathoutliers)
  if(length(outs[['filtered']]) == 0){
    print(paste("No outlier, save copy to", NMOUT))
  } else {
    if(verbose) {
      print(dim(seq))
      print( paste('remove', length(outs[['filtered']]), 'outliers of', dim(seq)[1],
                   'sequences:', signif(length(outs[['filtered']]) / dim(seq)[1] *100, 3), '%' ) )
    }
    seq <- seq[!(rownames(seq) %in% outs[['filtered']]), ]
    if(verbose){
      print(dim(seq))
      print( paste('write', NMOUT) )
    }
  }
  ape::write.dna(seq, file = NMOUT, format = 'fasta')
  return( NMOUT )
}

##---- find.clades.0 ----
#' Find clades of state = 0 descending from state = 1
#' @param tree phylo tree
#' @param tip.state 0/1 for UK/non-UK
#' @param verbose details
#' @param plots plot colored tree
#' @export
find.clades.0 <- function(tree, tip.state, verbose = TRUE, plots = FALSE){
  #require(ape)
  #require(phytools)
  ## states of nodes
  nnodes <- ape::Ntip(tree) + ape::Nnode(tree)
  state <- rep(0, length = nnodes)
  wh <- ape::which.edge( tree, tree$tip.label[tip.state == 1] ) # edges relating state=1 to their MRCA
  state[ unique(as.vector(tree$edge[wh,])) ] <- 1 # mark nodes of state 1
  if(plots){
    colo <- rep('black', ape::Nedge(tree))
    colo[wh] <- 'red'
    ape::plot.phylo(tree, show.tip.label = FALSE, use.edge.length = FALSE, edge.color = colo, edge.width = 2)
  }
  ##- roots = from 1 to 0 node transition
  which.root <- which(state[tree$edge[, 1]] == 1 & state[tree$edge[, 2]] == 0)
  roots <- tree$edge[which.root, 2] ## with tips
  if(state[tree$edge[1,1]] == 0){ roots <- c(tree$edge[1,1], roots)} ## add 'real' root if 0
  real_roots <- roots[roots > ape::Ntip(tree)] # only nodes
  ##- list of clades
  clades <- list()
  for (i in 1:length(real_roots)){
    non_clade <- vector()
    for (j in setdiff(roots, real_roots)){
      non_clade <- c(non_clade, tree$tip.label[phytools::getDescendants(tree, j)] )
    }
    clade <- ape::extract.clade(tree, real_roots[i]) # remove sub clades or tips
    clade <- ape::drop.tip(clade, tree$tip.label[which(state[1:length(tree$tip.label)] == 1)]) # remove state = 1
    clades[[i]] <- ape::drop.tip(clade, non_clade)
  }
  class(clades) <- "multiPhylo"
  if (verbose) print(paste('number of clades', length(clades)))
  return(clades)
}

##---- get_clades ----
#' Exclude treated patients and find clades
#' @param lsd_result LSD result as list from which [[2]] = dated tree
#' @param not_naive vector of non naive test index
#' @param plots Plot of colored tree(passed to \code{find.clades.0})
#' @param verbose details
#' @export
get_clades <- function(lsd_result, not_naive = NULL, plots = FALSE, verbose = TRUE){
  lsd_out <- readRDS( lsd_result )
  dated_tree <- lsd_out[[2]]
  if(!is.null(not_naive)){
    ##- exclude not naive UK sequences
    if(verbose) print(paste("excluding", sum(dated_tree$tip.label %in% not_naive), "experienced patients' sequences out of", length(dated_tree$tip.label)))
    dated_tree <- ape::drop.tip(dated_tree, not_naive)
  }
    ##- UK taxa are digits only
  tip.state <- as.numeric(grepl("[A-Z]", dated_tree$tip.label))
  st <- system.time(
    clades <- find.clades.0(tree = dated_tree, tip.state = tip.state, plots = plots)
  )
  if(verbose) print(st[3])
  return(clades)
}

##---- fast_tmrca ----
#' @title Calculate TMRCA for one sequence with few same-subtype external references and other-subtype outgroup
#' @description Used in \code{\link{do_tmrca_outliers_splits}}
#' @param pathlst Path to list as .rds containing [seq = full MSA, suk = label of sequence to test, sla = labels of reference sequences, og = labels of outgroup sequences)
#' @param outdir Directory for output
#' @param plots Logical to draw plots
# #' @export
fast_tmrca <- function(pathlst = "lst.rds", outdir = ".",  plots = FALSE){
  PID <- Sys.getpid()
  lst <- readRDS(pathlst)
  seq <- lst[['seq']] # full msa
  suk <- lst[['suk']] # vector of uk labels to test
  sla <- lst[['sla']] # vector of lanl label close to outgroup
  og <- lst[['og']] # vector of outgroup labels
  mrca_results <- lapply(suk, function(x){
    msa <- seq[c(sla, x, og),]
    d <- ape::dist.dna(msa, 'TN93')
    if(any(is.na(as.matrix(d)))){
      ## misaligned outlier at distance stage
      return(c(x, rep(-1, 8)))
    }
    tr <- ape::nj(d)
    tr <- ape::root(tr, og, resolve.root = TRUE)
    if(plots) plot_tree_col_og(tr, tiplab = TRUE)
    # rc <- c(paste0('mrca(', paste(include, collapse = ','), ')'), "b(1900,2000)")
    o <- lsd_fn(tr = tr,  outgroup = og, root.constraint = NULL, verbose = FALSE)
    metrics <- o[[1]]
    dated_tree <- o[[2]]
    tip_checked <- setNames(ifelse(length(x) == 0, 0L, x), 'tip_checked')
    at_root <- setNames(dated_tree$tip.label[ape::Ntip(dated_tree)], 'at_root')
    if(plots) plot_tree_col_og(dated_tree, tmrca =  metrics['tMRCA'], tiplab = TRUE)
    return(c(tip_checked, at_root, metrics))
  })
  outfile <- paste0(outdir, '/', 'tmrca_', PID, '.rds')
  saveRDS(mrca_results, outfile)
  print(paste("output saved to", outfile))
  return(outfile)
}

##---- merge_out_tmrca ----
#' @title Parse results form fast_tmrca
#' @description Used in \code{\link{do_tmrca_outliers_splits}}. Look for Grubbs outliers of TMRCA
#' @param indir Path to individual tmrca assesments
#' @param outfile Name of outfile containing outliers of tmrca
# #' @export
merge_out_tmrca <- function(indir = "data/MSM_CRF02AG_splits", outfile = "out_tmrca"){
  ##- names
  .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", indir )
  nm <- regmatches(indir, .m)
  OUT <- paste0(indir, "/", outfile, ".rds")
  ## parse split results
  splits <- list.files(path = indir, pattern = 'tmrca', full.names = TRUE)
  mrca_results <- unlist(lapply(splits, function(x){
    a <- readRDS(x)
    return(a)
  }), recursive = FALSE)
  .ch <- unfactorDataFrame(as.data.frame(do.call(rbind, mrca_results)))
  tmrcas <- setNames(.ch$tMRCA, unname(sapply(.ch[,1], function(x) unlist(strsplit(x,'\\.'))[2])))
  #head(tmrcas)
  if(FALSE) plot(tmrcas[tmrcas > 0])
  gru <- grubbs.flag(tmrcas)
  if(FALSE) plot(gru$X [gru$X > 0], col = gru[gru$X > 0,]$Outlier + 1)
  out.tmrca <- list(rownames(gru[gru$Outlier,]))
  names(out.tmrca) <- nm
  print(dput(out.tmrca))
  print(paste("save list outliers to", OUT))
  saveRDS(out.tmrca, OUT)
}

##---- redo_lsd ----
#' Filter LSD tree for tips producing MRCA ancestral to global sequences and re date new tree
#' @param pathlsdtree Path to LSD result
#' e.g. './BS_trees/MSM_CRF02AG/
#' ExaML_result.MSM_CRF02AG.finaltree.001_lsd.rds'
#' @param outgroup Outgroup option of LSD (LSD outputs a tree without outgroup)
#' @param ... Passed to hpc.cmd.funr (i.e. hpc.q = 'pqeelab', verbose = FALSE)
#' @export
redo_lsd <- function(pathlsdtree,
                     outgroup = TRUE,
                     ...){
  print('--- start redo LSD if problematic UK MRCA ---')
  PATH_ML_TREE <- sub(pattern = "_lsd.rds", replacement =  "", x = pathlsdtree)
  PATH_NEW_ML_TREE <- sub("finaltree", "post_lsd", PATH_ML_TREE)
  PATH_NEW_LSD_TREE <- sub("finaltree", "post_lsd", pathlsdtree)

  metrics.lsd(lst = pathlsdtree)
  outliers2 <- lsd.probs(lst = pathlsdtree)
  tr <- ape::read.tree(PATH_ML_TREE)

  if( length(outliers2[[1]]) == 0 ) {
    print('no outlier, just copy LSD tree with new name')
    file.copy(from = pathlsdtree, to = PATH_NEW_LSD_TREE) ## copy same tree
  } else {
    tr2 <- ape::drop.tip(tr, outliers2[[1]])
    print('These are tips producing an unexpected MRCA')
    print(dput(outliers2[[1]]))
    ape::write.tree(tr2, file = PATH_NEW_ML_TREE)

    print('--- New LSD estimation on filtered ML tree ---')
    lsd.examl.tree(TREE = PATH_NEW_ML_TREE, outgroup = TRUE, ...)

    ##- show results
    metrics.lsd(lst = PATH_NEW_LSD_TREE)
  }
}

##---- find_clades_slice ----
#' @title Extract clades slicing from root
#' @description Recursively slice big tree to pull out clades within given size bounds
#' @param tres Phylo tree
#' @param UB Upper bound of clade size
#' @param DSLICE Height above the root at which tree is sliced
find_clades_slice <- function(tres, UB = 1e3, DSLICE = 1){
  n <- sapply( tres, function(t) length(t$tip.label))
  repeat{
    cat ( 'begin iteration \n')
    print(date())
    .tres <- tres
    tres <- list()
    for (k in 1:length(.tres)){
      if ( n[k] < UB ){
        tres[[length(tres)+1]] <- .tres[[k]]
      } else{
        cat( 'slicing tree \n')
        print( .tres[[k]] )
        print(date())
        tstres <- phytools::treeSlice(.tres[[k]], DSLICE )
        cat( 'new trees \n')
        print(tstres )
        print(date())
        for (tre in tstres ){
          tres[[length(tres)+1]] <- tre
        }
      }
    }
    n <- sapply( tres, function(t) length(t$tip.label))

    cat('end of iteration; tree sizes: \n' )
    print( sort( n ))
    cat( '#####################################\n')

    if (all( n < UB )) break
  }
  return(tres)
  # saveRDS(tres, file = OFN )
}

##---- lsd_to_clades ----
#' Apply \code{\link{find_clades_slice}} on one LSD result file
#' Remove treated patients sequences and global sequences
#' @param max.clade.size Upper bound of clade size
#' @param rm.global remove global seqs
#' @inheritParams get_clades
lsd_to_clades <- function(lsd_result, max.clade.size = 1e3, not_naive = NULL, rm.global = TRUE, verbose = TRUE){
  lsd_out <- readRDS( lsd_result )
  tr <- lsd_out[[2]]
  if(!is.null(not_naive)){
    ##- exclude not naive UK sequences
    if(verbose) print(paste("excluding", sum(tr$tip.label %in% not_naive), "experienced patients' sequences out of", length(tr$tip.label)))
    tr <- ape::drop.tip(tr, not_naive)
  }
  ##- remove global seqs
  if (rm.global) {
    glo <-tr$tip.label[grep("[A-Z].*", tr$tip.label)]
    if(verbose) print(paste("remove", length(glo), "global sequences"))
    tr <- ape::drop.tip(tr, glo)
  }
  ##- clades extraction
  clades <- find_clades_slice(tres = list(tr), UB = max.clade.size)
  return(clades)
}
