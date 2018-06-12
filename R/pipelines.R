##- set path to external data with options
##- maybe keep a readme.txt in inst/data

##---- outgroup_choice ----
#' Selection of subtype reference sequence serving as outgroup
#' @param subtype subtype of seqs
outgroup_choice <- function(subtype){
  ifelse(subtype == "B", "C", "B")
}

##---- pipeline_blastdb ----
#' @title Series of functions to generate a BLAST database
#' @description Script for:
#' \describe{
#'   \item{\code{splitForDRM}}{splitting MSA}
#'   \item{\code{bigphylo_drm}}{drm removal}
#'   \item{\code{merge_split}}{reassembling MSA}
#'   \item{\code{make_blast_db}}{make BLAST database from the resulted MSA}
#'   }
#' @param path.msa Full path name to RDS or RDA containing DNAbin object
#' @param keep.splits Keep the tmp folder containing original and drm-free splits [default = TRUE]
#' @param make.blastdb Make BLAST database of final MSA [default = FALSE]
#' @param ... Options for \code{hpc.cmd.funr}, i.e. turning on/off verbose and choosing HPC queue (\code{hpc.q})
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @export
pipeline_blastdb <- function(path.msa = "data/MSM_CRF02AG.rds",
                           keep.splits = TRUE,
                           make.blastdb = TRUE,
                           ...) # verbose = FALSE
  {
  NAME.MSA <- sub("\\..*", "", basename(path.msa)) # no extension
  DIR.MSA <- dirname(path.msa)
  print(paste("##---- start", NAME.MSA, "----##"))

  ##- create folder
  TMP <- sub(".rds|.rda|.RData", "_splits", path.msa, ignore.case = TRUE)
  dir.create(TMP, showWarnings = FALSE)
  unlink(paste0(TMP,'/*'))
  NM.OUTFILE <- paste0(TMP, '/', NAME.MSA, "_nodrm.rds")

  ##- split MSA
  files <- splitForDRM(path.seqs = path.msa, path.out = TMP)

  ##- apply rm.drm
  ## script1 on HPC contains: list of list[script, jobid, name of nodrm split] for individual split
  script1 <- hpc.cmd.array(inputs = files,
                           f.name = "bigphylo_drm",
                           f.arg = paste0("infile=$N", " ", "wait=", 1),
                           submit = TRUE,
                           ...)

  ##- apply merge_splits
  script2 <- hpc.cmd.funr(f.name = "merge_splits",
                          f.arg = paste0('indir=', TMP, ' ', 'outname=', NAME.MSA),
                          submit = TRUE,
                          jobdepend = paste0('-W depend=afterok:', script1["jobid"]),
                          ...)
  jobid2 <- get_jobid(lst = script2)

  ##- clean
  #TOFILE <- paste(DIR.MSA, basename(NM.OUTFILE), sep = '/')
  #file.copy(NM.OUTFILE, TOFILE , overwrite = TRUE)
  #if (!keep.splits) unlink(TMP)

  ##- make database
  if (make.blastdb){
    script3 <- hpc.cmd.funr(f.name = "make_blast_db",
                            f.arg = paste0('pathfile=', NM.OUTFILE),
                            other.load = "ncbi-blast/2.2.24",
                            submit = TRUE,
                            jobdepend = paste0('-W depend=afterok:', jobid2),
                            ...)
  }
}

##---- pipeline_msa ----
#' @title Series of functions to generate a MSA ready for maximum likelihood tree
#' @description Script for:
#' \describe{
#'   \item{\code{splitForDRM}}{1. splitting MSA}
#'   \item{\code{bigphylo_drm}}{2a. removing drug resistance mutations}
#'   \item{\code{blast_splits}}{2b. blasting seq splits against corresponding database}
#'   \item{\code{merge_split}}{3. reassembling MSA}
#'   \item{\code{merge_blasts}}{4. merging unique BLAST hits}
#'   \item{\code{prep_and_align_blast}}{5. aligning seq and hits}
#'   \item{\code{align_mafft}}{6. adding outgroup sequences}
#'   }
#' @param path.msa Full path name to RDS or RDA containing DNAbin object
#' @param keep.splits Keep the tmp folder containing original and drm-free splits [default = TRUE]
#' @param do.blast Blast split sequences, merge unique outputs and align [default = FALSE]
#' @param pathdb Path to BLAST database [default = "blastDB"]
#' @param do.outgroup Add outgroup [default = FALSE]
#' @param nbysplit Size of split for \code{splitForDRM}
#' @param ... Options for \code{hpc.cmd.funr}, i.e. turning on/off verbose and choosing HPC queue (\code{hpc.q})
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @export
pipeline_msa <- function(path.msa = "data/MSM_CRF02AG.rds",
                         keep.splits = TRUE,
                         do.blast = TRUE,
                         do.outgroup = TRUE,
                         pathdb = "blastDB",
                         nbysplit = 200,
                         ...) # verbose = FALSE
{
  NAME.MSA <- sub("\\..*", "", basename(path.msa)) # no extension
  DIR.MSA <- dirname(path.msa)
  print(paste("##---- start", NAME.MSA, "----##"))
  # get subtype
  .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", NAME.MSA)
  ST <- unlist(strsplit(regmatches(NAME.MSA, .m), "_"))[2] # subtype
  ##- create folder
  TMP <- sub(".rds|.rda|.RData", "_splits", path.msa, ignore.case = TRUE)
  dir.create(TMP, showWarnings = FALSE)
  unlink(paste0(TMP,'/*'))
  NM.OUTFILE <- paste0(TMP, '/', NAME.MSA, "_nodrm.rds")
  ##- find DB based on subtype
  if(do.blast){
    if(!dir.exists(pathdb)) stop("need local BLAST database folder, see pipeline_blastdb")
    DB <- find_db(infile = path.msa, dbdir = pathdb)
    ##- select max_target_seqs: = 1 for large MSM_B tree, and 5 for others
    ntarget <- ifelse(ST == "B", 1, 5)
  }
  ##- split MSA
  print(paste("##---- 1. split msa:", NAME.MSA, "----##"))
  splits <- splitForDRM(path.seqs = path.msa, path.out = TMP, n.by.split = nbysplit)
  nsplits <- length(splits)
  ##- apply rm.drm and blast on splits
  ## script1 on HPC contains: list of list[script, jobid, name of nodrm split] for individual split
  print(paste("##---- 2. remove drm and blast:", NAME.MSA, "----##"))
  script1 <- hpc.cmd.array(inputs = splits,
                           f.name = "bigphylo_drm",
                           f.arg = paste0("infile=$N", " ", "wait=", 1),
                           submit = TRUE,
                           ...)

   ##- collect names output
    nodrm <- sapply(splits, function(x) sub(".rds|.rda", "_nodrm.rds", x, ignore.case = TRUE))
    blastout <- sapply(splits, function(x) sub(".rds|.rda", "_nodrm_blast.out", x, ignore.case = TRUE))

    ##- blast
    script1b <- hpc.cmd.array(inputs = nodrm,
                             f.name = "blast_splits",
                             f.arg = paste0("infile=$N", " ", "database=", DB, " ", "target_seqs=", ntarget),
                             other.load = "ncbi-blast/2.2.24",
                             submit = TRUE,
                             jobdepend = paste0('-W depend=afterok:', script1["jobid"]),
                             ...)

  ##- apply merge_splits
  print(paste("##---- 3. merge nodrm splits:", NAME.MSA, "----##"))
  script2 <- hpc.cmd.funr(f.name = "merge_splits",
                          f.arg = paste0("nsplits=", nsplits, " ", "indir=", TMP, " ", "outname=", NAME.MSA),
                          submit = TRUE,
                          jobdepend = paste0('-W depend=afterok:', script1["jobid"]),
                          ...)
  if (do.blast){
    ##- merge blast
    print(paste("##---- 4. merge blast splits:", NAME.MSA, "----##"))
    script5 <- hpc.cmd.funr(f.name = "merge_blasts",
                            f.arg = paste0("nblasts=", nsplits, " ", 'indir=', TMP, ' ', "outfile=", NM.OUTFILE),
                            submit = TRUE,
                            jobdepend = paste0('-W depend=afterok:', script1b["jobid"]),
                            ...)

    ##- do mafft
    print(paste("##---- 5. align seqs and blast hits:", NAME.MSA, "----##"))
    script6 <- hpc.cmd.funr(f.name = "prep_and_align_blast",
                            f.arg = paste0("infile=", NM.OUTFILE),
                            other.load = "mafft/7.271",
                            hpc.resource = LARGEPBS,
                            submit = TRUE,
                            jobdepend = paste0('-W depend=afterok:', script5["jobid"]),
                            ...)
  }
  if(do.outgroup){
    ##- add outgroup
    data("ref_hiv_msa", package = "tenbrit", envir = environment())
    OGST <- outgroup_choice(ST)
    OG <- paste0("data/", "outgroup_", OGST, ".fas")
    if(!file.exists(OG)) ape::write.dna(ref_hiv_msa[[OGST]], file = OG, format = 'fasta')
    UK <- sub("nodrm.rds", "nodrm_ref.fas", NM.OUTFILE)
    UKOG <- sub(".fas", "_og.fas", UK)
    print(paste("##---- 6. add", OGST, "outgroup:", NAME.MSA, "----##"))
    script7 <- hpc.cmd.funr(f.name = "align_mafft",
                            f.arg = paste0("fas1=", OG, " ", "fas2=", UK, " ", "out=", UKOG),
                            other.load = "mafft",
                            hpc.resource = LARGEPBS,
                            submit = TRUE,
                            jobdepend = paste0('-W depend=afterok:', script6["jobid"]),
                            ...)
    ##- TODO: copy everything back, use $TMPDIR
  }
  return(NULL)
}

##---- do_examl ----
#' @title Compute ExaML tree
#' @description Either in HPC or locally
#' @param FAS Path to fasta MSA
#' @param DIR Folder for input MSA and output tree
#' @param bsdebut Boostrap debut
#' @param bsend Bootstrap end
#' @param q HPC queue
#' @param mem memory
#' @return If on HPC, output PBS_JOBID number(s)
#' @export
do_examl <- function(FAS = "data/MSM_CRF02AG_splits/MSM_CRF02AG_nodrm_ref_og.fas",
                     DIR = "data/ML_TREE0",
                     bsdebut = 0,
                     bsend = 2,
                     mem = "32GB",
                     q = NA # 'pqeelab'
                     ){
  jobid <- NULL
  dir.create(DIR, showWarnings= FALSE)
  FULLDIR <- normalizePath(DIR) # full name
  ##- HPC
  if( nchar( Sys.which("qsub")) ){
  ## copy as *.R
  INFILE <- sub('.fas', '.R', basename(FAS))
  seq <- ape::read.dna(FAS, format = 'fasta')
  save(seq, file = paste(DIR, INFILE, sep = '/'))
  ##- run
  wd <- setwd(DIR) ## return current
  jobid <- big.phylo::pipeline.ExaML.bootstrap.per.proc(indir=FULLDIR, infile=INFILE, bs.from=bsdebut, bs.to=bsend, bs.n=bsend, hpc.walltime=24, hpc.q=q, hpc.mem=mem, hpc.nproc=16)
  on.exit(setwd(wd)) ## back home no matter what
  } else {
  tofile = paste(DIR, sub('\\.[a-z]{1,}$', '.fasta', basename(FAS), perl = TRUE), sep = '/')
  file.copy(from = FAS, to = tofile)
  INFILE <- basename(tofile)
  ##- local machine (fasta input)
  cmd <- big.phylo::cmd.examl.bootstrap.on.one.machine(indir=FULLDIR, infile=INFILE, bs.from=bsdebut, bs.to=bsend)
  }
  return(jobid)
}

##---- do_rtt ----
#' RTT filtering
#' @inheritParams rm_rtt_outliers
#' @param ... Passed to hpc.cmd.funr (i.e. hpc.q = 'pqeelab', verbose = FALSE)
#' @export
do_rtt <- function(TREE = "data/ML_TREE0/ExaML_result.MSM_CRF02AG_nodrm_ref_G.finaltree.000", ...)
{
  hpc.cmd.funr(f.name = "rm_rtt_outliers",
               f.arg = paste0("TREE=", TREE),
               submit = TRUE,
               ...)
}

##---- do_lsd ----
#' LSD on HPC
#' @inheritParams lsd.examl.tree
#' @param ... Passed to hpc.cmd.funr (i.e. hpc.q = 'pqeelab', verbose = FALSE)
#' @export
do_lsd <- function(TREE = "data/ML_TREE0/ExaML_result.MSM_CRF02AG_nodrm_ref_G.finaltree.000", outgroup = TRUE,
                           ...){
hpc.cmd.funr(f.name = "lsd.examl.tree",
             f.arg = paste0("TREE=", TREE, " ", "outgroup=", outgroup),
             #hpc.resource = SMALLPBS,
             submit = TRUE,
             #jobdepend = "",
             ...)
}

##---- do_tmrca_outliers ----
#' @title Finding tmrca outliers on HPC
#' @description See find_tmrca_outliers
#' @param PATHMSA Path of msa file
#' @param ... Passed to hpc.cmd.funr (i.e. hpc.q = 'pqeelab', verbose = FALSE)
#' @export
do_tmrca_outliers <- function(PATHMSA = 'data/ML_TREE0/MSM_CRF02AG_nodrm_ref_og.R',
                   ...){
  hpc.cmd.funr(f.name = "find_tmrca_outliers",
               f.arg = paste0("pathmsa=", PATHMSA),
               other.load = "ncbi-blast/2.2.24",
               hpc.resource = SMALLPBS,
               hpc.time = "#PBS -l walltime=24:00:00",
               submit = TRUE,
               #jobdepend = "",
               ...)
}


##---- do_tmrca_outliers_splits ----
#' @title Finding tmrca outliers on HPC with splitting step
#' @description See find_tmrca_outliers
#' @param pathmsa Path of msa file
#' e.g. "data/MSM_CRF02AG_splits/MSM_CRF02AG_nodrm_ref_og.fas"
#' @param ... Passed to hpc.cmd.funr (i.e. hpc.q = 'pqeelab', verbose = FALSE)
#' @export
do_tmrca_outliers_splits <- function(
  pathmsa,
  ...){
  load_msa(pathmsa) ## as seq
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
  ##- splits
  n.by.split <- 100
  sp <- split(new_uk_lbl, ceiling(seq_along(new_uk_lbl) / n.by.split  ))
  inputs <- lapply(1:length(sp), function(i){
    lst <- list(seq = seq, suk = sp[[i]], sla = include, og = og)
    pathlst <- paste0(dirname(pathmsa), '/', 'lst_',i, '.rds')
    saveRDS(lst, pathlst)
    return(pathlst)
  })
  print(inputs)
  print('clean previous results')
  to_rm <- list.files(path = dirname(pathmsa), pattern = 'tmrca', full.names = TRUE)
  file.remove(to_rm)
  ##- run
  script1 <- lapply(inputs, function(x){
    hpc.cmd.funr(f.name = 'fast_tmrca',
                 f.arg = paste0("pathlst=", x, " ", "outdir=", dirname(pathmsa)),
                 submit = TRUE,
                 ...)
  })
  ##- merge
  jobid1 <- get_jobid(lst = script1)
  script2 <- hpc.cmd.funr(f.name = "merge_out_tmrca",
                          f.arg = paste0("indir=", dirname(pathmsa)),
                          submit = TRUE,
                          jobdepend = paste0('-W depend=afterok:', jobid1),
                          ...)
  return(inputs)
  ## TODO: clean lst and tmrca
  }

##---- select_cohort ----
#' Make subset MSA for samples collected as of year x
#'
#' From year contained in \code{dateres} in \code{df}
#' @param pathmsa path to RDS file of MSA
#' @param year first year of selection
#' @export
select_cohort <- function(pathmsa, year = 2010){
  NM0 <- sub(".rds","", pathmsa)
  NM <- paste0(NM0,"_",year,".rds")
  data(df, envir = environment())
  yres <- as.numeric(substr(df$dateres, 1,4))
  subid <- df$testindex[yres >= year]
  a <- readRDS(pathmsa)
  b <- a[names(a) %in% subid]
  print(paste("save subset MSA of length", length(b), "in", NM))
  saveRDS(b, file = NM)
  return(length(b))
}



