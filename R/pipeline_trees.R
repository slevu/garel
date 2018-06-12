LSTFAS <- list.files(pattern = "_nodrm_ref_og.fas$", recursive = TRUE, full.names=TRUE)
(pathmsa <- LSTFAS[4])

##---- pipeline.tree ----
#' Pipeline for tree0 and filtering
#' @param pathmsa path to fasta
#' @param dirtree dir of output
#' @param Q HPC queue
#' @export
  pipeline_tree0 <- function(
    pathmsa = "./data/MSM_CRF02AG_splits/MSM_CRF02AG_nodrm_ref_og.fas",
    dirtree = "data/MLTREE0",
    Q = NA){ # "pqeelab"
    ##- names
    .m <- regexpr("(MSM|HET)_(CRF02AG|A|B|C)", pathmsa )
    NM <- regmatches(pathmsa, .m)
    ##- shorten name copy
    dir.create(dirtree, showWarnings= FALSE)
    MSA <- paste0(dirtree, "/", NM, ".fas")
    TREE0 <- paste(dirtree, paste0("ExaML_result.", NM, ".finaltree.000"), sep = "/")
    FLTDTREE0 <- sub("\\.finaltree", "_fltd\\.finaltree", TREE0)
    LSDTREE0 <- paste0(FLTDTREE0, "_lsd.rds")
    FLTDOUT <- sub(".fas", "_fltd.outliers.rds", MSA)
    FLTDMSA <- sub(".fas", "_fltd.fas", MSA)
    file.copy(pathmsa, MSA)

    print("##---- 1. first tree ----")
    # system('rm data/MLTREE0/*')
    jobid1 <- do_examl(FAS = MSA, DIR = dirtree, bsdebut=0, bsend=0, q = Q)
    jobid <- unlist(jobid1)
    njobid <- as.numeric(regmatches(jobid, regexpr("[[:digit:]]+", jobid)))

    print("##---- 2. filter outliers ----")
    script2 <- hpc.cmd.funr(f.name = "filter_tree_outliers",
                 f.arg = paste0("TREE=", TREE0),
                 submit = TRUE,
                 jobdepend = paste0('-W depend=afterok:', njobid),
                 hpc.q = Q)

    jobid2 <- get_jobid(lst = script2, qstat = FALSE)

    print("##---- 3. remove outliers from MSA ----")
    script3 <- hpc.cmd.funr(f.name = "rm_outliers_msa",
                            f.arg = paste0("pathmsa=", MSA, " ", "pathoutliers=", FLTDOUT),
                            submit = TRUE,
                            jobdepend = paste0('-W depend=afterok:', jobid2),
                            hpc.q = Q)

    print("##---- 4. first LSD ----")
    script4 <- do_lsd(TREE = FLTDTREE0,
                      outgroup = TRUE,
                      hpc.q = Q,
                      jobdepend = paste0('-W depend=afterok:', jobid2))
    return(FLTDMSA)
  }

##---- pipeline_bstrees  ----
#  FLTDMSAS <- list.files(pattern = "_fltd.fas$", recursive = TRUE, full.names=TRUE)
#  (pathmsa <- FLTDMSAS[4])
#' Pipeline for bootstrap trees
#' @inheritParams pipeline_tree0
#' @param dirbstree output dir
#' @param bs_only turn off LSD steps
#' @param BSD debut bootstrap
#' @param BSE end bootstrap
#' @export
pipeline_bstrees <- function(
  pathmsa = "data/MLTREE0/MSM_CRF02AG_fltd.fas",
  dirbstree = "data/BSTREES",
  bs_only = FALSE,
  BSD=1,
  BSE=2,
  Q = NA){ # "pqeelab"
  jobid1 <- jobid2  <- NULL # to trigger job dependency
  ## names
  n <- sprintf("%03d", BSD:BSE) # "001" "002"
  NM <- sub(".fas", "", basename(pathmsa))
  BSTREES <- paste(dirbstree, paste0("ExaML_result.", NM, ".finaltree.", n), sep = "/")
  BSLSD <- paste0(BSTREES, "_lsd.rds")
  POSTLSD <- sub("finaltree", "post_lsd", BSLSD)

  print("##---- 5. Bootstrap trees ----")
  if( all(file.exists(BSTREES)) ){
    print('All BS trees already exist, do nothing')
  } else {
    jobid1 <- do_examl(FAS = pathmsa, DIR = dirbstree, bsdebut = BSD, bsend = BSE, q = Q)
    .j <- unlist(jobid1)
    njobid <- as.numeric(regmatches(.j, regexpr("[[:digit:]]+", .j)))
    dep1 <- cbind(njobid, BSTREES)
    print(dep1)
  }

  if(!bs_only){
    print("##---- 6. First LSD ----")
    if( all(file.exists(BSLSD)) ){
      print('All LSD trees already exist, do nothing')
    } else {
      script2 <- lapply(BSTREES, function(x){
        Sys.sleep(.2)
        do_lsd(TREE = x,
               outgroup = TRUE,
               jobdepend = if(is.null(jobid1)){NULL} else {paste0('-W depend=afterok:', dep1[BSTREES == x, "njobid"]) },
               hpc.q = Q)
      })

      jobid2 <- get_jobid(lst = script2, array = FALSE, qstat = FALSE, verbose = FALSE)
      dep2 <- cbind(jobid2, BSLSD)
      print(dep2)
    }


    print("##---- 7. Second LSD ----")
    if( all(file.exists(POSTLSD)) ){
      print('All POST LSD trees already exist, do nothing')
    } else {
    script3 <- lapply(BSLSD, function(x){
      Sys.sleep(.2)
      hpc.cmd.funr(f.name = "redo_lsd",
                   f.arg = paste0("pathlsdtree=", x, " ", "outgroup=", TRUE),
                   jobdepend = if(is.null(jobid2)) {NULL} else {paste0('-W depend=afterok:', dep2[BSLSD == x, "jobid2"]) },
                   submit = TRUE,
                   hpc.q = Q)
    })
  }
  }
  # plot_all_lsd_tree(POSTLSD)
}


##---- extract_clades ----
#' Extract clades from list of LSD bootstrap trees
#' @param lsd_results List of LSD results pathname (\code{*.rds})
#' @param dirclades Where to save clades
#' @param max.clade.size For \code{\link{find_clades_slice}}, maximum clade size
#' @param verbose Details
#' @return List of clades filenames
#' @export
extract_clades <- function(lsd_results, dirclades = "data/CLADES", max.clade.size = 1e3, verbose = TRUE){
  dir.create(dirclades, showWarnings= FALSE)
  not_naive <- get_notnaive()
  d <- lapply(lsd_results, function(x){
    ## names
    NM <- sub(GROUP_NAME_PTRN, "\\1", x) # name of group
    N <- sub(".*(\\d{3}).*", "\\1", x) # number bootstrap
    CLADE <- paste0(dirclades, "/", "clades_", NM, "_", N, ".rds")
    #list.clade.trees.naive <- get_clades(lsd_result = x, not_naive = not_naive, verbose = verbose)
    list.clade.trees.naive <- lsd_to_clades(lsd_result = x, max.clade.size = max.clade.size, not_naive = not_naive, rm.global = TRUE, verbose = verbose)
    ## save
    if(verbose) print(paste("saving to", CLADE))
    saveRDS(list.clade.trees.naive, file = CLADE)
    return(CLADE)
  })
  return(unlist(d))
}
