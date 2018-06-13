PKG <- "garel"
RLOAD <- "anaconda3/personal" #"R/3.3.3"
OUTSH <- "code"
TIMEPBS <- "#PBS -l walltime=24:00:00"
SMALLPBS <- "#PBS -l select=1:ncpus=1:mem=1gb"
LARGEPBS <- "#PBS -l select=1:ncpus=16:mem=8gb"
#PATHRS <- system.file("ext", "do_funr.Rscript", package = PKG, mustWork = TRUE)

####---- hpc.cmd.funr ----
#' @title Create cmd with do_funr.Rscript
#' @description Execute any R function as Rscript, either by wrapping for HPC or as a bash script.
#' @param f.name Function name in quotes [default = norm_test], example function
#' @param f.arg Function arguments in the form \code{"arg1=x arg2=y"}
#' @param r.load R version to load with \code{module load}
#' @param other.load Other software to load with \code{module load}
#' @param hpc.time Job resource selection: time [default = TIMEPBS]
#' @param hpc.resource Job resource selection: core / node / memory [default = SMALLPBS]
#' @param hpc.q PBS queue name [default = NA]
#' @param submit Run either with \code{qsub} or \code{bash} [default = FALSE]
#' @param jobdepend In the form \code{depend=afterok:previous_jobs}
#' @param verbose Write stuff
#' @return List containing script to execute with \code{qsub} or \code{bash}. If \code{submit = TRUE}, run the script. If \code{submit = TRUE} and on HPC, output the PBS_JOBID number
#' @examples
#' hpc.cmd.funr("norm_test", "n=3", submit = FALSE)
#' @export
hpc.cmd.funr <- function(f.name = "norm_test",
                         f.arg = "n=5",
                         r.load = RLOAD,
                         other.load = "",
                         hpc.time = TIMEPBS,
                         hpc.resource = SMALLPBS,
                         hpc.q = NA,
                         submit = FALSE,
                         jobdepend = NULL,
                         verbose = TRUE)
{
  #browser()
  PATHRS <- system.file("ext", "do_funr.Rscript", package = PKG, mustWork = TRUE)
  cmd <- paste("Rscript", PATHRS, f.name, f.arg)
  if (verbose) {print(cmd)}
  bang <- "#!/usr/bin/env bash"

  ##- HPC
  if( nchar( Sys.which("qsub") ) ){
    ##- dir for log files
    dirlog <- paste0("stdout", "/", f.name, ".log")
    dir.create("stdout", showWarnings = FALSE)
    dir.create(dirlog, showWarnings = FALSE)
    ##- header
    if (!is.na(hpc.q)) q <- paste("#PBS -q", hpc.q) else q <- ""
    hpc.head <- paste(bang, hpc.time, hpc.resource, paste("#PBS -o", dirlog, "-j oe"), q, "cd $PBS_O_WORKDIR", paste("module load", r.load, other.load), sep ="\n ")
    fullscript <- paste(hpc.head, cmd, sep = "\n")
    path.qsub <- tempfile(pattern = paste0(f.name,"_"), fileext = ".sub") # tmpdir = "."
    writeLines(fullscript, con = path.qsub)
    if(verbose) print(paste("Writing qsub script in", path.qsub))
    args <- ifelse(!is.null(jobdepend), paste(jobdepend, path.qsub), path.qsub)
    if(submit){
      jobid <- system2("qsub", args, stdout = TRUE)
      #extract first numeric part of jobid (e.i. "543220.cx1")
      njobid <- as.numeric(regmatches(jobid, regexpr("[[:digit:]]+", jobid)))
      return(list(script = fullscript, jobid = njobid))
    } else {
      return(list(script = fullscript))
    }
  }
  ##- local machine
  else {
    fullscript <- paste(bang, cmd, sep = "\n")
    path.sh <- tempfile(pattern = paste0(f.name,"_"), fileext = ".sh") # tmpdir = "."
    writeLines(fullscript, con = path.sh)
    if(verbose) print(paste("Writing bash script in", path.sh))
    if(submit){
      system(paste("bash", path.sh))
    }
    return(list(script = fullscript))
  }
}

####---- hpc.cmd.array ----
#' Create HPC command with qsub array jobs
#'
#' See \code{\link{hpc.cmd.funr}}.
#' If not on HPC, bash scripts are generated using lapply().
#' @param inputs R list of inputs or infiles, passed to an bash array and indexed as \code{$PBS_ARRAY_INDEX}
#' @param f.name Function name in quotes [default = norm_test], example function
#' @param f.arg Function arguments in the form \code{"arg1=x arg2=y"}
#' @param r.load R version to load with \code{module load}
#' @param other.load Other software to load with \code{module load}
#' @param hpc.time Job resource selection: time [default = TIMEPBS]
#' @param hpc.resource Job resource selection: core / node / memory [default = SMALLPBS]
#' @param hpc.q PBS queue name [default = NA]
#' @param submit Run either with \code{qsub} or \code{bash} [default = FALSE]
#' @param jobdepend In the form \code{depend=afterok:previous_jobs}
#' @param verbose Write stuff
#' @return List containing script to execute with \code{qsub} or \code{bash}. If \code{submit = TRUE}, run the script. If \code{submit = TRUE} and on HPC, output the PBS_JOBID number
#' @examples
#' hpc.cmd.array(inputs = list(3,5), f.name = "norm_test", f.arg = "n=$N", submit = TRUE)
#' @seealso \code{\link{hpc.cmd.funr}}
#' @export
hpc.cmd.array <- function(inputs = list(3,5),
                          f.name = "norm_test",
                          f.arg = "n=$N",
                          r.load = RLOAD,
                          other.load = "",
                          hpc.time = TIMEPBS,
                          hpc.resource = SMALLPBS,
                          hpc.q = NA,
                          submit = FALSE,
                          jobdepend = NULL,
                          verbose = TRUE)
{
  #browser()
  PATHRS <- system.file("ext", "do_funr.Rscript", package = PKG, mustWork = TRUE)
  bang <- "#!/usr/bin/env bash"

  ##- HPC
  if( nchar( Sys.which("qsub") ) ){
  ins <- paste(paste0("IN=(",paste(inputs, collapse = ' '),")"), "N=${IN[$PBS_ARRAY_INDEX - 1]}", sep ="\n ") # zero indexed array
  cmd <- paste(ins, paste("Rscript", PATHRS, f.name, f.arg), sep ="\n ")
  if (verbose) {print(cmd)}
      ##- dir for log files
    dirlog <- paste0("stdout", "/", f.name, ".log")
    dir.create("stdout", showWarnings = FALSE)
    dir.create(dirlog, showWarnings = FALSE)
    ##- header
    if (!is.na(hpc.q)) q <- paste("#PBS -q", hpc.q) else q <- ""
    hpc.head <- paste(bang, hpc.time, hpc.resource, paste("#PBS -o", dirlog, "-j oe"), q, "cd $PBS_O_WORKDIR", paste("module load", r.load, other.load), sep ="\n ")
    fullscript <- paste(hpc.head, cmd, sep = "\n")
    path.qsub <- tempfile(pattern = paste0(f.name,"_"), fileext = ".sub") # tmpdir = "."
    writeLines(fullscript, con = path.qsub)
    if(verbose) print(paste("Writing qsub script in", path.qsub))
    array <- paste0("-J ", "1-", length(inputs))
    args <- ifelse(!is.null(jobdepend), paste(jobdepend, array, path.qsub), paste(array, path.qsub))
    if(submit){
      jobid <- system2("qsub", args, stdout = TRUE) # system("qstat -t")
      #extract first numeric part of jobid (e.i. "543220.cx1")
      njobid <- as.numeric(regmatches(jobid, regexpr("[[:digit:]]+", jobid)))
      return(list(script = fullscript, jobid = paste0(njobid, "[]")))
    } else {
      return(list(script = fullscript))
    }
  }
  ##- local machine
  else {
    lscripts <- lapply(inputs, function(x){
      ins <- paste0("N=", x)
      cmd <- paste(ins, paste("Rscript", PATHRS, f.name, f.arg), sep ="\n ")
      fullscript <- paste(bang, cmd, sep = "\n")
      path.sh <- tempfile(pattern = paste0(f.name,"_"), fileext = ".sh") # tmpdir = "."
      writeLines(fullscript, con = path.sh)
      if(verbose) print(paste("Writing bash script in", path.sh))
      if(submit){
        system(paste("bash", path.sh))
      }
      return(fullscript)
    })
    return(lscripts)
  }
}

##---- get_jobid ----
#' Get PBS_JOBID from qsub query
#'
#' Check depth of list to find \code{jobid} object
#' @param lst List (or list of list) of qsub
#' @param array array jobs notation
#' @param qstat print qstat
#' @param verbose print stuff
#' @return Numeric of JOBID in form: single integer, array notation x:y, or vector if array = FALSE
get_jobid <- function(lst = "script", array = TRUE, qstat = FALSE, verbose = TRUE){
  if (nchar(Sys.which("qsub")) > 0) { # if HPC, record JOBID
    if (qstat) system('qstat')
    if (list.depth(lst) == 1) {
      jobs <- lst[['jobid']]
    } else {
      jobs <- sapply(lst, function(x) x[["jobid"]])
      if (verbose) print(paste(length(jobs), "jobs launched"))
      if (array) jobs <- paste0(min(jobs), ":", max(jobs))
    }
    if (verbose) print(paste0("previous_jobs=", jobs))
    return(jobs)
  } else {
    return(NA)
  }
}

