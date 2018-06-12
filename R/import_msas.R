##-- Functions to process data, starting from data-raw
GROUP_NAME_PTRN <- ".*((MSM|HET)_(CRF02AG|A|B|C)).*"

##---- stratify_msa ----
#' Extract MSAs stratified by main group/subtype
#'
#' @return Save MSAs as RDS in \code{data}
#' @export

stratify_msa <- function(){
  dir.create("data", showWarnings = FALSE)
  print("copy MSAs stratified by main group/subtype in './data/'")
  for(i in levels(tenbrit::df$tran2)[-length(levels(tenbrit::df$tran2))] ){
    for(j in levels(tenbrit::df$subtype)[-length(levels(tenbrit::df$subtype))] ){
      selec <- tenbrit::df[tenbrit::df$tran2 == i & tenbrit::df$subtype == j, 'testindex']
      seqs <- tenbrit::s950[ names(tenbrit::s950) %in% selec]
      #browser()
      print(paste(length(selec), 'sequences for', i, j, 'as', class(seqs)))
      saveRDS(seqs, file = paste0('./data/', i, '_', j, '.rds') )
    }
  }
}

#' Load LANL sequences in data for processing
#'
#' Save RDS instead of rda in package installation ?
#' @export
lanl_msa <- function(){
  dir.create("data", showWarnings = FALSE)
  d <- data(package = 'tenbrit')
  nm <- d$results[, "Item"]
  lst <- nm[grepl("^db_.*$", nm)]
  data(list = lst, package = 'tenbrit', envir = environment() )
  lapply(lst, function(x) {
    print(paste("Save", paste0("data/",x,".rds")))
    saveRDS(get(x), file = paste0("data/",x,".rds"))
    })
}





