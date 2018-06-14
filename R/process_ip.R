##---- mean_ip_per_bs ----
#' Mean infector probabilities by bootstrap
#' @param pathip List by group of infector probability files
#' @export
mean_ip_per_bs <- function(pathip, outfile = NULL){
  ## regroup path by subtype then by BS
  IPS <- name.list.group(pathip, regroup = TRUE)
  dd <- lapply(IPS, function(g){
    m <- as.numeric(sapply(g, function(x) unlist(strsplit(basename(x), "_"))[4]) ) ## 001, 002, etc.
    lapply(unique(m), function(i){
      g[i == m]
    })
  })
  ##- read all in a list: big !: list_ip_bs[[subtype]][[BS]][[list(n, bind df)]]
  system.time(
    list_ip_bs <- lapply(dd, function(g){
      lapply(g, function(bs){
        l <- lapply(bs, function(x){
          ip <- readRDS(x)
          ip[['W']]
        })
        list(n = length(l), DT = data.table::rbindlist(l))
      })
    })
  )

  ##- mean ip by donor-recip pair
  system.time(
    list_mean_ips_per_bs <- lapply(list_ip_bs, function(g){
      lapply(g, function(x){
        n <- x[['n']]
        DT <- x[['DT']]
        DTa <- DT[, sum(ip)/n, by=list(donor, recip)][order(donor, recip)]
        names(DTa)[names(DTa)=='V1'] <- 'ip'
        return(as.data.frame(DTa))
      })

    })
    # just new mean variable # do not use mean(ip) but rather consider not existing pairs as zero infector prob
  )
  ##- save
  if(!is.null(outfile)){
    print(paste("save", outfile))
    saveRDS(list_mean_ips_per_bs, file = outfile)
  }
  list_mean_ips_per_bs
}
#mean_ip_per_bs(pathip = IP)

##---- ip_by_var ----
#' Edge list of infector probabilities by patients variable
#' @param iplist list of ip by patient
#' @param df patient data
#' @param var variable
#' @export
ip_by_var <- function(iplist, df, var){
  #df <- as.data.frame(df) # prevent data.table error
  donor_var <- df[match(iplist[, 'donor'], df[,'testindex']), var]
  recip_var <- df[match(iplist[, 'recip'], df[,'testindex']), var]
  new_df <- data.frame(donor_var, recip_var, as.numeric(iplist$ip), stringsAsFactors = FALSE)
  colnames(new_df) <- c(paste0('donor', '_', substitute(var)),
                        paste0('recip', '_', substitute(var)),
                        'ip')
  return(new_df)
}
