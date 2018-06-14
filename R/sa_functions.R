#' Compute infector probabilities on clades
#' @param tr phylo tree
#' @param parms_SA Parameters from \code{\link{get_sa_parms_range}}
#' @param p proportion of subtype represented (useless)
#' @return dataframe(donor, recip, ip)
#' @details \code{phydynR::phylo.source.attribution.hiv.msm} wraps \code{phydynR::phylo.source.attribution.multiDeme.fgy} that runs \code{phydynR::sourceAttribMultiDemeCpp} (mode = 1) or \code{phydynR::sourceAttribMultiDemeCpp2} (mode = 2)
#' @importFrom ape reorder.phylo
#' @export
    sa_by_clades <- function(tr, parms_SA =  get_sa_parms_range(), p = 1){
      require(ape) ## temp fix for ape::reorder.phylo in phylo.source.attribution.multiDeme.fgy ## importFrom is not enough ?
      CD4s <- parms_SA[['CD4s']][tr$tip.label]
      STs <- parms_SA[['STs']][tr$tip.label]
      EHIs <- parms_SA[['EHIs']][tr$tip.label]
      MH <- parms_SA[['MH']]
      # apply prop of subtype p
      PLWHIV <- parms_SA[['PLWHIV']] * p
      NEWINF <- parms_SA[['NEWINF']] * p

      ##- run SA
      W <- phydynR::phylo.source.attribution.hiv.msm(
        tree = tr,
        sampleTimes = STs, # years
        cd4s = CD4s,
        ehi = EHIs,
        numberPeopleLivingWithHIV = PLWHIV,
        numberNewInfectionsPerYear = NEWINF,
        maxHeight = MH,
        res = 1e3,
        treeErrorTol = Inf,
        minEdgeLength = 1/52,
        mode = 2
      )

      ##- return a DF
      return(data.frame(donor = W[[1]],
                          recip = W[[2]],
                          ip = as.numeric(W[[3]]),
                           stringsAsFactors = FALSE))
    }

#' Compute infector probabilities and bind all clade results from tree
#' @param pathclade Path to list of clades
#' @param parms_SA Parameters from \code{\link{get_sa_parms_range}}
#' @param outdir Create directory in which to save results
#' @param verbose print
#' @return List of [(inc, prev), dataframe(donor, recip, ip)]
#' @export
 get_ip <- function(pathclade, parms_SA =  get_sa_parms_range(), outdir = "data/IPS", verbose = TRUE){
   ##- names
   PID <- Sys.getpid()
   dir.create(outdir, showWarnings = FALSE)
   bn <- sub("\\.rds", paste0("_", PID, "\\.rds"), sub("clades_", "ip_", basename(pathclade)) )
   IP <- paste(outdir, bn, sep = "/")
   ##-
   if(!file.exists(IP)){
     subtype <- get_subtype(pathclade)
     # prop of subtype
     p <- parms_SA[['p_msm']][subtype]

     clades <- readRDS(pathclade)
     w.clades <- lapply(clades, function(x){
       sa_by_clades(tr = x,
                    parms_SA = parms_SA,
                    p = p)
     })
     w.cl <- do.call(rbind, w.clades)
     lst <- list(parms = parms_SA[c("MH","NEWINF", "PLWHIV")], W = w.cl)
     ##- save
     if(verbose) print(paste('Save', IP))
     saveRDS(lst, file = IP)
   } else {
     if(verbose) print(paste('File exists:', IP))
   }
      return(IP)
 }


 ##---- get_sa_parms_range ----
 #' Default parameters for source attribution in MSM.
 #'
 #' Surveillance data taken from Brown, 2017 and Yin, 2014.
 #' \itemize{
 #' \item For incidence in MSM:
 #' "Following this adjustment, the estimated number of infections acquired per year rose from around 2,200 infections (95\% credible interval (CrI)1,800 to 2,500) in 2007 to a peak of 2,800 (CrI 2,300 to 3,200) in 2012 before falling to 1,700 (CrI 900 to 2,700) in 2016"
 #' \item Let's assume in the period preceding last sequence (2013) an incidence of 2800 [2300 - 3200].
 #' \item For prevalence, from Yin et al. 2014, MSM living with HIV in 2013 = 43,500 [40,200 - 48,200]
 #' \item Let's further assume that incidence and prevalence are normally distributed and credible interval bounds are equal to 95\% confidence interval bounds, so that \code{sd = (up - lo)/(2 * qnorm(.975))}
 #' \item Hence, values of incidence / prevalence can be drawn from normal distributions
 #' }
 #' @param filename Whether and where to save file
 #' @references
 #' 1. Brown AE, Kirwan P, Chau C, Khawam J, Gill ON, Delpech VC. Towards elimination of HIV transmission AIDS and HIV related deaths in the UK - 2017 report. Public Health England; 2017.
 #' \url{https://www.gov.uk/government/uploads/system/uploads/attachment_data/file/662306/Towards_elimination_of_HIV_transmission_AIDS_and_HIV_related_deaths_in_the_UK.pdf}
 #' 2. Yin Z, Brown AE, Hughes G, Nardone A, Gill ON, Delpech VC, et al. HIV  in  the  United  Kingdom  2014  Report:  data  to  end  2013. London: Public Health England; 2014.
 #' \url{https://www.gov.uk/government/uploads/system/uploads/attachment_data/file/401662/2014_PHE_HIV_annual_report_draft_Final_07-01-2015.pdf}

 #' @export
 get_sa_parms_range <- function(filename = NA){

     data(list = "df", package = 'tenbrit', envir = environment())
     ## for MSM
     sa_df <- df[df$tran2 == "MSM",]
     ##- Maximum height
     MH <- 20
     ##- incidence, prevalence: draw from norm
     prev <- setNames(c(40200, 43500, 48200), c('lo', 'mu', 'up'))
     inc <- setNames(c(2300, 2800, 3200), c('lo', 'mu', 'up'))
     draw_norm_ci <- function(n, var, plots = FALSE){
       lo <- sort(var)[1]
       mu <- sort(var)[2]
       up <- sort(var)[3]
       sd <- (up - lo)/(2 * stats::qnorm(0.975))
       x <- rnorm(n = n, mean = mu, sd = sd)
       if(plots){
         graphics::hist(x, breaks = 30)
         graphics::abline(v = mu, col = 'red')
         graphics::abline(v = lo, col = 'blue')
         graphics::abline(v = up, col = 'blue')
       }
       return(x)
     }
     NEWINF <- draw_norm_ci(n = 1, var = inc)
     PLWHIV <- draw_norm_ci(n = 1, var = prev)
     CD4s <- setNames(sa_df$cd4_365, sa_df$testindex)
     EHIs <- setNames(sa_df$stage == 1, sa_df$testindex)
     STs <- setNames(as.numeric(sa_df$dateres + 14)/365.25 + 1970, sa_df$testindex)
     # proportion subtype in MSM
     tt <- table(df$tran2, df$subtype, useNA='ifany')
     p_msm <- as.matrix(prop.table(tt,1))[1,]
     # head(CD4s); head(EHIs); table(EHIs); head(STs)
     parms_SA <- list(MH = MH, NEWINF = NEWINF,
                      PLWHIV = PLWHIV, STs = STs,
                      CD4s = CD4s, EHIs = EHIs,
                      p_msm = p_msm)
     if(!is.na(filename)){
     ##- save
     print(paste("save", filename))
     saveRDS(parms_SA, file = filename)
   }
  return(parms_SA)
 }
