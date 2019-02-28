
##---- high_contrast ----
high_contrast <- function(x) {
  (abs(x)^(1/3)) * sign(x)
}

##---- plot_w_reg ----
#' XY bubble plot of weighted regression
#' @param df dataframe with x,y, weights
#' @param x x
#' @param y y
#' @param w weights
#' @param title title passed to ggplot
#' @param labx name x
#' @param laby name y
#' @param labz name of legend
#' @param lims optional imposed limits for x & y axis
#' @param alpha logical, control transparency
#' @param thr_w Threshold of weights above which values are plotted. If \code{NULL}, plot weights above \code{qt} percentile. Set to \code{-Inf} to plot everything
#' @param qt probability for percentile threshold [default=0.975]
#' @param regline logical, calculate and plot regression line
#' @param type bubble (plot) or hexplot
#' @param FUN function for \code{\link[ggplot2]{stat_summary_hex}}
#' @param nbins for hexplot
#' @param collo colour low level
#' @param colhi colour low level
#' @param diverging diverging colour gradient with center 0
#' @param verbose details
#' @import ggplot2
#' @importFrom stats lm quantile coef
#' @export
plot_w_reg <- function(df, x = "donor_v", y = "recip_v", w = "ip",
                       title = "", labx = "", laby = "", labz= "",
                       lims = NULL, alpha = TRUE,
                       thr_w = NULL, qt = 0.975, regline = TRUE,
                       FUN = function(x) high_contrast(sum(x)),
                       nbins = 30,
                       collo = "white", colhi = "black",
                       diverging = FALSE,
                       type = c("bubble", "hexplot"),
                       verbose = TRUE){
  type <- match.arg(type)
  if(is.null(thr_w)){
    thr_w <- quantile(df[, w], probs = qt)
    if(verbose) print(paste("plot ip >", thr_w))
  }
  ##- regression (with all values)
  if(regline){
    coefs <- coef(lm(df[, y] ~ df[, x], weights = df[, w]))
    print(coefs)
  }
  print(title)
  ##- subset
  df <- df[df[, w] > thr_w, ]
  ka <- ifelse(alpha, nrow(df)^(1/6), 1) # transparency factor # log10(nrow(df))
  if(is.null(lims)){
    limx <- c(min(df[, x], na.rm = T), max(df[, x], na.rm = T))
    limy <- c(min(df[, y], na.rm = T), max(df[, y], na.rm = T))
  } else limx <- limy <- lims
  ##- plot with selection of minimum ip
  if (type == "bubble"){
    p <- ggplot(df, aes_string(x = x, y = y, size = w)) +
      geom_point(alpha = 1/ka) +
      scale_size_area(max_size = 10)
  } else if (type == "hexplot"){
    p <- ggplot(data = df, aes_string(x = x, y = y)) +
      stat_summary_hex(aes_string(z = w), #stat_summary_2d
                       bins = nbins, alpha = 1,
                       fun = FUN )
      #fun = function(x) sum(high.contrast(x)) )+
      #fun = sum) +
      #fun = function(x) sum(log10(x))) +
      #fun = function(x) mean(log10(x))) +
      #fun = median) +
      # + theme(legend.position="none")
  } else stop('Need type bubble or hexplot')
 if(diverging){
   p <- p + scale_fill_gradient2()
 } else {
   p <- p + scale_fill_continuous(low = collo, high = colhi, name = labz)
   }
  p <- p +
  labs(title = title, x = labx, y = laby, size = w) +
    xlim(limx) + ylim(limy) +
    theme_linedraw() #+ theme_light() # + theme_classic()
  if (regline){p <- p +
    geom_abline(intercept = coefs[1], slope = coefs[2], colour = 'blue') +
    geom_abline(intercept = 0, slope = 1, size = .2, linetype = 2)
  }
  return(p)
}

##---- plot_list ----
#' Make a list of plots and arrange in a grid
#' @param lst Named list of dataframes. Names are passed to \code{title =} of \code{funplot}
#' @param funplot Plotting function
#' @param ncols Number of columns if multiple plots in grid [default = 2]
#' @param ... further options to pass to \code{funplot}
#' @export
plot_list <- function(lst, funplot, ncols = 2, ...){
  list_of_plots <- lapply(1:length(lst), function(i) funplot(lst[[i]], title = names(lst)[i], ...))
  do.call(gridExtra::grid.arrange, c( list_of_plots, ncol = ifelse(length(list_of_plots) == 1, 1, ncols) ))
}

##---- mat_to_bar ----
#' Bar plot from mixing matrix
#'
#' Representation of matrices standardized by row / column in horizontal / vertical bar density
#' @param m matrix
#' @param title title
#' @param x_lab rows label
#' @param y_lab columns label
#' @param z_lab values label
#' @param horiz horizontal TRUE (default), vertical FALSE
#' @param values show values in or above bar, in black or white [default=FALSE]
#' @param rotate_y y-axis labels vertical for vertical plot [default=TRUE]
#' @param abbrev optionally pass \code{abbreviate} to axes labels [default=identity]
#' @importFrom stats as.formula
#' @examples
#' k <- 4
#' m <- matrix(1:k^2/k^2, nrow = k, byrow = TRUE)
#' rownames(m) <- paste0("R", 1:k)
#' colnames(m) <- paste0("C", 1:k)
#' m
#' mat_to_bar(m = m, horiz = TRUE, values = TRUE)
#' mat_to_bar(m = m, horiz = FALSE, values = TRUE)
#' @import ggplot2
#' @export
mat_to_bar <- function(m, title = NULL,
                       x_lab = "recipient",
                       y_lab = "donor",
                       z_lab = "ip",
                       horiz = TRUE,
                       values = FALSE,
                       rotate_y = TRUE,
                       abbrev = identity){
  ##- dimnames
  if(is.null(rownames(m))) rownames(m) <- paste0("R", 1:nrow(m))
  if(is.null(colnames(m))) colnames(m) <- paste0("C", 1:ncol(m))
  dimnames(m) <- lapply(dimnames(m), abbrev)
  ##- melt
  df <- as.data.frame(m)
  df[, x_lab] <- factor(rownames(m), levels = rev(rownames(m))) # rownames(m)
  d <- reshape2::melt(df, id.vars = x_lab, variable.name = y_lab, value.name = z_lab)

  if(horiz){
    ## for values, arrange matrix
    ## flip matrix
    flip_mat <- function(m){
      apply(m, 2, rev)
    }
    vals <- as.vector(m) # as.vector(t( flip_mat(m) ))
    jus <- ifelse(vals < max(vals)/3, -0.5, 1.5)
    cols <- ifelse(vals < max(vals)/3, "black", "white")

    p  <-  ggplot(data = d, aes_string(x = y_lab, y = z_lab)) +
      theme_minimal()
    p <- p + geom_bar(stat = "identity")
    p <- p + facet_grid(stats::as.formula(paste(x_lab, "~ .")), switch = "y") +
      theme(strip.background = element_blank()) +
      scale_y_continuous(breaks = NULL, name = x_lab)
    if(values) p <- p + geom_text(aes(label = round(vals*100)),
                                  vjust = jus, color = cols)
  } else {
    ## for values, arrange matrix
    vals <- as.vector( m )
    jus <- ifelse(vals < max(vals)/3, -0.5, 1.5)
    cols <- ifelse(vals < max(vals)/3, "black", "white")

    p  <-  ggplot(data = d, aes_string(x = x_lab, y = z_lab)) +
      theme_minimal()
    p <- p + geom_bar(stat = "identity") + coord_flip()
    p <- p + facet_grid(as.formula(paste(". ~", y_lab)), switch = "x") +
      theme(strip.background = element_blank()) +
      scale_y_continuous(breaks = NULL, name = y_lab) + # hack for x name of strip
      scale_x_discrete(limits = rev(levels(df[, x_lab])))
    if(rotate_y){p <- p +  theme(axis.text.y = element_text(angle = 90, hjust = 0.5))}
    if(values) p <- p + geom_text(aes(label = round(vals*100)),
                                  hjust = jus, color = cols)
  }
  return(p + ggtitle(title)) # theme(aspect.ratio=1)
}
