# Implement patchwork bugfix in older version here
# From https://github.com/thomasp85/patchwork/issues/325

add_strips <- function(gt) {
  panel_loc <- find_panel(gt)
  strip_pos <- switch(
    find_strip_pos(gt),
    inside = 0,
    outside = 2
  )
  if (!any(grepl('strip-b', gt$layout$name))) {
    gt <- gtable_add_rows(gt, unit(0, 'mm'), panel_loc$b + strip_pos)
  } else if (strip_pos == 2 && !any(gt$layout$b == panel_loc$b + 2)) {
    # Merge the strip-gap height into the axis and remove it. Only performed if
    # an axis exist
    gt$heights[panel_loc$b + 1] <- sum(gt$heights[panel_loc$b + c(1, 2)])
    gt <- gt[-(panel_loc$b + 2), ]
  }
  if (!any(grepl('strip-t', gt$layout$name))) {
    gt <- gtable_add_rows(gt, unit(0, 'mm'), panel_loc$t - 1 - strip_pos)
  } else if (strip_pos == 2 && !any(gt$layout$t == panel_loc$t - 2)) {
    gt$heights[panel_loc$t - 1] <- sum(gt$heights[panel_loc$t - c(1, 2)])
    gt <- gt[-(panel_loc$t - 2), ]
  }
  if (!any(grepl('strip-r', gt$layout$name))) {
    gt <- gtable_add_cols(gt, unit(0, 'mm'), panel_loc$r + strip_pos)
  } else if (strip_pos == 2 && !any(gt$layout$r == panel_loc$r + 2)) {
    gt$widths[panel_loc$r + 1] <- sum(gt$widths[panel_loc$r + c(1, 2)])
    gt <- gt[, -(panel_loc$r + 2)]
  }
  if (!any(grepl('strip-l', gt$layout$name))) {
    gt <- gtable_add_cols(gt, unit(0, 'mm'), panel_loc$l - 1 - strip_pos)
  } else if (strip_pos == 2) {
    gt$widths[panel_loc$l - 1] <- sum(gt$widths[panel_loc$l - c(1, 2)])
    gt <- gt[, -(panel_loc$l - 2)]
  }
  gt
}

environment(add_strips) <- asNamespace("patchwork")
assignInNamespace("add_strips", add_strips, ns = "patchwork")


# From modelsummary:::fmt_labels_md
# Authored by Lukas Wallrich

fmt_labels_md <- function (tab, position = c("both", "row", "column"))
{
  out <- tab
  if (match.arg(position) %in% c("both", "row")) {
    out <- gt::fmt_markdown(out, columns = 1)
  }
  if (match.arg(position) %in% c("both", "column")) {
    f <- function(x) stats::setNames(lapply(names(x$`_data`),
                                            gt::md), names(x$`_data`))
    out <- gt::cols_label(out, .list = f(out))
  }
  return(out)
}


# From https://raw.githubusercontent.com/MathiasHarrer/dmetar/master/R/mlm.variance.distribution.R
mlm.variance.distribution = var.comp = function(x){

  m = x

  # Check class
  if (!(class(m)[1] %in% c("rma.mv", "rma"))){
    stop("x must be of class 'rma.mv'.")
  }

  # Check for three level model
  if (m$sigma2s != 2){
    stop("The model you provided does not seem to be a three-level model. This function can only be used for three-level models.")
  }

  # Check for right specification (nested model)
  if (sum(grepl("/", as.character(m$random[[1]]))) < 1){
    stop("Model must contain nested random effects. Did you use the '~ 1 | cluster/effect-within-cluster' notation in 'random'? See ?metafor::rma.mv for more details.")
  }

  # Get variance diagonal and calculate total variance
  n = m$k.eff
  vector.inv.var = 1/(diag(m$V))
  sum.inv.var = sum(vector.inv.var)
  sum.sq.inv.var = (sum.inv.var)^2
  vector.inv.var.sq = 1/(diag(m$V)^2)
  sum.inv.var.sq = sum(vector.inv.var.sq)
  num = (n-1)*sum.inv.var
  den = sum.sq.inv.var - sum.inv.var.sq
  est.samp.var = num/den

  # Calculate variance proportions
  level1=((est.samp.var)/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
  level2=((m$sigma2[2])/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
  level3=((m$sigma2[1])/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)

  # Prepare df for return
  Level=c("Level 1", "Level 2", "Level 3")
  Variance=c(level1, level2, level3)
  df.res=data.frame(Variance)
  colnames(df.res) = c("% of total variance")
  rownames(df.res) = Level
  I2 = c("---", round(Variance[2:3], 2))
  df.res = as.data.frame(cbind(df.res, I2))

  totalI2 = Variance[2] + Variance[3]


  # Generate plot
  df1 = data.frame("Level" = c("Sampling Error", "Total Heterogeneity"),
                   "Variance" = c(df.res[1,1], df.res[2,1]+df.res[3,1]),
                   "Type" = rep(1,2))

  df2 = data.frame("Level" = rownames(df.res),
                   "Variance" = df.res[,1],
                   "Type" = rep(2,3))

  df = as.data.frame(rbind(df1, df2))


  g = ggplot(df, aes(fill=Level, y=Variance, x=as.factor(Type))) +
    coord_cartesian(ylim = c(0,1), clip = "off") +
    geom_bar(stat="identity", position="fill", width = 1, color="black") +
    scale_y_continuous(labels = scales::percent)+
    theme(axis.title.x=element_blank(),
          axis.text.y = element_text(color="black"),
          axis.line.y = element_blank(),
          axis.title.y=element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_line(lineend = "round"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.background = element_rect(linetype="solid",
                                           colour ="black"),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"cm"),
          axis.ticks.length=unit(.25, "cm"),
          plot.margin = unit(c(1,3,1,1), "lines")) +
    scale_fill_manual(values = c("darkseagreen3", "deepskyblue3", "darkseagreen2",
                                 "deepskyblue1", "deepskyblue2")) +

    # Add Annotation

    # Total Variance
    annotate("text", x = 1.5, y = 1.05,
             label = paste("Total Variance:",
                           round(m$sigma2[1]+m$sigma2[2]+est.samp.var, 3))) +

    # Sampling Error
    annotate("text", x = 1, y = (df[1,2]/2+df[2,2])/100,
             label = paste("Sampling Error Variance: \n", round(est.samp.var, 3)), size = 3) +

    # Total I2
    annotate("text", x = 1, y = ((df[2,2])/100)/2-0.02,
             label = bquote("Total"~italic(I)^2*":"~.(round(df[2,2],2))*"%"), size = 3) +
    annotate("text", x = 1, y = ((df[2,2])/100)/2+0.05,
             label = paste("Variance not attributable \n to sampling error: \n", round(m$sigma2[1]+m$sigma2[2],3)), size = 3) +

    # Level 1
    annotate("text", x = 2, y = (df[1,2]/2+df[2,2])/100, label = paste("Level 1: \n",
                                                                       round(df$Variance[3],2), "%", sep=""), size = 3) +

    # Level 2
    annotate("text", x = 2, y = (df[5,2]+(df[4,2]/2))/100,
             label = bquote(italic(I)[Level2]^2*":"~.(round(df[4,2],2))*"%"), size = 3) +

    # Level 3
    annotate("text", x = 2, y = (df[5,2]/2)/100,
             label = bquote(italic(I)[Level3]^2*":"~.(round(df[5,2],2))*"%"), size = 3)

  returnlist = list(results = df.res,
                    totalI2 = totalI2,
                    plot = g)
  class(returnlist) = c("mlm.variance.distribution", "list")

  invisible(returnlist)

  returnlist

}

# From lukaswallrich/rNuggets

#' Create mice predictorMatrix with relevant factors (extends mice::quickpred)
#'
#' \code{\link[mice]{quickpred}} creates a predictor matrix that suggests which variables should
#' be used for multiple imputation. However, it does not yield meaningful information
#' for (unordered) factors. This function returns the same as mice::quickpred for
#' numeric and logical variables and binary factors, but excludes character
#' vectors and tests the predictive power of each level of factors with more
#' than two levels separately (through dummy coding). If the correlation exceeds `mincor` for at least one
#' level, the factor variable is used as a predictor.
#'
#' @param data Matrix or data frame with incomplete data.
#' @param mincor A scalar, numeric vector (of size \code{ncol(data))} or numeric
#' matrix (square, of size \code{ncol(data)} specifying the minimum
#' threshold(s) against which the absolute correlation in the data is compared.
#' @param minpuc A scalar, vector (of size \code{ncol(data))} or matrix (square,
#' of size \code{ncol(data)} specifying the minimum threshold(s) for the
#' proportion of usable cases.
#' @param include A string or a vector of strings containing one or more
#' variable names from \code{names(data)}. Variables specified are always
#' included as a predictor.
#' @param exclude A string or a vector of strings containing one or more
#' variable names from \code{names(data)}. Variables specified are always
#' excluded as a predictor.
#' @param method A string specifying the type of correlation. Use
#' \code{'pearson'} (default), \code{'kendall'} or \code{'spearman'}. Can be
#' abbreviated.
#' @return A square binary matrix of size \code{ncol(data)}.
#' @source This function is based on the mice::quickpred function as available in mice v 3.11.4, written by Stef van Buuren. The code for dummy coding factor variables is based on psych::dummy.code, written by William Revelle
#' @seealso \code{\link{mice::quickpred}}
#' @export

quickpred_ext <- function (data, mincor = 0.1, minpuc = 0, include = "", exclude = "",
                           method = "pearson")
{
  data <- mice:::check.dataform(data)
  nvar <- ncol(data)

  #Identify characters and factors
  chr <- names(data)[sapply(data, is.character)]
  chr_which <- NULL

  if (length(chr) > 0) {
    message(paste0("Data contains character variable(s): ", paste(chr, collapse = " "),
                   ". These will not be used as predictors."))
    chr_which <- which(names(data) %in% chr)
  }

  fct <- names(data)[sapply(data, nlevels )>2]
  fct_which <- NULL

  if (length(fct) > 0) {
    message(paste0("Data contains factor variable(s) with more than 2 levels: ", paste(fct, collapse = " "),
                   ". mincor argument will be tested for each level."))
    fct_which <- which(names(data) %in% fct)
  }

  # initialize
  predictorMatrix <- matrix(0, nrow = nvar, ncol = nvar, dimnames = list(names(data),
                                                                         names(data)))
  x <- data.matrix(data)
  r <- !is.na(x)

  # include predictors with
  # 1) pairwise correlation among data
  ## Without factors
  if (length(fct) == 0) {

    suppressWarnings(v <- abs(stats::cor(x, use = "pairwise.complete.obs",
                                         method = method)))
    v[is.na(v)] <- 0
  } else {
    ## Correction for factors
    x_dummies <- x

    fct_start <- c()
    fct_end <- c()

    for (i in seq_along(fct_which)) {
      fct_start[i] <- ncol(x_dummies) + 1
      d <- dummy_code(x[,fct_which[i]])
      x_dummies <- cbind(x_dummies, d)
      fct_end[i] <- ncol(x_dummies)
    }
    suppressWarnings(v <- abs(stats::cor(x_dummies, use = "pairwise.complete.obs",
                                         method = method)))

    for (i in seq_along(fct_which)) {
      v[,fct_which[i]] <- do.call(pmax, c(data.frame(v[, fct_start[i]:fct_end[i]]), na.rm = TRUE))
    }
    v <- v[1:ncol(x), 1:ncol(x)]
    v[is.na(v)] <- 0
  }
  # 2) pairwise correlation of data with response indicator higher than mincor
  suppressWarnings(u <- abs(stats::cor(y = x, x = r, use = "pairwise.complete.obs",
                                       method = method)))
  u[is.na(u)] <- 0
  maxc <- pmax(v, u)
  predictorMatrix[maxc > mincor] <- 1

  # exclude predictors with a percentage usable cases below minpuc
  p <- mice::md.pairs(data)
  puc <- p$mr/(p$mr + p$mm)
  predictorMatrix[puc < minpuc] <- 0

  #exclude character vars
  predictorMatrix[, chr_which] <- 0

  # exclude predictors listed in the exclude argument
  yz <- pmatch(exclude, names(data))
  predictorMatrix[, yz] <- 0

  # include predictors listed in the include argument
  yz <- pmatch(include, names(data))
  predictorMatrix[, yz] <- 1

  # some final processing
  diag(predictorMatrix) <- 0
  predictorMatrix[colSums(!r) == 0, ] <- 0

  return(predictorMatrix)
}

#' Dummy code variable
#'
#' Simplified from psych::dummy.code and
#' changed to return data.frame
#'
#' @param x A vector to be transformed into dummy codes.
#' @source Simplified from psych::dummy.code, written by William Revelle


dummy_code <- function(x) {
  t <- table(x)
  lt <- length(t)
  n.obs <- length(x)
  new <- matrix(0,nrow=n.obs,ncol=lt)
  new[is.na(x),] <- NA
  xlev <- factor(x,levels=names(t))

  for (i in 1:n.obs) {new[i,xlev[i]] <- 1}

  as.data.frame(new)
}


for(i in seq_along(names(df))) {
  n_rows <- nrow(df)
  perc_missing <- 15 # percentage missing data
  row_missing <- sample(1:n_rows, sample(1:n_rows, round(perc_missing/100 * n_rows,0))) # sample randomly x% of rows
  col_missing <- i # define column
  df[row_missing, col_missing] <- NA # assign missing values
}
