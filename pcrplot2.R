lmmatrix1<-function (y, cov) {
  p <- matrix(NA, ncol = ncol(cov), nrow = ncol(y))
  for (j in 1:ncol(cov)) {
    x <- cov[, j]
    for (i in 1:ncol(y)) {
      fit <- summary(lm(y[, i] ~ x, na.action = na.omit))
      f <- fit$fstatistic
      p[i, j] <- pf(f["value"], f["numdf"],
                    f["dendf"], lower.tail = FALSE)
    }
  }
  colnames(p) <- names(cov)
  p
}
pcrplot2<-function (beta, cov, npc = 50)
{
  if (!is.matrix(beta)) {
    stop("beta is not a data matirx")
  }
  if (!is.data.frame(cov)) {
    stop("cov is not a data frame")
  }
  if (ncol(beta) != nrow(cov)) {
    stop("number of columns in beta is not equal to number of rows in cov")
  }
  cat("Analysis is running, please wait...!", "\n")
  npc <- min(ncol(beta), npc)
  svd <- prcomp(t(beta), center = TRUE, scale = TRUE, retx = TRUE)#,
  #na.action = "na.omit")
  screeplot(svd, npc, type = "barplot")
  eigenvalue <- svd[["sdev"]]^2
  prop <- (sum(eigenvalue[1:npc])/sum(eigenvalue)) * 100
  cat("Top ", npc, " principal components can explain ",
      prop, "% of data \n    variation", "\n")
  prop <- (sum(eigenvalue[1:1])/sum(eigenvalue)) * 100
  cat("Top ", "1", " principal components can explain ",
      prop, "% of data \n    variation", "\n")
  prop <- (sum(eigenvalue[1:2])/sum(eigenvalue)) * 100
  cat("Top ", "2", " principal components can explain ",
      prop, "% of data \n    variation", "\n")
  prop <- (sum(eigenvalue[1:3])/sum(eigenvalue)) * 100
  cat("Top ", "3", " principal components can explain ",
      prop, "% of data \n    variation", "\n")
  prop <- (sum(eigenvalue[1:4])/sum(eigenvalue)) * 100
  cat("Top ", "4", " principal components can explain ",
      prop, "% of data \n    variation", "\n")
  p <- lmmatrix1(svd$x[, 1:npc], cov)
  yaxis <- colnames(p)
  title = "Principal Component Regression Analysis"
  xmax=npc
  plot(1, xlim = c(0, xmax), ylim = c(0, length(yaxis) + 1),
       type = "n", bty = "n", axes = FALSE, xlab = "Principal Component",
       ylab = "", main = title)
  axis(1, at = c(1:xmax), pos = 0.5, las = 1, lwd = 3)
  for (i in 1:length(yaxis)) {
    text(0.3, i, yaxis[i], xpd = TRUE, adj = 1)
  }
  for (i in 1:ncol(p)) {
    for (j in 1:nrow(p)) {
      pp <- p[j, i]
      colcode <- "white"
      if (pp <= 1e-09) {
        colcode = "darkred"
      }
      else if (pp <= 1e-04) {
        colcode = "red"
      }
      else if (pp <= 0.01) {
        colcode = "orange"
      }
      else if (pp <= 0.05) {
        colcode = "pink"
      }
      polygon(c(j - 0.5, j - 0.5, j + 0.5, j + 0.5), c(i -
                                                         0.5, i + 0.5, i + 0.5, i - 0.5), col = colcode,
              border = NA)
    }
  }
  legend("topright", c("<0.05", "<0.01",
                       "<10E-5", "<10E-10"), col = c("pink",
                                                     "orange", "red", "darkred"), pch = 15,
         pt.cex = 2, bty = "o", horiz = TRUE, xpd = TRUE)
}

