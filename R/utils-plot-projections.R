#' Projection plot
#'
#' Produce the graph of the observed and projected age-standardized rates and numbers
#'
#' @inheritParams canproj
#' @param site.asr Data frame with age-standardized rates and total numbers.
#' @param sex Specify (`F` or `M`) for defining colour used in graph.
#' @param ma Coefficient for left y-axis scale.
#' @param mr Coefficient for right y-axis scale.
#' @param starty Starting calendar year for observation.
#'
#' @keywords internal
projplot <- function(
  site.asr,
  sex = NULL,
  ma = 2,
  mr = 1.02,
  starty = 1986,
  startp = 2011
) {
  case <- site.asr[, 2]
  asr <- site.asr[, 1]
  n <- dim(site.asr)[1] #total number of years
  n1 <- startp - starty #n1 in observations
  n2 <- n - n1 #n2 in projections
  maxn <- max(case)
  maxr <- max(asr)
  if (is.null(sex)) {
    mycol <- c("red4", "red2", "red4", "red2")
  } else if (sex == "M") {
    mycol <- c("darkblue", "blue", "darkblue", "blue2")
  } else {
    mycol <- c("hotpink", "lightpink", "hotpink", "lightpink2")
  }
  graphics::par(mar = c(4, 4, 2, 4), mgp = c(3, 1, 0), cex = 0.7)
  graphics::barplot(
    c(case[1:n1], rep(0, n2)),
    col = mycol[1],
    space = 1,
    ylim = c(0, ma * maxn),
    cex.axis = 1,
    xlab = "",
    las = 1,
    xaxs = "r",
    border = NA
  )
  graphics::abline(h = 0)
  graphics::par(mgp = c(3, 1, 0))
  graphics::barplot(
    c(rep(0, n1), case[(n1 + 1):n]),
    add = T,
    col = mycol[2],
    space = 1,
    cex.axis = 1,
    las = 1,
    xaxs = "r",
    border = NA
  )
  graphics::par(new = T, mar = c(4, 4, 2, 4), mgp = c(3, 1, 0), cex = 1)
  graphics::plot(
    c(starty:(starty + n - 1)),
    c(asr[1:n1], rep("", n2)),
    type = "o",
    pch = 20,
    ylim = c(0, maxr),
    xaxt = "n",
    yaxt = "n",
    bty = "n",
    lwd = 3,
    ylab = "",
    col = mycol[3],
    xlab = "",
    cex.axis = 1
  )
  graphics::lines(
    c(starty:(starty + n - 1)),
    c(rep("", (n1 - 1)), asr[n1:n]),
    type = "l",
    pch = 20,
    ylim = c(0, maxr),
    xaxt = "n",
    yaxt = "n",
    bty = "n",
    lwd = 3,
    lty = 2,
    col = mycol[4],
    cex.axis = 1
  )
  graphics::axis(
    4,
    at = seq(0, (mr * maxr), ceiling(maxr / 6)),
    las = 1,
    cex.axis = 1
  )
  graphics::mtext(side = 4, "", line = 2.5, cex = 1)
}
