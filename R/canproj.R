#' CANPROJ
#'
#' Cancer incidence and mortality projections.
#'
#' @param cdat 19(age groups)*N(years) historical cancer data, 15<=N<=125.
#' @param pdat 19(age groups)*(N+M)(years) observed and projected population, 5<=M<=25.
#' @param ncase minimum number of cancer cases/deaths per year for splitting data.
#' @param nagg number of years for data aggregation (by years), default: 1-annual data
#' @param startestage user defined age groups for modeling can be input here
#' @param cuttrd Degenerating percent of trends per year after 5 years (shortp=0) or the first projection year.
#' @param projfor Specify "incidence" or "mortality" if want ASR as criteria for nagg.
#' @param newcohort assign new cohort effect as 0 (NULL) or the last estimated
#'   cohort effect (T), default is 0, use "T" only if having evidence on negative new cohort effect.
#' @param Ave5 Ave5=T invokes the 5 year average method when age-only model is selected.
#' @param sum5 When the 5-year average method is used, sum5=NULL call the 5-year
#'   period based rate,otherwise, average the 5 rates in the 5 years for each age group.
#' @param methods user required projection method can be specified by ADPC
#'   models: "nordpred" or "adpc-nb"; age-cohort models: "ac", "ac-nb"; age-period models: "age-trd", "com-trd"; and average: "ave5".
#' @param linkfunc Link function, default is power5, can be log, sqrt and identity.
#' @param pD Trend selecting criteria of p-value of drift (linear trend) term.
#' @param pGOF Model selection criteria of p-value of goodness-of-fit.
#' @param standpop The weights (proportions) of 19 age groups in a standard population.
#' @param startp The start calendar year of projection, e.g. 2009.
#'
#' @md
#'
#' @return A `list()`.
#'
#' @export
canproj <- function(
  cdat,
  pdat,
  startp,
  projfor = "incidence",
  nagg = NULL,
  ncase = NULL,
  startestage = NULL,
  newcohort = NULL,
  Ave5 = NULL,
  sum5 = NULL,
  methods = NULL,
  linkfunc = "power5",
  cuttrd = 0.04,
  shortp = 0,
  pD = 0.05,
  pGOF = 0.05,
  standpop = ca11
) {
  # Check data:
  if (dim(cdat)[1] != 19 || dim(pdat)[1] != 19) {
    stop(
      "\"cdat\" and \"pdat\" must have data for 19 age groups by row in this version"
    )
  }
  if (dim(cdat)[2] > dim(pdat)[2]) {
    stop("\"pdat\" must include information about all years in \"cdat\"")
  }
  if (dim(pdat)[2] == dim(cdat)[2]) {
    stop("\"pdat\" must include information in projection years")
  }
  if ((dim(pdat)[2] - dim(cdat)[2]) > 30) {
    stop("Package can not project more than 30 years")
  }
  ## The five year average method is used due to small annual number
  if (mean(apply(cdat, 2, sum)[(dim(cdat)[2] - 9):(dim(cdat)[2])]) < 3) {
    out <- ave5proj(cdat, pdat, startp = startp, sum5 = sum5)
    outasp <- ave5proj.getproj(pdat, out, standpop = NULL)
    outann <- ave5proj.getproj(pdat, out, standpop = standpop)
    mod <- "average5"
    pdPC <- NA
  } else {
    # Define number of cases for data splitting:
    if (is.null(ncase)) {
      if (projfor == "incidence") {
        ncase <- 1
      } else {
        ncase <- 3 / 5
      }
    }

    # Define number of years for data aggregation:
    if (is.null(nagg)) {
      masr <- mean(obasr(cdat, pdat)[, 1])
      if (projfor == "incidence") {
        if (masr > 15) {
          nagg <- 1
        } else if (masr > 10) {
          nagg <- 2
        } else if (masr > 5) {
          nagg <- 3
        } else {
          nagg <- 4
        }
      } else {
        if (masr > 9) {
          nagg <- 1
        } else if (masr > 6) {
          nagg <- 2
        } else if (masr > 3) {
          nagg <- 3
        } else {
          nagg <- 4
        }
      }
    }

    ## aggregating data by 5 years:
    aggdata <- datagg(cdat, pdat, 5)
    cases <- aggdata$cases
    pyr <- aggdata$pyr

    nototper <- dim(pyr)[2]
    noobsper <- dim(cases)[2]
    nonewpred <- nototper - noobsper

    # Preparing parameters for APC method:
    n5case <- 5 * ncase
    cuttrend <- rep(shortp, nonewpred)
    for (i in 1:nonewpred) {
      cuttrend[i] <- shortp + (i - 1) * 5 * cuttrd
    }
    cuttrend[cuttrend > 1] <- 1

    # when projection method is not specified:
    if (is.null(methods)) {
      # Performing adpc modeling:
      mod0 <- adpcproj(
        cdat,
        pdat,
        projfor = projfor,
        n5case = n5case,
        newcohort = newcohort,
        pGOF = 0,
        cuttrd = 0.04,
        shortp = 0,
        linkfunc = linkfunc
      )
      pv1 <- summary(mod0$glm)$coef["Period", 4]
      options(warn = -1)
      pv2 <- drop1(mod0$glm, test = "F")["as.factor(Cohort)", 5]
      pv3 <- drop1(mod0$glm, test = "F")["as.factor(Period)", 5]
      options(warn = 0)
      pdPC <- c(pv1, pv3, pv2)
      # Selecting method by decision tree:
      if (pv2 < pGOF) {
        if (pv1 < pD | pv3 < pGOF) {
          out <- adpcproj(
            cdat,
            pdat,
            projfor = projfor,
            n5case = n5case,
            startestage = startestage,
            newcohort = newcohort,
            pGOF = pGOF,
            cuttrd = 0.04,
            shortp = 0,
            linkfunc = linkfunc
          )
          r0 <- adpcproj.getpred(out, incidence = T)
          outasp <- asrpy(r0, cdat, pdat, startp = startp, nagg = 5)
          outann <- asry(outasp, pdat, standpop = standpop)
          mod <- "ADPC"
        } else {
          out <- acproj(
            cdat,
            pdat,
            projfor = projfor,
            n5case = n5case,
            startestage = startestage,
            cuttrd = 0.04,
            shortp = 0,
            pGOF = pGOF,
            linkfunc = linkfunc
          )
          outasp <- acproj.getproj(cdat, pdat, startp = startp, out)
          outann <- acproj.getproj(
            cdat,
            pdat,
            startp = startp,
            out,
            standpop = standpop
          )
          mod <- "AC"
        }
      } else {
        out <- hybdproj(
          cdat,
          pdat,
          projfor = projfor,
          nagg = nagg,
          ncase = ncase,
          cuttrd = cuttrd,
          shortp = shortp,
          linkfunc = linkfunc,
          pD = pD,
          pGOF = pGOF
        )
        outasp <- hybdproj.getproj(
          cdat,
          pdat,
          startp = startp,
          out,
          Ave5 = Ave5,
          sum5 = sum5
        )
        outann <- asry(outasp, pdat, standpop = standpop)
        mod <- "Hybrid"
      }
    } else {
      # when specified method: "nordpred", "adpc-nb", "ac-poi", "ac-nb",
      #                        "age-trd-nb", "age-trd-poi", "com-trd", "age-only", "ave5":
      pdPC <- NA
      if (methods == "nordpred") {
        out <- adpcproj(
          cdat,
          pdat,
          projfor = projfor,
          n5case = n5case,
          startestage = startestage,
          newcohort = newcohort,
          pGOF = 0,
          cuttrd = 0.04,
          shortp = 0,
          linkfunc = linkfunc
        )
        r0 <- adpcproj.getpred(out, incidence = T)
        outasp <- asrpy(r0, cdat, pdat, startp = startp, nagg = 5)
        outann <- asry(outasp, pdat, standpop = standpop)
        mod <- "nordpred"
      }
      if (methods == "adpc-nb") {
        out <- adpcproj(
          cdat,
          pdat,
          projfor = projfor,
          n5case = n5case,
          startestage = startestage,
          newcohort = newcohort,
          pGOF = 1,
          cuttrd = 0.04,
          shortp = 0,
          linkfunc = linkfunc
        )
        r0 <- adpcproj.getpred(out, incidence = T)
        outasp <- asrpy(r0, cdat, pdat, startp = startp, nagg = 5)
        outann <- asry(outasp, pdat, standpop = standpop)
        mod <- "adpc-nb"
      }
      if (methods == "ac-poi") {
        out <- acproj(
          cdat,
          pdat,
          projfor = projfor,
          n5case = n5case,
          startestage = startestage,
          cuttrd = 0.04,
          shortp = 0,
          pGOF = 0,
          linkfunc = linkfunc
        )
        outasp <- acproj.getproj(cdat, pdat, startp = startp, out)
        outann <- acproj.getproj(
          cdat,
          pdat,
          startp = startp,
          out,
          standpop = standpop
        )
        mod <- "ac-poi"
      }
      if (methods == "ac-nb") {
        out <- acproj(
          cdat,
          pdat,
          projfor = projfor,
          n5case = n5case,
          startestage = startestage,
          cuttrd = 0.04,
          shortp = 0,
          pGOF = 1,
          linkfunc = linkfunc
        )
        outasp <- acproj.getproj(cdat, pdat, startp = startp, out)
        outann <- acproj.getproj(
          cdat,
          pdat,
          startp = startp,
          out,
          standpop = standpop
        )
        mod <- "ac-nb"
      }
      if (methods == "age-trd-nb") {
        out <- hybdproj(
          cdat,
          pdat,
          projfor = projfor,
          nagg = nagg,
          ncase = ncase,
          cuttrd = cuttrd,
          shortp = shortp,
          linkfunc = linkfunc,
          pD = 1,
          pGOF = 1
        )
        outasp <- hybdproj.getproj(
          cdat,
          pdat,
          startp = startp,
          out,
          Ave5 = Ave5
        )
        outann <- asry(outasp, pdat, standpop = standpop)
        mod <- "a-s-nb"
      }
      if (methods == "age-trd-poi") {
        out <- hybdproj(
          cdat,
          pdat,
          projfor = projfor,
          nagg = nagg,
          ncase = ncase,
          cuttrd = cuttrd,
          shortp = shortp,
          linkfunc = linkfunc,
          pD = 0,
          pGOF = 1
        )
        outasp <- hybdproj.getproj(
          cdat,
          pdat,
          startp = startp,
          out,
          Ave5 = Ave5
        )
        outann <- asry(outasp, pdat, standpop = standpop)
        mod <- "a-s-poi"
      }
      if (methods == "com-trd") {
        out <- hybdproj(
          cdat,
          pdat,
          projfor = projfor,
          nagg = nagg,
          ncase = ncase,
          cuttrd = cuttrd,
          shortp = shortp,
          linkfunc = linkfunc,
          pD = 1,
          pGOF = 0
        )
        outasp <- hybdproj.getproj(
          cdat,
          pdat,
          startp = startp,
          out,
          Ave5 = Ave5
        )
        outann <- asry(outasp, pdat, standpop = standpop)
        mod <- "c-t"
      }
      if (methods == "age-only") {
        out <- hybdproj(
          cdat,
          pdat,
          projfor = projfor,
          nagg = nagg,
          ncase = ncase,
          cuttrd = cuttrd,
          shortp = shortp,
          linkfunc = linkfunc,
          pD = 0,
          pGOF = 0
        )
        outasp <- hybdproj.getproj(
          cdat,
          pdat,
          startp = startp,
          out,
          Ave5 = Ave5,
          sum5 = sum5
        )
        outann <- asry(outasp, pdat, standpop = standpop)
        mod <- "average"
      }
      if (methods == "ave5") {
        out <- ave5proj(cdat, pdat, startp = startp, sum5 = sum5)
        outasp <- ave5proj.getproj(pdat, out, standpop = NULL)
        outann <- ave5proj.getproj(pdat, out, standpop = standpop)
        mod <- "average5"
      }
    }
  }
  obsy <- dim(cdat)[2]
  resut <- list(
    annproj = outann,
    agsproj = outasp,
    method = mod,
    out = out,
    obsy = obsy,
    pdPC = pdPC
  )
  class(resut) <- "canproj"
  attr(resut, "Call") <- sys.call()
  return(resut)
}


#' @rdname canproj
#' @export
canproj.getproj <- function(canproj.object, standpop = NULL) {
  if (is.null(standpop)) {
    return(canproj.object$agsproj)
  } else {
    return(canproj.object$annproj)
  }
}


#' @rdname canproj
#' @export
summary.canproj <- function(canproj.object) {
  summary(canproj.object$out)
  invisible(canproj.object)
}


#' @rdname canproj
#' @export
glm.canproj <- function(canproj.object) {
  summary(canproj.object$out$glm)
}


#' @rdname canproj
#' @export
plot.canproj <- function(
  canproj.object,
  standpop,
  startplot = 1,
  xlab = "Calendar Year",
  ylab = "Rates",
  main = "",
  labels = NULL,
  ylim = NULL,
  lty = c(1, 3),
  col = c(1, 1),
  new = T,
  ...
) {
  if (class(canproj.object) != "canproj") {
    stop("Variable \"canproj.object\" must be of type \"canproj\"")
  }
  # Reading & formating data:
  indat <- canproj.getproj(canproj.object, standpop = standpop)
  indata <- indat[, 1]
  indata <- indata[startplot:length(indata)]

  # Seting internal variables:
  nopredy <- length(indata) - canproj.object$obsy
  if (is.null(labels)) {
    labels <- row.names(indat)
  }
  # Create plots:
  maxx <- length(indata)
  if (new) {
    maxy <- max(indata) * (21 / 20)
    if (is.null(ylim)) {
      ylim <- c(0, maxy)
    }
    plot(c(1, maxx), ylim, type = "n", ylab = ylab, xlab = xlab, axes = F, ...)
    axis(2)
    axis(1, at = 1:maxx, labels = labels)
    box()
    title(main)
  }
  lines(
    1:(maxx - nopredy),
    indata[1:(maxx - nopredy)],
    type = "o",
    pch = 20,
    lty = lty[1],
    col = col[1],
    ...
  )
  lines(
    (maxx - nopredy):maxx,
    indata[(maxx - nopredy):maxx],
    lty = lty[2],
    col = col[2],
    ...
  )
  # Returning object as invisible
  invisible(canproj.object)
}
