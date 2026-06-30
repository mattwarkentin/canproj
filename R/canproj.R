#' CANPROJ
#'
#' Cancer incidence and mortality projections.
#'
#' @param cdat 19(age groups)*N(years) historical cancer data, 15<=N<=125.
#' @param pdat 19(age groups)*(N+M)(years) observed and projected population, 5<=M<=25.
#' @param startp The start calendar year of projection, e.g. 2009.
#' @param standpop The weights (proportions) of 19 age groups in a standard population.
#' @param projfor Specify "incidence" or "mortality" if want ASR as criteria for nagg.
#' @param nagg Number of years for data aggregation (by years), default: 1-annual data
#' @param ncase Minimum number of cancer cases/deaths per year for splitting data.
#' @param startage Youngest age group to include in the GLM. Default (`NULL`) picks based on mean cases.
#' @param newcohort Assign new cohort effect as 0 (NULL) or the last estimated
#'   cohort effect (T), default is 0, use "T" only if having evidence on negative new cohort effect.
#' @param Ave5 Ave5=T invokes the 5 year average method when age-only model is selected.
#' @param sum5 When the 5-year average method is used, sum5=NULL call the 5-year
#'   period based rate, otherwise, average the 5 rates in the 5 years for each age group.
#' @param methods User required projection method can be specified by ADPC
#'   models: "nordpred" or "adpc-nb"; age-cohort models: "ac", "ac-nb"; age-period models: "age-trd", "com-trd"; and average: "ave5".
#' @param linkfunc Link function, default is power5, can be log, sqrt and identity.
#' @param cuttrd Degenerating percent of trends per year after 5 years (shortp=0) or the first projection year.
#' @param shortp Attenuation percent of drift term or slope for the first 5 projection years.
#' @param pD Trend selecting criteria of p-value of drift (linear trend) term.
#' @param pGOF Model selection criteria of p-value of goodness-of-fit.
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
  standpop,
  projfor = "incidence",
  nagg = NULL,
  ncase = NULL,
  startage = NULL,
  newcohort = NULL,
  Ave5 = NULL,
  sum5 = NULL,
  methods = NULL,
  linkfunc = "power5",
  cuttrd = 0.04,
  shortp = 0,
  pD = 0.05,
  pGOF = 0.05
) {
  S7::check_is_S7(standpop, StandardPopulation)

  if (dim(cdat)[2] > dim(pdat)[2]) {
    stop("\"pdat\" must include information about all years in \"cdat\"")
  }
  if (dim(pdat)[2] == dim(cdat)[2]) {
    stop("\"pdat\" must include information in projection years")
  }
  if ((dim(pdat)[2] - dim(cdat)[2]) > 30) {
    stop("Package can not project more than 30 years")
  }

  if (mean(apply(cdat, 2, sum)[(dim(cdat)[2] - 9):(dim(cdat)[2])]) < 3) {
    out <- ave5proj(cdat, pdat, startp = startp, sum5 = sum5)
    outasp <- ave5proj_get_projections(pdat, out, standpop = NULL)
    outann <- ave5proj_get_projections(pdat, out, standpop = standpop)
    mod <- "average5"
    pdPC <- NA
  } else {
    if (is.null(ncase)) {
      if (projfor == "incidence") {
        ncase <- 1
      } else {
        ncase <- 3 / 5
      }
    }

    if (is.null(nagg)) {
      nagg <- select_nagg(cdat, pdat, standpop, projfor)
    }

    aggdata <- aggregate_data(cdat, pdat, 5)
    cases <- aggdata$cases
    pyr <- aggdata$pyr
    nonewpred <- ncol(pyr) - ncol(cases)

    n5case <- 5 * ncase
    cuttrend <- rep(shortp, nonewpred)
    for (i in 1:nonewpred) {
      cuttrend[i] <- shortp + (i - 1) * 5 * cuttrd
    }
    cuttrend[cuttrend > 1] <- 1

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
      pv2 <- suppressWarnings(stats::drop1(mod0$glm, test = "F")[
        "as.factor(Cohort)",
        5
      ])
      pv3 <- suppressWarnings(stats::drop1(mod0$glm, test = "F")[
        "as.factor(Period)",
        5
      ])

      pdPC <- c(pv1, pv3, pv2)

      if (pv2 < pGOF) {
        if (pv1 < pD | pv3 < pGOF) {
          out <- adpcproj(
            cdat,
            pdat,
            projfor = projfor,
            n5case = n5case,
            startage = startage,
            newcohort = newcohort,
            pGOF = pGOF,
            cuttrd = 0.04,
            shortp = 0,
            linkfunc = linkfunc
          )

          r0 <- adpcproj.getpred(out, incidence = T)
          outasp <- interpolate_age_specific_rates(
            r0,
            cdat,
            pdat,
            startp = startp,
            nagg = 5
          )
          outann <- standardize_annual_rates(outasp, pdat, stdpop = standpop)
          mod <- "ADPC"
        } else {
          out <- acproj(
            cdat,
            pdat,
            projfor = projfor,
            n5case = n5case,
            startage = startage,
            cuttrd = 0.04,
            shortp = 0,
            pGOF = pGOF,
            linkfunc = linkfunc
          )

          outasp <- acproj_get_projections(cdat, pdat, startp = startp, out)
          outann <- acproj_get_projections(
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
          standpop,
          projfor = projfor,
          nagg = nagg,
          ncase = ncase,
          cuttrd = cuttrd,
          shortp = shortp,
          linkfunc = linkfunc,
          pD = pD,
          pGOF = pGOF
        )

        outasp <- hybdproj_get_projections(
          cdat,
          pdat,
          startp = startp,
          out,
          Ave5 = Ave5,
          sum5 = sum5
        )

        outann <- standardize_annual_rates(outasp, pdat, stdpop = standpop)
        mod <- "Hybrid"
      }
    } else {
      pdPC <- NA

      if (methods == "nordpred") {
        out <- adpcproj(
          cdat,
          pdat,
          projfor = projfor,
          n5case = n5case,
          startage = startage,
          newcohort = newcohort,
          pGOF = 0,
          cuttrd = 0.04,
          shortp = 0,
          linkfunc = linkfunc
        )

        r0 <- adpcproj.getpred(out, incidence = T)
        outasp <- interpolate_age_specific_rates(
          r0,
          cdat,
          pdat,
          startp = startp,
          nagg = 5
        )
        outann <- standardize_annual_rates(outasp, pdat, stdpop = standpop)
        mod <- "nordpred"
      } else if (methods == "adpc-nb") {
        out <- adpcproj(
          cdat,
          pdat,
          projfor = projfor,
          n5case = n5case,
          startage = startage,
          newcohort = newcohort,
          pGOF = 1,
          cuttrd = 0.04,
          shortp = 0,
          linkfunc = linkfunc
        )

        r0 <- adpcproj.getpred(out, incidence = T)
        outasp <- interpolate_age_specific_rates(
          r0,
          cdat,
          pdat,
          startp = startp,
          nagg = 5
        )
        outann <- standardize_annual_rates(outasp, pdat, stdpop = standpop)
        mod <- "adpc-nb"
      } else if (methods == "ac-poi") {
        out <- acproj(
          cdat,
          pdat,
          projfor = projfor,
          n5case = n5case,
          startage = startage,
          cuttrd = 0.04,
          shortp = 0,
          pGOF = 0,
          linkfunc = linkfunc
        )

        outasp <- acproj_get_projections(cdat, pdat, startp = startp, out)
        outann <- acproj_get_projections(
          cdat,
          pdat,
          startp = startp,
          out,
          standpop = standpop
        )
        mod <- "ac-poi"
      } else if (methods == "ac-nb") {
        out <- acproj(
          cdat,
          pdat,
          projfor = projfor,
          n5case = n5case,
          startage = startage,
          cuttrd = 0.04,
          shortp = 0,
          pGOF = 1,
          linkfunc = linkfunc
        )

        outasp <- acproj_get_projections(cdat, pdat, startp = startp, out)
        outann <- acproj_get_projections(
          cdat,
          pdat,
          startp = startp,
          out,
          standpop = standpop
        )
        mod <- "ac-nb"
      } else if (methods == "age-trd-nb") {
        out <- hybdproj(
          cdat,
          pdat,
          standpop,
          projfor = projfor,
          nagg = nagg,
          ncase = ncase,
          cuttrd = cuttrd,
          shortp = shortp,
          linkfunc = linkfunc,
          pD = 1,
          pGOF = 1
        )

        outasp <- hybdproj_get_projections(
          cdat,
          pdat,
          startp = startp,
          out,
          Ave5 = Ave5
        )
        outann <- standardize_annual_rates(outasp, pdat, stdpop = standpop)
        mod <- "a-s-nb"
      } else if (methods == "age-trd-poi") {
        out <- hybdproj(
          cdat,
          pdat,
          standpop,
          projfor = projfor,
          nagg = nagg,
          ncase = ncase,
          cuttrd = cuttrd,
          shortp = shortp,
          linkfunc = linkfunc,
          pD = 0,
          pGOF = 1
        )

        outasp <- hybdproj_get_projections(
          cdat,
          pdat,
          startp = startp,
          out,
          Ave5 = Ave5
        )
        outann <- standardize_annual_rates(outasp, pdat, stdpop = standpop)
        mod <- "a-s-poi"
      } else if (methods == "com-trd") {
        out <- hybdproj(
          cdat,
          pdat,
          standpop,
          projfor = projfor,
          nagg = nagg,
          ncase = ncase,
          cuttrd = cuttrd,
          shortp = shortp,
          linkfunc = linkfunc,
          pD = 1,
          pGOF = 0
        )

        outasp <- hybdproj_get_projections(
          cdat,
          pdat,
          startp = startp,
          out,
          Ave5 = Ave5
        )
        outann <- standardize_annual_rates(outasp, pdat, stdpop = standpop)
        mod <- "c-t"
      } else if (methods == "age-only") {
        out <- hybdproj(
          cdat,
          pdat,
          standpop,
          projfor = projfor,
          nagg = nagg,
          ncase = ncase,
          cuttrd = cuttrd,
          shortp = shortp,
          linkfunc = linkfunc,
          pD = 0,
          pGOF = 0
        )

        outasp <- hybdproj_get_projections(
          cdat,
          pdat,
          startp = startp,
          out,
          Ave5 = Ave5,
          sum5 = sum5
        )
        outann <- standardize_annual_rates(outasp, pdat, stdpop = standpop)
        mod <- "average"
      } else if (methods == "ave5") {
        out <- ave5proj(cdat, pdat, startp = startp, sum5 = sum5)
        outasp <- ave5proj_get_projections(pdat, out, standpop = NULL)
        outann <- ave5proj_get_projections(pdat, out, standpop = standpop)
        mod <- "average5"
      }
    }
  }

  obsy <- ncol(cdat)
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


#' canproj_get_projections
#'
#' Extract projection results.
#'
#' @param canproj.object An object based on the 'canproj()' function.
#' @inheritParams canproj
#'
#' @return A `data.frame()`
#'
#' @export
canproj_get_projections <- function(canproj.object, standpop = NULL) {
  if (is.null(standpop)) {
    return(canproj.object$agsproj)
  } else {
    return(canproj.object$annproj)
  }
}


#' summary.canproj
#'
#' Summarize information on projection method used.
#'
#' @param object An object based on the 'canproj()' function.
#' @param ... Other parameters.
#'
#' @return An information table describing the method used.
#'
#' @export
summary.canproj <- function(object, ...) {
  summary(object$out)
  invisible(object)
}


#' glm.canproj
#'
#' Summarize estimations from the final model.
#'
#' @inheritParams canproj_get_projections
#'
#' @return A summary table from the `glm.object`.
#'
#' @export
glm.canproj <- function(canproj.object) {
  summary(canproj.object$out$glm)
}

#' plot.canproj
#'
#' Create the graph of the observed and projected age-standardized rates from a
#' canproj object.
#'
#' @inheritParams canproj
#' @param x An object based on the 'canproj()' function.
#' @param startplot Start for usage of data. Default (`1`) uses all years,
#'  increasing this number removes the oldest years from the graph.
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param main Title for graph
#' @param labels Labels for age groups. `NULL` matches canproj labels.
#' @param ylim y-axis limit
#' @param lty Line type. Applies to oberved rates and predicted rates, respectively.
#' @param col Line colour.Applies to observed rates and predicted rates, respectively.
#' @param new Produce graph in new window (`T`), or in the same window (`F`).
#' @param ... Other parameters
#' @export
plot.canproj <- function(
  x,
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
  if (!inherits(x, "canproj")) {
    stop("Variable \"x\" must be of type \"canproj\"")
  }
  # Reading & formatting data:
  indat <- canproj_get_projections(x, standpop = standpop)
  indata <- indat[, 1]
  indata <- indata[startplot:length(indata)]

  # Setting internal variables:
  nopredy <- length(indata) - x$obsy
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
    graphics::plot(
      c(1, maxx),
      ylim,
      type = "n",
      ylab = ylab,
      xlab = xlab,
      axes = F,
      ...
    )
    graphics::axis(2)
    graphics::axis(1, at = 1:maxx, labels = labels)
    graphics::box()
    graphics::title(main)
  }
  graphics::lines(
    1:(maxx - nopredy),
    indata[1:(maxx - nopredy)],
    type = "o",
    pch = 20,
    lty = lty[1],
    col = col[1],
    ...
  )
  graphics::lines(
    (maxx - nopredy):maxx,
    indata[(maxx - nopredy):maxx],
    lty = lty[2],
    col = col[2],
    ...
  )
  # Returning object as invisible
  invisible(x)
}
