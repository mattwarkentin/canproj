#' Canproj: Cancer Projections
#'
#' Cancer incidence and mortality projections using the Canproj method. In
#'   brief, Canproj uses a decision-tree approach to determine the optimal model
#'   to project cancer data into the future based on historical trends. More
#'   information about this approach can be found here:
#'   <https://doi.org/10.24095/hpcdp.40.9.02>.
#'
#' @param cdat (age groups) * N (years) historical cancer data, 15<=N<=125.
#' @param pdat (age groups) * (N + M) (years) observed and projected population,
#'   5<=M<=25.
#' @param startp The start calendar year of projection (e.g., 2009).
#' @param standpop A `StandardPopulation` object that provides the weights
#'   (proportions) for each age groups in a standard population.
#' @param projfor Specify `"incidence"` or `"mortality"` if you want ASR as a
#'   criteria for `nagg`.
#' @param nagg Number of years for data aggregation (by years). Default is
#'   1-annual data.
#' @param ncase Minimum number of cancer cases/deaths per year for splitting
#'   data.
#' @param startage Youngest age group to include in the GLM. Default (`NULL`)
#'   picks based on mean cases.
#' @param newcohort Assign new cohort effect as `0` (`FALSE`) or the last
#'   estimated cohort effect (`TRUE`), default is `0`, use `TRUE` only if having
#'   evidence on negative new cohort effect.
#' @param Ave5 `Ave5 = TRUE` invokes the 5-year average method when age-only
#'   model is selected.
#' @param sum5 When the 5-year average method is used, `sum5 = TRUE` call the
#'   5-year period based rate, otherwise (`sum5 = FALSE`), average the 5 rates
#'   in the 5 years for each age group.
#' @param methods User required projection method can be specified by ADPC
#'   models: `"nordpred"` or `"adpc-nb"`; age-cohort models: `"ac"`, `"ac-nb"`;
#'   age-period models: `"age-trd"`, `"com-trd"`; and average: `"ave5"`.
#' @param linkfunc Link function. Default is `"power5"`. Can be one of `"log"`,
#'   `"sqrt"`, or `"identity"`.
#' @param cuttrd Degenerating percent of trends per year after 5 years
#'   (`shortp = 0`) or the first projection year.
#' @param shortp Attenuation percent of drift term or slope for the first 5
#'   projection years.
#' @param pD Trend selecting criteria of p-value of drift (linear trend) term.
#' @param pGOF Model selection criteria of p-value of goodness-of-fit.
#'
#' @param x A object of class `"canproj"` to print.
#' @param ... Not currently used.
#'
#' @md
#'
#' @return A named-`list` with class `"canproj"`. The `list` contains the
#'   following:
#'
#'   - `annproj`: A `matrix` of age-standardized rates and case counts.
#'   - `agsproj`: A `data.frame` of case counts for each age group and year.
#'   - `method`: The chosen method for projection. One of `"average"`,
#'     `"average5"`, `"ADPC"`, `"AC"`, `"Hybrid"`, `"nordpred"`, `"adpc-nb"`,
#'     `"ac-poi"`, `"ac-nb"`, `"a-s-nb"`, `"a-s-poi"`, or `"c-t"`.
#'   - `out`: The output from function call for the chosen `method`
#'     (e.g., `acproj()`, `adpcproj()`, `ave5proj()`, `hybdproj()`).
#'   - `obsy`: Number of observed years of cancer data from `cdat`.
#'   - `projy`: Number of projected years based on `pdat`.
#'   - `pdPC`: A vector of p-values for the drift, period, and cohort effects
#'     in the `"adpc"` model.
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
  newcohort = FALSE,
  Ave5 = FALSE,
  sum5 = TRUE,
  methods = NULL,
  linkfunc = "power5",
  cuttrd = 0.04,
  shortp = 0,
  pD = 0.05,
  pGOF = 0.05
) {
  S7::check_is_S7(standpop, StandardPopulation)

  if (ncol(cdat) < 15) {
    rlang::abort("\"cdat\" must have at least 15 years")
  }

  if (ncol(pdat) < 20) {
    rlang::abort("\"pdat\" must have at least 20 years")
  }

  if ((ncol(pdat) - ncol(cdat)) > 30) {
    rlang::abort("Package can not project more than 30 years")
  }

  if (!is.null(ncase) & !inherits(ncase, "numeric")) {
    rlang::abort("\"ncase\" must be of type \"numeric\" or NULL")
  }

  if (!inherits(Ave5, "logical")) {
    rlang::abort("\"Ave5\" must be of type \"logical\"")
  }

  if (!inherits(sum5, "logical")) {
    rlang::abort("\"sum5\" must be of type \"logical\"")
  }

  validate_projection_inputs(
    cdat,
    pdat,
    startp = startp,
    projfor = projfor,
    nagg = nagg,
    ncase = ncase,
    startage = startage,
    newcohort = newcohort,
    linkfunc = linkfunc,
    cuttrd = cuttrd,
    shortp = shortp,
    pD = pD,
    pGOF = pGOF
  )

  if (
    !is.null(methods) &&
      !(methods %in%
        list(
          "nordpred",
          "adpc-nb",
          "ac-poi",
          "ac-nb",
          "age-trd-poi",
          "age-trd-nb",
          "com-trd",
          "age-only",
          "ave5"
        ))
  ) {
    rlang::abort(
      "`methods` must be one of \"nordpred\", \"adpc-nb\", \"ac-poi\", \"ac-nb\", \"age-trd-poi\", \"age-trd-nb\", \"com-trd\", \"age-only\", \"ave5\", or NULL."
    )
  }

  if (mean(apply(cdat, 2, sum)[(ncol(cdat) - 9):(ncol(cdat))]) < 3) {
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

          r0 <- adpcproj_get_predictions(out, incidence = TRUE)
          outasp <- interpolate_age_specific_rates(
            r0,
            cdat,
            pdat,
            startp = startp,
            nagg = 5
          )
          outann <- standardize_annual_rates(
            outasp,
            pdat,
            startp,
            stdpop = standpop
          )
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

        outann <- standardize_annual_rates(
          outasp,
          pdat,
          startp,
          stdpop = standpop
        )
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

        r0 <- adpcproj_get_predictions(out, incidence = TRUE)
        outasp <- interpolate_age_specific_rates(
          r0,
          cdat,
          pdat,
          startp = startp,
          nagg = 5
        )
        outann <- standardize_annual_rates(
          outasp,
          pdat,
          startp,
          stdpop = standpop
        )
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

        r0 <- adpcproj_get_predictions(out, incidence = TRUE)
        outasp <- interpolate_age_specific_rates(
          r0,
          cdat,
          pdat,
          startp = startp,
          nagg = 5
        )
        outann <- standardize_annual_rates(
          outasp,
          pdat,
          startp,
          stdpop = standpop
        )
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
        outann <- standardize_annual_rates(
          outasp,
          pdat,
          startp,
          stdpop = standpop
        )
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
        outann <- standardize_annual_rates(
          outasp,
          pdat,
          startp,
          stdpop = standpop
        )
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
        outann <- standardize_annual_rates(
          outasp,
          pdat,
          startp,
          stdpop = standpop
        )
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
        outann <- standardize_annual_rates(
          outasp,
          pdat,
          startp,
          stdpop = standpop
        )
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
  projy <- ncol(pdat) - ncol(cdat)
  result <- list(
    annproj = outann,
    agsproj = outasp,
    method = mod,
    out = out,
    obsy = obsy,
    projy = projy,
    pdPC = pdPC
  )
  class(result) <- c("canproj", class(result))
  attr(result, "Call") <- sys.call()
  result
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
  if (!inherits(canproj.object, "canproj")) {
    rlang::abort("Variable \"canproj.object\" must be of type \"canproj\"")
  }

  if (!is.null(standpop)) {
    S7::check_is_S7(standpop, StandardPopulation)
  }

  if (is.null(standpop)) {
    return(canproj.object$agsproj)
  } else {
    return(canproj.object$annproj)
  }
}


#' Print method
#' @rdname canproj
#' @export
print.canproj <- function(x, ...) {
  cli::cat_line(glue::glue("Canproj ({x$method} approach)"))
  cli::cli_alert("Observed data: {x$obsy} years")
  cli::cli_alert("Projected data: {x$projy} years")
  invisible(x)
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
#' @param lty Line type. Applies to observed rates and predicted rates, respectively.
#' @param col Line colour.Applies to observed rates and predicted rates, respectively.
#' @param ... Other parameters
#'
#' @importFrom rlang .data
#'
#' @export
plot.canproj <- function(
  x,
  standpop,
  startplot = 1,
  xlab = "Calendar Year",
  ylab = "Rate per 100,000 people",
  main = "",
  lty = c(1, 2),
  col = c("black", "black"),
  ...
) {
  if (!inherits(x, "canproj")) {
    rlang::abort("Variable \"x\" must be of type \"canproj\"")
  }

  indat <- canproj_get_projections(x, standpop = standpop)
  indata <- indat[, 1]
  data <- as.data.frame(indata)
  data$year <- as.numeric(names(indata))

  observed <- data[startplot:x$obsy, ]
  projected <- data[(x$obsy + 1):length(indata), ]

  plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = observed,
      mapping = ggplot2::aes(x = .data$year, y = .data[[2]]),
      colour = col[1],
      size = 2
    ) +
    ggplot2::geom_line(
      data = observed,
      mapping = ggplot2::aes(x = .data$year, y = .data[[2]]),
      colour = col[1],
      linetype = lty[1],
      linewidth = 1
    ) +
    ggplot2::geom_line(
      data = projected,
      mapping = ggplot2::aes(x = .data$year, y = .data[[2]]),
      colour = col[2],
      linetype = lty[2],
      linewidth = 1
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::ggtitle(main) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black"),
      plot.title = ggplot2::element_text(size = 15)
    )

  return(plot)
}
