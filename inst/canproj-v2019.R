############################################################################## 
#    CANPROJ:  R functions for projection of cancer incidence/mortality      #
# For Cancer Projections Network (C-Proj)                                    #
# Sponsored by Canadian Partnership Against Cancer (CPAC)                    #
# Combining modified Nordpred, AC-model and Hybrid packages:                 #
# (i) Modified Nordpred: it is a age-period-cohort with drift term GLM       #
#     revised the `nordpred' by introducing negative binomial distribution   #
#     when lack of fit appears from Nordpred, additional link functions of   #
#     sqrt and identity, and settings of startestage and startuseage;        #
# (ii)Modified Hybrid: the final model can be selected from age-GLM,         #
#     common-trend GLM, age-specific GLM and negative binomial based GLM     # 
#     revised the original `hybrid' by adding cut-trend parameter            #
#     and power 5 link function                                              #
# (iii) Age-cohort model based projection method                             #
# Annual rates and numbers are produced                                      #
# Code modified to adapt to 19 age groups (separate 85+ into 85_89 and 90+)  #
# Written by: Dr. Zhenguo QIU, 2011                                          #
# Modified by: Yibing Ruan, 2019
# License:	GNU version I                                                    #
############################################################################## 

## call MASS library to get glm.nb
library(MASS)

## Control function of canproj:
canproj <- function(cdat, pdat, startp, projfor="incidence", 
             nagg=NULL, ncase=NULL, startestage=NULL, newcohort=NULL, Ave5=NULL, sum5=NULL, methods=NULL,   
             linkfunc="power5", cuttrd=0.04, shortp=0, pD=0.05, pGOF=0.05, standpop=ca11) {
#cdat:     19(age groups)*N(years) historical cancer data, 15<=N<=125
#pdat:     19(age groups)*(N+M)(years) observed and projected population, 5<=M<=25
#ncase:    minimum number of cancer cases/deaths per year for splitting data.
#nagg:     number of years for data aggregation (by years), default: 1-annual data
#startestage: user defined age groups for modeling can be input here
#cuttrd:   degenerating percent of trends per year after 5 years (shortp=0) 
#          or the first projection year.
#projfor:  specify "incidence" or "mortality" if want ASR as criteria for nagg
#newcohort:assign new cohort effect as 0 (NULL) or the last estimated cohort effect (T),
#          default is 0, use "T" only if having evidence on negative new cohort effect 
#Ave5:     Ave5=T invokes the 5 year average method when age-only model is selected
#sum5:     when the 5-year average method is used, sum5=NULL call the 5-year period based rate,
#          otherwise, average the 5 rates in the 5 years for each age group.
#methods:  user required projection method can be specified by ADPC models: "nordpred" or "adpc-nb"; 
#          age-cohort models: "ac", "ac-nb"; age-period models: "age-trd", "com-trd"; and average: "ave5"
#linkfunc: link function, default is power5, can be log, sqrt and identity
#pD:       trend selecting criteria of p-value of drift (linear trend) term
#pGOF:     model selection criteria of p-value of goodness-of-fit
#standpop: the weights (proportions) of 19 age groups in a standard population
#startp:   the start calendar year of projection, e.g. 2009
  
# Check data:
    if ( dim(cdat)[1]!=19 || dim(pdat)[1]!=19 ) {
      stop("\"cdat\" and \"pdat\" must have data for 19 age groups by row in this version")	
    }
    if ( dim(cdat)[2] > dim(pdat)[2]) {
      stop("\"pdat\" must include information about all years in \"cdat\"")	
    }
    if (dim(pdat)[2]==dim(cdat)[2]) {
      stop("\"pdat\" must include information in projection years")	
    }
    if ((dim(pdat)[2]-dim(cdat)[2]) > 30) {
      stop("Package can not project more than 30 years")	
    }
## The five year average method is used due to small annual number
  if (mean(apply(cdat, 2, sum)[(dim(cdat)[2] - 9):(dim(cdat)[2])]) < 3) {      	
      out <- ave5proj(cdat, pdat, startp=startp, sum5=sum5)
      outasp <- ave5proj.getproj(pdat, out, standpop=NULL)
      outann <- ave5proj.getproj(pdat, out, standpop=standpop)
      mod <- "average5"
      pdPC <- NA
  } else {

# Define number of cases for data splitting:
    if (is.null(ncase)) {
       if (projfor=="incidence") {
          ncase <- 1
       } else {
          ncase <- 3/5
       } 
    }

# Define number of years for data aggregation:
    if (is.null(nagg)) {
      masr <- mean(obasr(cdat, pdat)[,1])
      if (projfor=="incidence") {   
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

  nototper  <- dim(pyr)[2]
  noobsper  <- dim(cases)[2]
  nonewpred <- nototper - noobsper

# Preparing parameters for APC method:
    n5case <- 5*ncase
    cuttrend <- rep(shortp, nonewpred)
    for (i in 1:nonewpred) {
        cuttrend[i] <- shortp + (i-1)*5*cuttrd
    }
    cuttrend[cuttrend > 1] <- 1

# when projection method is not specified:
 if (is.null(methods)) {
# Performing adpc modeling:
      mod0 <- adpcproj(cdat, pdat, projfor=projfor, n5case=n5case,  
                       newcohort=newcohort, pGOF=0, cuttrd=0.04, shortp=0, linkfunc=linkfunc)
      pv1 <- summary(mod0$glm)$coef["Period", 4]
  options(warn=-1)
      pv2 <- drop1(mod0$glm, test="F")["as.factor(Cohort)", 5]
      pv3 <- drop1(mod0$glm, test="F")["as.factor(Period)", 5]
  options(warn=0)
      pdPC <- c(pv1, pv3, pv2)
# Selecting method by decision tree:
      if (pv2 < pGOF) {
         if (pv1 < pD | pv3 < pGOF) {
            out <- adpcproj(cdat, pdat, projfor=projfor, n5case=n5case, startestage=startestage, 
                            newcohort=newcohort, pGOF=pGOF, cuttrd=0.04, shortp=0, linkfunc=linkfunc)
            r0 <- adpcproj.getpred(out, incidence=T)
            outasp <- asrpy(r0, cdat, pdat, startp=startp, nagg=5)
            outann <- asry(outasp, pdat, standpop=standpop)
            mod <- "ADPC"
         } else {
            out <- acproj(cdat, pdat, projfor=projfor, n5case=n5case, startestage=startestage,
                            cuttrd=0.04, shortp=0, pGOF=pGOF, linkfunc=linkfunc)
            outasp <- acproj.getproj(cdat, pdat, startp=startp, out)
            outann <- acproj.getproj(cdat, pdat, startp=startp, out, standpop=standpop)
            mod <- "AC"
         }
      } else {
          out <- hybdproj(cdat, pdat, projfor=projfor, nagg=nagg, ncase=ncase, 
                          cuttrd=cuttrd, shortp=shortp, linkfunc=linkfunc, pD=pD, pGOF=pGOF)
          outasp <- hybdproj.getproj(cdat, pdat, startp=startp, out, Ave5=Ave5, sum5=sum5)
          outann <- asry(outasp, pdat, standpop=standpop)
          mod <- "Hybrid"
      }
 } else { 
# when specified method: "nordpred", "adpc-nb", "ac-poi", "ac-nb", 
#                        "age-trd-nb", "age-trd-poi", "com-trd", "age-only", "ave5":
      pdPC <- NA
      if (methods=="nordpred") {
            out <- adpcproj(cdat, pdat, projfor=projfor, n5case=n5case, startestage=startestage, 
                            newcohort=newcohort, pGOF=0, cuttrd=0.04, shortp=0, linkfunc=linkfunc )
            r0 <- adpcproj.getpred(out, incidence=T)
            outasp <- asrpy(r0, cdat, pdat, startp=startp, nagg=5)
            outann <- asry(outasp, pdat, standpop=standpop)
            mod <- "nordpred"
      } 
      if (methods=="adpc-nb") {
            out <- adpcproj(cdat, pdat, projfor=projfor, n5case=n5case, startestage=startestage, 
                            newcohort=newcohort, pGOF=1, cuttrd=0.04, shortp=0, linkfunc=linkfunc )
            r0 <- adpcproj.getpred(out, incidence=T)
            outasp <- asrpy(r0, cdat, pdat, startp=startp, nagg=5)
            outann <- asry(outasp, pdat, standpop=standpop)
            mod <- "adpc-nb"
      }
      if (methods=="ac-poi") {
            out <- acproj(cdat, pdat, projfor=projfor, n5case=n5case, startestage=startestage,
                            cuttrd=0.04, shortp=0, pGOF=0, linkfunc=linkfunc)
            outasp <- acproj.getproj(cdat, pdat, startp=startp, out)
            outann <- acproj.getproj(cdat, pdat, startp=startp, out, standpop=standpop)
            mod <- "ac-poi"
      }
      if (methods=="ac-nb") {
            out <- acproj(cdat, pdat, projfor=projfor, n5case=n5case, startestage=startestage,
                            cuttrd=0.04, shortp=0, pGOF=1, linkfunc=linkfunc)
            outasp <- acproj.getproj(cdat, pdat, startp=startp, out)
            outann <- acproj.getproj(cdat, pdat, startp=startp, out, standpop=standpop)
            mod <- "ac-nb"
      }
      if (methods=="age-trd-nb") {
            out <- hybdproj(cdat, pdat, projfor=projfor, nagg=nagg, ncase=ncase, 
                          cuttrd=cuttrd, shortp=shortp, linkfunc=linkfunc, pD=1, pGOF=1)
            outasp <- hybdproj.getproj(cdat, pdat, startp=startp, out, Ave5=Ave5)
            outann <- asry(outasp, pdat, standpop=standpop)
            mod <- "a-s-nb"
      }
      if (methods=="age-trd-poi") {
            out <- hybdproj(cdat, pdat, projfor=projfor, nagg=nagg, ncase=ncase, 
                          cuttrd=cuttrd, shortp=shortp, linkfunc=linkfunc, pD=0, pGOF=1)
            outasp <- hybdproj.getproj(cdat, pdat, startp=startp, out, Ave5=Ave5)
            outann <- asry(outasp, pdat, standpop=standpop)
            mod <- "a-s-poi"
      }
      if (methods=="com-trd") {
            out <- hybdproj(cdat, pdat, projfor=projfor, nagg=nagg, ncase=ncase, 
                          cuttrd=cuttrd, shortp=shortp, linkfunc=linkfunc, pD=1, pGOF=0)
            outasp <- hybdproj.getproj(cdat, pdat, startp=startp, out, Ave5=Ave5)
            outann <- asry(outasp, pdat, standpop=standpop)
            mod <- "c-t"
      }
      if (methods=="age-only") {
            out <- hybdproj(cdat, pdat, projfor=projfor, nagg=nagg, ncase=ncase, 
                          cuttrd=cuttrd, shortp=shortp, linkfunc=linkfunc, pD=0, pGOF=0)
            outasp <- hybdproj.getproj(cdat, pdat, startp=startp, out, Ave5=Ave5, sum5=sum5)
            outann <- asry(outasp, pdat, standpop=standpop)
            mod <- "average"
      }
      if (methods=="ave5") {
            out <- ave5proj(cdat, pdat, startp=startp, sum5=sum5)
            outasp <- ave5proj.getproj(pdat, out, standpop=NULL)
            outann <- ave5proj.getproj(pdat, out, standpop=standpop)
            mod <- "average5"
      }
  }
 }
      obsy <- dim(cdat)[2]
      resut <- list(annproj=outann, agsproj=outasp, method=mod, out=out, obsy=obsy, pdPC=pdPC)
      class(resut) <- "canproj"
      attr(resut,"Call") <- sys.call()
      return(resut)   
}

canproj.getproj <- function(canproj.object, standpop=NULL) {
  if (is.null(standpop)) {
    return(canproj.object$agsproj)
  } else {
    return(canproj.object$annproj)
  }
}

summary.canproj <- function(canproj.object) {
     summary(canproj.object$out)
     invisible(canproj.object)
}

glm.canproj <- function(canproj.object) {
     summary(canproj.object$out$glm)
}

plot.canproj <- function(canproj.object, standpop, startplot=1, xlab="Calendar Year",ylab="Rates",
                main="",labels=NULL,ylim=NULL,lty=c(1,3),col=c(1,1),new=T,...) {
  if (class(canproj.object)!="canproj") {
    stop("Variable \"canproj.object\" must be of type \"canproj\"")	
  } 
  # Reading & formating data:
  indat <- canproj.getproj(canproj.object, standpop=standpop)
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
    maxy <- max(indata)*(21/20)	
    if (is.null(ylim)) {
      ylim <- c(0,maxy)	
    }
    plot(c(1,maxx),ylim,type="n",ylab=ylab,xlab=xlab,axes=F,...)
    axis(2)
    axis(1,at=1:maxx,labels=labels)
    box()
    title(main)
  }
  lines(1:(maxx-nopredy),indata[1:(maxx-nopredy)],type="o",pch=20,lty=lty[1],col=col[1],...)
  lines((maxx-nopredy):maxx,indata[(maxx-nopredy):maxx],lty=lty[2],col=col[2],...)
  # Returning object as invisible
  invisible(canproj.object)
}



############################################################################## 
#    ADPCPROJ:  R functions for projection of cancer incidence/mortality     #
# Revising nordpred and introducing negative binomial distribution           #
# when lack of fit appears from nordpred, additional link functions of       #
# sqrt and identity, and settings of startestage and startuseage             #
# Written by: Dr. Zhenguo QIU, 2011                                          #
# License:	GNU version 1                                                    #
############################################################################## 

adpcproj <- function(cdat, pdat, projfor="incidence", n5case=NULL, noperiods=NULL, recent=NULL, 
              startestage=NULL, newcohort=NULL, pGOF=0.05, cuttrd=0.04, shortp=0, linkfunc="power5") {  
# Define number of cases for data splitting:
    if (is.null(n5case)) {
       if (projfor=="incidence") {
          n5case <- 5
       } else {
          n5case <- 3
       } 
    }

## aggregating data by 5 years:
   aggdata <- datagg(cdat, pdat, 5)
   cases <- aggdata$cases
   pyr <- aggdata$pyr
## Setting startestage and startuseage:
   if (is.null(startestage)) {
       dat <- as.matrix(cases)
       iage <- 1
       while(mean(dat[iage, ]) < n5case | dat[iage, 1] == 0) {
            iage <- iage + 1
       }
       startestage <- iage   
   } 
   if (is.null(newcohort)) {
      startuseage <- startestage
   } else {
      startuseage <- startestage + 1 # assign newcohort effects as the last estimated cohort effect
   }                          

## Setting default and checking data:
  # Number of periods in observed data
  percases <- dim(cases)[2]   
  if (percases < 3) { 
     stop("Minimum number of period is 3 (15 years) in \"cases\"") 
  }
   
## Setting number of periods for projection base
  # List possible candidates for number of periods to base predictions on
  # Default is 4:6 if available
  if (is.null(noperiods)) { 
     noperiods <- c(min(percases,4):min(percases,6)) 
  }
 
  # Choose number of periods by cutting stepwise execution of the
  # highest candidate number (i.e. cutting the most ancient periods)
  noperiods <- sort(noperiods)
  while(length(noperiods) > 1) {
    maxnoperiod <- max(noperiods)
    glmn <- adpcproj.estimate(cases, pyr, maxnoperiod, startestage)$glm
    pval <- 1 - pchisq(glmn$deviance, glmn$df.residual)
    if (pval < 0.01) { 
        noperiods <- noperiods[1:(length(noperiods) - 1)]
    } else {
        noperiods <- maxnoperiod
    }
  }      
  noperiod <- noperiods
  
## Setting status for recent (whether to use recent trend or whole trend)
  if (is.null(recent)) {
    recent <- adpcproj.estimate(cases, pyr, noperiod, startestage)$suggestionrecent
  }
  
## Perform estimation and prediction:
  est  <- adpcproj.estimate(cases=cases, pyr=pyr, noperiod=noperiod, pGOF=pGOF,
                            startestage=startestage, linkfunc=linkfunc)
  pred <- adpcproj.prediction(adpcproj.estimate.object=est, startuseage=startuseage,
                              recent=recent, cuttrd=0.04, shortp=0)
  return(pred)
}

## Fitting models:
adpcproj.estimate <- function(cases, pyr, noperiod, startestage, pGOF=0.05, linkfunc="power5") {
## Checking data
  if ( dim(cases)[1]!=19 || dim(pyr)[1]!=19 ) {
    stop("\"cases\" and \"pyr\" must have data for 19 age groups")	
  }

  if ( dim(cases)[2]>dim(pyr)[2]) {
    stop("\"pyr\" must include information about all periods in \"cases\"")	
  }
  
  if (dim(pyr)[2]==dim(cases)[2]) {
    stop("\"pyr\" must include information on future rates")	
  }

  if ((dim(pyr)[2]-dim(cases)[2])>6) {
    stop("Package can not project more than 6 periods")	
  }

  if ((dim(cases)[2]-noperiod)<0) {
    stop("More periods specified in \"noperiod\" than columns in \"cases\"")   	
  }

  if (noperiod<3) {
    stop("\"noperiod\" must be 3 or larger")   	
  }  

## Setting internal variables:
  dnoperiods <- dim(cases)[2]
  dnoagegr   <- dim(cases)[1]

## Transform dataformat:
  ageno    <- rep(1:dnoagegr,dnoperiods)
  periodno <- sort(rep(1:dnoperiods,dnoagegr))
  cohort   <- max(ageno)-ageno+periodno
  y        <- c(as.matrix(pyr[,1:dnoperiods]))
  apcdata <- data.frame(Age=ageno, Cohort=cohort, Period=periodno,
                        Cases=c(as.matrix(cases)), y=y)

## Selecting data for regression:
  apcdata <- apcdata[apcdata$Age >= startestage, ]
  apcdata <- apcdata[apcdata$Period > (dnoperiods - noperiod), ]
  

## Setting contrast for age+drift+period+cohort
options(contrasts=c("contr.treatment","contr.poly")) 

## Estimation:
  if (linkfunc=="power5") {
## Creation of power5 link:
 # Setting population variable 
    y <- apcdata$y
 # Make power5 link function for poisson family:
    power5link <- poisson()
    power5link$link <- "0.2 root link Poisson family"
    power5link$linkfun <- function(mu)  { (mu/y)^0.2 }
    power5link$linkinv <- function(eta) { pmax(.Machine$double.eps, y*eta^5) }
    power5link$mu.eta <- function(eta)  { pmax(.Machine$double.eps, 5*y*eta^4) }
    res.glm <- glm(Cases~as.factor(Age)+Period+as.factor(Period)+as.factor(Cohort) - 1,
                   data=apcdata, family=power5link)
  } else  if (linkfunc=="log") {
    res.glm <- glm(Cases~as.factor(Age)+Period+as.factor(Period)+as.factor(Cohort)+ offset(log(y))-1,
                   data=apcdata, family=poisson(link=log))  	
  } else  if (linkfunc=="sqrt") {
    res.glm <- glm(Cases/y ~ as.factor(Age)+Period+as.factor(Period)+as.factor(Cohort) -1,
                   data=apcdata, family=poisson(link=sqrt))  
  } else  if (linkfunc=="identity") {
    res.glm <- glm(Cases/y ~ as.factor(Age)+Period+as.factor(Period)+as.factor(Cohort) -1,
                   data=apcdata, family=poisson(link=identity))  
  } else {
    stop("Unknown \"linkfunc\"")	
  }
  pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df.residual)
  distn <- "Poisson"

## change model to negative binomial glm when lack of fit:
  if (pvalue < pGOF) {
    distn <- "Negative-binomial"
options(warn=-1)
  # Estimation:
    if (linkfunc=="power5") {
      y <- apcdata$y
    # Make power5 link function for negative binomial family:
      glmnb <- glm.nb(Cases ~as.factor(Age)+Period+as.factor(Period)+as.factor(Cohort) - 1
                      + offset(log(y)), data=apcdata, link=log)
      theta <- as.numeric(theta.md(apcdata$Cases, fitted(glmnb), dfr = df.residual(glmnb)))
      nbpower5link <- negative.binomial(theta)
      nbpower5link$link <- "0.2 root link negative.binomial(theta) family"
      nbpower5link$linkfun <- function(mu)  { (mu/y)^0.2 }
      nbpower5link$linkinv <- function(eta) { pmax(.Machine$double.eps, y*eta^5) }
      nbpower5link$mu.eta <- function(eta)  { pmax(.Machine$double.eps, 5*y*eta^4) }
      res.glm <- glm(Cases ~ as.factor(Age)+Period+as.factor(Period)+as.factor(Cohort) - 1,
                                   data=apcdata, family=nbpower5link)
    } else  if (linkfunc=="log") {
      res.glm <- glm.nb(Cases ~as.factor(Age)+Period+as.factor(Period)+as.factor(Cohort) - 1
                 + offset(log(y)), data=apcdata, link=log)  	
    } else  if (linkfunc=="sqrt") {
      res.glm <- glm.nb(Cases/y ~ as.factor(Age)+Period+as.factor(Period)+as.factor(Cohort) -1,
                                   data=apcdata, link=sqrt)  
    } else  if (linkfunc=="identity") {
      res.glm <- glm.nb(Cases/y ~ as.factor(Age)+Period+as.factor(Period)+as.factor(Cohort) -1,
                                   data=apcdata, link=identity)  
    } else {
      stop("Unknown \"linkfunc\"")	
    }
  # updating p-value of goodness of fit and distribution function:
    pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df.residual)
options(warn=0)
  }

## setting suggestion for 'recent' (whether to use recent trend or whole trend)
  if (distn=="poisson") {
      mod1 <- glm(Cases ~ as.factor(Age) + Period + as.factor(Cohort) + offset(log(y)) -1,
                  data=apcdata, family=poisson)
      mod2 <- glm(Cases ~ as.factor(Age)+Period +I(Period^2)+as.factor(Cohort) + offset(log(y)) -1,
                  data=apcdata, family=poisson)
      pdiff <- 1 - pchisq((mod1$deviance-mod2$deviance), (mod1$df.residual-mod2$df.residual))

      if (pdiff < 0.05) {
          suggestionrecent <- T
      } else {
          suggestionrecent <- F
      }   
  } else {
options(warn=-1)
      mod1 <- glm.nb(Cases ~ as.factor(Age) + Period + as.factor(Cohort) + offset(log(y)) -1,
                     data=apcdata, link=log)
      mod2 <- glm.nb(Cases ~ as.factor(Age)+Period +I(Period^2)+as.factor(Cohort) + offset(log(y)) -1,
                     data=apcdata, link=log)
options(warn=0)
      pdiff <- 1 - pchisq((mod2$twologlik-mod1$twologlik), (mod1$df.residual-mod2$df.residual))
      if (pdiff < 0.05) {
          suggestionrecent <- T
      } else {
          suggestionrecent <- F
      }   
  }

#options(contrasts("unordered", "ordered"))  

## Set class and return results
  res <- list(glm=res.glm, 
              cases=cases, pyr=pyr, 
              noperiod=noperiod, linkfunc=linkfunc, startestage=startestage, 
              distribution=distn, gofpvalue=pvalue,
              suggestionrecent=suggestionrecent, pvaluerecent=pdiff)
  class(res) <- "adpcproj.estimate"
  attr(res, "Call") <- sys.call()
  return(res)
}


## Project age-specific rates:
adpcproj.prediction <- function(adpcproj.estimate.object,
                                startuseage, recent, shortp=0, cuttrd=0.04) {
## running conditions:  
  if (class(adpcproj.estimate.object)!="adpcproj.estimate") {
    stop("Variable \"adpcproj.estimate.object\" must be of type \"adpcproj.estimate\"")	
  } 

  if (adpcproj.estimate.object$startestage > startuseage) {
    stop("\"startuseage\" is set too high compared to \"startestage\"") 	
  }

## Setting local variables:
  cases     <-  adpcproj.estimate.object$cases
  pyr	    <-  adpcproj.estimate.object$pyr
  noperiod  <- adpcproj.estimate.object$noperiod
  nototper  <- dim(pyr)[2]
  noobsper  <- dim(cases)[2]
  nonewpred <- nototper - noobsper

  cuttrend <- rep(shortp, nonewpred)
  for (i in 1:nonewpred) {
        cuttrend[i] <- shortp + (i-1)*5*cuttrd
  }
  cuttrend[cuttrend > 1] <- 1
  if (length(cuttrend) > nonewpred) {
    cuttrend <- cuttrend[1:nonewpred] 
  }
  if (length(cuttrend)<(dim(adpcproj.estimate.object$pyr)[2]
                       -dim(adpcproj.estimate.object$cases)[2])){
    err <- paste("\"cuttrend\" must always be at least the same length as")
    err <- paste(err,"the number of periods with population forecasts")
    stop(err)
  }
  
  if (is.data.frame(pyr)) {
    years <- names(pyr)
  } else {
    if (is.null(dimnames(pyr))) {
      years <- paste("Periode",1:nototper)	
    } else {
      years <- dimnames(pyr)[[2]]
    }	
  }  

## Making data object:
  datatable <- matrix(NA, 19, nototper)
  # fill in observed cases:
  datatable[ , 1:(nototper-nonewpred)] <- as.matrix(cases)
  datatable <- data.frame(datatable)
  row.names(datatable) <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
                            "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89", "90+")
  names(datatable) <- years
 
## Calculate predictions in number of cases:
  # For young agegroups with little data, use average from last two periods:
  for (age in 1:(startuseage-1)) {
    obsinc <- cases[age,(noobsper-1):noobsper]/pyr[age,(noobsper-1):noobsper]
    if (sum(is.na(obsinc))) {
      obsinc[is.na(obsinc)] <- 0 	
    }
    datatable[age,(noobsper+1):nototper] <- ((obsinc[,1]+obsinc[,2])/2)*pyr[age,(noobsper+1):nototper]
  }

  # For old agegroups, use trend from model estimates:
  for (age in startuseage:19) {
    startestage  <- adpcproj.estimate.object$startestage
    coefficients <- adpcproj.estimate.object$glm$coefficients

  # Cohort index: No. agegoups - age + period
    coh       <- (19-startestage) - (age-startestage) + (noperiod+1:nonewpred)

    noages    <- 19-startestage+1
    driftmp   <- cumsum(1-cuttrend)
    cohfind   <- noages + (noperiod-1) + 1 + (coh-1)
    maxcoh    <- 19 - startuseage + noperiod 
    agepar    <- as.numeric(coefficients[age-startestage+1])
    driftfind <- pmatch("Period", attributes(coefficients)$names)
    driftpar  <- as.numeric(coefficients[driftfind])
    
    cohpar <- rep(NA, length(coh))
    for (i in 1:length(coh)) {
      if (coh[i] < maxcoh) {
        cohpar[i] <- as.numeric(coefficients[cohfind[i]])
      } else {
        # For new young cohorts in the future, define the effect as the last cohort:
        cohpar[i] <- as.numeric(coefficients[length(coefficients)-(startuseage-startestage)])
#       cohpar[i][is.na(cohpar[i])] <- as.numeric(coefficients[length(coefficients)-1])
        cohpar[i][is.na(cohpar[i])] <- 0
      }
    } 
  
  # Getting recent drift (D_last) estimate:  
    if (recent) {
      # localize the last period effect (note: p.first=p.last)
      lpfind <- driftfind + noperiod-2
      # find the estimate of the last two period effect (P-1)
      lppar  <- as.numeric(coefficients[lpfind])
      driftrecent <- driftpar - lppar
    }

  # Project age-specific rates:    
    if (adpcproj.estimate.object$linkfunc=="power5") { 
      if (recent) {
        rate <- (agepar+driftpar*noobsper+driftrecent*driftmp+cohpar)^5      
      } else {
        rate <- (agepar+driftpar*(noobsper+driftmp)+cohpar)^5
      }
    } else  if (adpcproj.estimate.object$linkfunc=="log") {
      if (recent) {
        rate <- exp(agepar+driftpar*noobsper+driftrecent*driftmp+cohpar)      
      } else {
        rate <- exp(agepar+driftpar*(noobsper+driftmp)+cohpar)
      }
    } else  if (adpcproj.estimate.object$linkfunc=="sqrt") {
      if (recent) {
        rate <- (agepar+driftpar*noobsper+driftrecent*driftmp+cohpar)^2      
      } else {
        rate <- (agepar+driftpar*(noobsper+driftmp)+cohpar)^2
      }
    } else { # identity link:
      if (recent) {
        rate <- (agepar+driftpar*noobsper+driftrecent*driftmp+cohpar)      
      } else {
        rate <- (agepar+driftpar*(noobsper+driftmp)+cohpar)
      }
    }
    datatable[age,(noobsper+1):nototper] <- rate*pyr[age,(noobsper+1):nototper]
  }

  # Structure and return results:
  res <- list(predictions=datatable,
              pyr=pyr, nopred=nonewpred, noperiod=adpcproj.estimate.object$noperiod,
              gofpvalue=adpcproj.estimate.object$gofpvalue, recent=recent,
              pvaluerecent=adpcproj.estimate.object$pvaluerecent, cuttrd=cuttrd, shortp=shortp,
              cuttrend=cuttrend, distribution=adpcproj.estimate.object$distribution,
              startuseage=startuseage, startestage=adpcproj.estimate.object$startestage,
              glm=adpcproj.estimate.object$glm)
  class(res) <- "adpcproj"
  attr(res,"Call") <- sys.call()
  return(res)
}


## Summary Projection Results by 5-year period(ASR and total number):
adpcproj.getpred <- function(adpcproj.object,
                             incidence=T, standpop=NULL, excludeobs=F, byage, agegroups="all") {  
## Setting defaults:
  if (missing(byage)) {
     byage <- ifelse(is.null(standpop), T, F)
  }
  
## Checking input:
  if (class(adpcproj.object)!="adpcproj") {
     stop("Variable \"adpcproj.object\" must be of type \"adpcproj\"")	
  } 
  
  if ((!is.null(standpop)) && (!incidence)) {
     stop("\"standpop\" should only be used with incidence predictions (incidence=T)") 	
  }
 
  if (!is.null(standpop)) {
    if (round(sum(standpop),5)!=1) {
       stop("\"standpop\" must be of sum 1") 	
    }
    if ((length(standpop)!=length(agegroups)) && (agegroups[1]!="all")) {
       stop("\"standpop\" must be the same length as \"agegroups\"") 	
    }
    if (byage) {
       stop("\"standpop\" is only valid for \"byage=T\"") 
    }
  }
 
## Seting local data:
  datatable <- adpcproj.object$predictions
  pyr       <- data.frame(adpcproj.object$pyr)

## Secting agegroups:
  if (agegroups[1]!="all") {
     datatable <- datatable[agegroups, ]
     pyr       <- pyr[agegroups, ]
  }

## If needed; Standardize data and Collapse agegroups
  if (!is.null(standpop)) {
     datainc <- (datatable/pyr)*100000
     if (sum(is.na(datainc)) > 0) {
        datainc[is.na(datainc)] <- 0
     } 	
    res <- apply(datainc*standpop, 2, sum) 
  } else {
    if (!byage) {
       datatable <- apply(datatable, 2, sum)
       pyr       <- apply(pyr, 2, sum) 
    }
    if (incidence) {
       res <- (datatable/pyr)*100000
       if (sum(is.na(res)) > 0) {
          res[is.na(res)] <- 0
       } 	
    }  else {
          res <- datatable
    } 
  }

## Select data:
  if (excludeobs) {
     if (is.matrix(res)) {
        predstart <- dim(res)[2] - adpcproj.object$nopred + 1
        res <- res[ , predstart:(predstart+adpcproj.object$nopred-1)]
     }  else {
          predstart <- length(res) - adpcproj.object$nopred + 1
          res <- res[predstart:(predstart+adpcproj.object$nopred-1)]
     }
  } 
  
## Return data: 
  return(res)
}


## Summary Annual Projection Results:
adpcproj.getproj <- function(cdat, pdat, startp, adpcproj.object, standpop=NULL) {
  r0 <- adpcproj.getpred(adpcproj.object, incidence=T)
  outasp <- asrpy(r0, cdat, pdat, startp=startp, nagg=5)
  if (is.null(standpop)) {
    return(outasp)
  } else {
    outann <- asry(outasp, pdat, standpop=standpop)
    return(outann)
  }
}


summary.adpcproj <- function(adpcproj.object, printpred=F, printcall=F, digits=0) {
  method <- "Age-drift-Period-Cohort Model"
  if (class(adpcproj.object)!="adpcproj") {
    stop("Variable \"adpcproj.object\" must be of type \"adpcproj\"")	
  } 
  # Setting internal variables:
  obsto <- names(adpcproj.object$predictions)[dim(adpcproj.object$predictions)[2]-adpcproj.object$nopred]

  if (!is.null(adpcproj.object$pvaluerecent)) {
    precent <- round(adpcproj.object$pvaluerecent,4)
  } else { precent <- NA }
   
  if (!is.null(adpcproj.object$gofpvalue)) {
    gofpvalue <- round(adpcproj.object$gofpvalue,4)
  } else { gofpvalue <- NA }
      
  # Print information about object: 
  if (printpred) {
    cat("Observed and predicted values:")
    cat("(observations up to",obsto,")\n")
    print(round(as.matrix(adpcproj.object$predictions),digits=digits))
    cat("\n")
  }
  cat("\nPrediction done with:\n")
    
  moptions <- matrix(NA,10,2)
  moptions[,1] <- c("Method:",
                    "Number of periods predicted (nopred):",
                    "Trend used in predictions (cuttrend):",
                    "Number of periods used in estimate (noperiod):",
                    "Distribution function of regression:",
                    "P-value for goodness of fit:",
                    "Used recent (recent):", "P-value for recent:",
                    "First age group used (startuseage):","First age group estimated (startestage):")        
  moptions[,2] <- c(method,
                    adpcproj.object$nopred, paste(adpcproj.object$cuttrend, collapse=" , "),
                    adpcproj.object$noperiod, 
                    adpcproj.object$distribution,
                    gofpvalue,
                    adpcproj.object$recent, precent,
                    adpcproj.object$startuseage, adpcproj.object$startestage)
  maxl <- max(nchar(moptions[,1]))
    
  for (i in 1:dim(moptions)[1]) {
    spaces <- rep (" ", maxl-nchar(moptions[i,1])+2)
    cat(moptions[i,1], spaces, moptions[i,2], "\n", sep="")	
  }
    
  if (printcall) {
    cat("\n  Call: ")
    dput(attr(adpcproj.object, "Call"))
  }
  invisible(adpcproj.object)
}

glm.adpcproj <- function(adpcproj.object) {
     summary(adpcproj.object$glm)
}

plot.adpcproj<- function(cdat, pdat, startp, adpcproj.object, standpop, startplot=1, 
                         xlab="Calendar Year", ylab="Rates",
                main="",labels=NULL,ylim=NULL,lty=c(1,3),col=c(1,1),new=T,...) {
  if (class(adpcproj.object)!="adpcproj") {
    stop("Variable \"adpcproj.object\" must be of type \"adpcproj\"")	
  } 
  
  # Reading & formating data:
  indat <- adpcproj.getproj(cdat, pdat, startp=startp, adpcproj.object, standpop=standpop)
  indata <- indat[, 1]
  indata <- indata[startplot:length(indata)]

  # Seting internal variables:
  obsy <- dim(cdat)[2]
  nopredy <- length(indata) - obsy 
  if (is.null(labels)) {
    labels <- row.names(indat)
  }

  # Create plots:
  maxx <- length(indata)
  if (new) { 	
    maxy <- max(indata)*(21/20)	
    if (is.null(ylim)) {
      ylim <- c(0,maxy)	
    }
    plot(c(1,maxx),ylim,type="n",ylab=ylab,xlab=xlab,axes=F,...)
    axis(2)
    axis(1,at=1:maxx,labels=labels)
    box()
    title(main)
  }
  lines(1:(maxx-nopredy),indata[1:(maxx-nopredy)],type="o",pch=20,lty=lty[1],col=col[1],...)
  lines((maxx-nopredy):maxx,indata[(maxx-nopredy):maxx],lty=lty[2],col=col[2],...)

  # Returning object as invisible
  invisible(adpcproj.object)
}



############################################################################## 
#      ACPROJ:  R functions for projection of cancer incidence/mortality     #
# Revising and combining nordpred and Osmond's to extrapolation cohort       #
# when no drift appears from nordpred.                                       #
# Written by: Dr. Zhenguo QIU, 2011                                          #
# License:	GNU version 1                                                    #
############################################################################## 

acproj <- function(cdat, pdat, projfor="incidence", n5case=NULL, noperiods=NULL, startestage=NULL,
                   cuttrd=0.04, shortp=0, pGOF=0.05, linkfunc="power5") {  
# Define number of cases for data splitting:
    if (is.null(n5case)) {
       if (projfor=="incidence") {
          n5case <- 5
       } else {
          n5case <- 3
       } 
    }

## aggregating data by 5 years:
   aggdata <- datagg(cdat, pdat, 5)
   cases <- aggdata$cases
   pyr <- aggdata$pyr
## Setting startestage:
   if (is.null(startestage)) {
       dat <- as.matrix(cases)
       iage <- 1
       while(mean(dat[iage, ]) < n5case | dat[iage, 1] == 0) {
            iage <- iage + 1
       }
       startestage <- iage    
   }
## Setting default and checking data:
  # Number of periods in observed data
  percases <- dim(cases)[2]   
  if (percases < 3) { 
     stop("Minimum number of period is 3 (15 years) in \"cases\"") 
  }
   
## Setting number of periods for projection base
  noperiod <- percases
    
## Perform estimation and prediction:
  est  <- acproj.estimate(cases=cases, pyr=pyr, noperiod=noperiod, pGOF=pGOF,
                            startestage=startestage, linkfunc=linkfunc)
  pred <- acproj.prediction(acproj.estimate.object=est, cuttrd=0.04, shortp=0)
  return(pred)
}

## Fitting models:
acproj.estimate <- function(cases, pyr, noperiod, startestage, pGOF=0.05, linkfunc="power5") {
## Checking data
  if ( dim(cases)[1]!=19 || dim(pyr)[1]!=19 ) {
    stop("\"cases\" and \"pyr\" must have data for 19 age groups")	
  }

  if ( dim(cases)[2]>dim(pyr)[2]) {
    stop("\"pyr\" must include information about all periods in \"cases\"")	
  }
  
  if (dim(pyr)[2]==dim(cases)[2]) {
    stop("\"pyr\" must include information on future rates")	
  }

  if ((dim(pyr)[2]-dim(cases)[2])>6) {
    stop("Package can not project more than 6 periods")	
  }

  if ((dim(cases)[2]-noperiod)<0) {
    stop("More periods specified in \"noperiod\" than columns in \"cases\"")   	
  }

  if (noperiod<3) {
    stop("\"noperiod\" must be 3 or larger")   	
  }  

## Setting internal variables:
  dnoperiods <- dim(cases)[2]
  dnoagegr   <- dim(cases)[1]

## Transform dataformat:
  ageno    <- rep(1:dnoagegr,dnoperiods)
  periodno <- sort(rep(1:dnoperiods,dnoagegr))
  cohort   <- max(ageno)-ageno+periodno
  y        <- c(as.matrix(pyr[,1:dnoperiods]))
  apcdata <- data.frame(Age=ageno, Cohort=cohort, Period=periodno,
                        Cases=c(as.matrix(cases)), y=y)

## Selecting data for regression:
  apcdata <- apcdata[apcdata$Age >= startestage, ]
  apcdata <- apcdata[apcdata$Period > (dnoperiods - noperiod), ]
  maxc <- max(apcdata$Cohort)
  midc <- ceiling(maxc/2)
  
## Estimation:
  if (linkfunc=="power5") {
## Creation of power5 link:
 # Setting population variable 
    y <- apcdata$y
 # Make power5 link function for poisson family:
    power5link <- poisson()
    power5link$link <- "0.2 root link Poisson family"
    power5link$linkfun <- function(mu)  { (mu/y)^0.2 }
    power5link$linkinv <- function(eta) { pmax(.Machine$double.eps, y*eta^5) }
    power5link$mu.eta <- function(eta)  { pmax(.Machine$double.eps, 5*y*eta^4) }
    res.glm <- glm(Cases ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
                           data=apcdata, family=power5link)
  } else  if (linkfunc=="log") {
    res.glm <- glm(Cases ~ factor(Age) + relevel(factor(Cohort), midc)+ offset(log(y))-1,
                           data=apcdata, family=poisson(link=log))  	
  } else  if (linkfunc=="sqrt") {
    res.glm <- glm(Cases/y ~ factor(Age) + relevel(factor(Cohort), midc) -1,
                             data=apcdata, family=poisson(link=sqrt))  
  } else  if (linkfunc=="identity") {
    res.glm <- glm(Cases/y ~ factor(Age) + relevel(factor(Cohort), midc) -1,
                             data=apcdata, family=poisson(link=identity))  
  } else {
    stop("Unknown \"linkfunc\"")	
  }
  pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df.residual)
  distn <- "Poisson"

## change model to negative binomial glm when lack of fit:
  if (pvalue < pGOF) {
    distn <- "Negative-binomial"
  options(warn=-1)
  # Estimation:
    if (linkfunc=="power5") {
      y <- apcdata$y
    # Make power5 link function for negative binomial family:
      glmnb <- glm.nb(Cases ~ factor(Age) + relevel(factor(Cohort), midc) + offset(log(y))- 1, 
                              data=apcdata, link=log)
      theta <- as.numeric(theta.md(apcdata$Cases, fitted(glmnb), dfr = df.residual(glmnb)))
      nbpower5link <- negative.binomial(theta)
      nbpower5link$link <- "0.2 root link negative.binomial(theta) family"
      nbpower5link$linkfun <- function(mu)  { (mu/y)^0.2 }
      nbpower5link$linkinv <- function(eta) { pmax(.Machine$double.eps, y*eta^5) }
      nbpower5link$mu.eta <- function(eta)  { pmax(.Machine$double.eps, 5*y*eta^4) }
      res.glm <- glm(Cases ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
                             data=apcdata, family=nbpower5link)
    } else  if (linkfunc=="log") {
      res.glm <- glm.nb(Cases ~ factor(Age) + relevel(factor(Cohort), midc) + offset(log(y))- 1, 
                                data=apcdata, link=log)  	
    } else  if (linkfunc=="sqrt") {
      res.glm <- glm.nb(Cases/y ~ factor(Age) + relevel(factor(Cohort), midc) -1,
                                  data=apcdata, link=sqrt)  
    } else  if (linkfunc=="identity") {
      res.glm <- glm.nb(Cases/y ~ factor(Age) + relevel(factor(Cohort), midc) -1,
                                  data=apcdata, link=identity)  
    } else {
      stop("Unknown \"linkfunc\"")	
    }
  # updating p-value of goodness of fit and distribution function:
    pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df.residual)
  options(warn=0)
  }
  
## Set class and return results
  res <- list(glm=res.glm, 
              cases=cases, pyr=pyr, maxc=maxc, midc=midc,
              noperiod=noperiod, linkfunc=linkfunc, startestage=startestage, 
              distribution=distn, gofpvalue=pvalue)
  class(res) <- "acproj.estimate"
  attr(res, "Call") <- sys.call()
  return(res)
}

## Project age-specific rates:
acproj.prediction <- function(acproj.estimate.object, cuttrd=0.04, shortp=0 ) {
## running conditions:  
  if (class(acproj.estimate.object)!="acproj.estimate") {
    stop("Variable \"acproj.estimate.object\" must be of type \"acproj.estimate\"")	
  } 

## Setting local variables:
  cases     <-  acproj.estimate.object$cases
  pyr	    <-  acproj.estimate.object$pyr
  noperiod  <- acproj.estimate.object$noperiod
  startestage <- acproj.estimate.object$startestage
  maxc <- acproj.estimate.object$maxc
  midc <- acproj.estimate.object$midc
  nototper  <- dim(pyr)[2]
  noobsper  <- dim(cases)[2]
  nonewpred <- nototper - noobsper

  cuttrend <- rep(shortp, nonewpred)
  for (i in 1:nonewpred) {
        cuttrend[i] <- shortp + (i-1)*5*cuttrd
  }
  cuttrend[cuttrend > 1] <- 1
  if (length(cuttrend) > nonewpred) {
    cuttrend <- cuttrend[1:nonewpred] 
  }

  if (is.data.frame(pyr)) {
    years <- names(pyr)
  } else {
    if (is.null(dimnames(pyr))) {
      years <- paste("Periode",1:nototper)	
    } else {
      years <- dimnames(pyr)[[2]]
    }	
  }  

## Making data object:
  datatable <- matrix(NA, 19, nototper)
  # fill in observed cases:
  datatable[ , 1:(nototper-nonewpred)] <- as.matrix(cases)
  datatable <- data.frame(datatable)
  row.names(datatable) <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
                            "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+")
  names(datatable) <- years
 
## Calculate predictions in number of cases:
  # For young agegroups with little data, use average from last two periods:
  for (age in 1:(startestage-1)) {
    obsinc <- cases[age,(noobsper-1):noobsper]/pyr[age,(noobsper-1):noobsper]
    if (sum(is.na(obsinc))) {
      obsinc[is.na(obsinc)] <- 0 	
    }
    datatable[age,(noobsper+1):nototper] <- ((obsinc[,1]+obsinc[,2])/2)*pyr[age,(noobsper+1):nototper]
  }

  # For old agegroups, use AC model estimates:
  coefficients <- acproj.estimate.object$glm$coefficients
  # Age effects:
  noages    <- 19-startestage+1
  ageff <- coefficients[1:noages]
  # Cohort effects:
  coeff <- c(coefficients[(noages+1):(noages+midc-1)],0,coefficients[(noages+midc):length(coefficients)])
  lncoh <- length(coeff)
  # For new young cohorts in "nonewpred" by linear regression without the first and last 4 cohorts:
  ceff <- coeff[-c(1,2, lncoh-1,lncoh)]
  lnceff <- length(ceff)
  tt <- 1:lnceff 
  coeflm <- coef(lm(ceff ~ tt-1))
  ncoeff <- coeflm*c((lncoh-3):(lncoh-2), (cumsum(1- cuttrend)+lncoh-2) )
  coheff <- c(coeff[1:(lncoh-2)], ncoeff) # replace the last 2 cohort by estimates

  for (age in startestage:19) {
  # Age effect:
    agepar    <- as.numeric(ageff[age-startestage+1])
  # Cohort index: No. agegoups - age + period
    coh       <- (19-startestage) - (age-startestage) + (noperiod+1:nonewpred)
    cohpar <- coheff[coh]
  # Project age-specific rates:    
    if (acproj.estimate.object$linkfunc=="power5") { 
        rate <- (agepar + cohpar)^5      
     } else  if (acproj.estimate.object$linkfunc=="log") {
        rate <- exp(agepar + cohpar)      
     } else  if (acproj.estimate.object$linkfunc=="sqrt") {
         rate <- (agepar + cohpar)^2      
     } else { # identity link:
        rate <- (agepar + cohpar)      
     }
    datatable[age,(noobsper+1):nototper] <- rate*pyr[age,(noobsper+1):nototper]
  }

  # Structure and return results:
  res <- list(predictions=datatable,
              pyr=pyr, nopred=nonewpred, noperiod=acproj.estimate.object$noperiod,
              cuttrd=cuttrd, shortp=shortp, cuttrend=cuttrend, 
              gofpvalue=acproj.estimate.object$gofpvalue, 
              distribution=acproj.estimate.object$distribution,
              startestage=acproj.estimate.object$startestage,
              glm=acproj.estimate.object$glm)
  class(res) <- "acproj"
  attr(res,"Call") <- sys.call()
  return(res)
}


## Summary Projection Results by 5-year Period(ASR and total number):
acproj.getpred <- function(acproj.object,
                           incidence=T, standpop=NULL, excludeobs=F, byage, agegroups="all") {  
## Setting defaults:
  if (missing(byage)) {
     byage <- ifelse(is.null(standpop), T, F)
  }
  
## Checking input:
  if (class(acproj.object)!="acproj") {
     stop("Variable \"acproj.object\" must be of type \"acproj\"")	
  } 
  
  if ((!is.null(standpop)) && (!incidence)) {
     stop("\"standpop\" should only be used with incidence predictions (incidence=T)") 	
  }
 
  if (!is.null(standpop)) {
    if (round(sum(standpop),5)!=1) {
       stop("\"standpop\" must be of sum 1") 	
    }
    if ((length(standpop)!=length(agegroups)) && (agegroups[1]!="all")) {
       stop("\"standpop\" must be the same length as \"agegroups\"") 	
    }
    if (byage) {
       stop("\"standpop\" is only valid for \"byage=T\"") 
    }
  }
 
## Seting local data:
  datatable <- acproj.object$predictions
  pyr       <- data.frame(acproj.object$pyr)

## Secting agegroups:
  if (agegroups[1]!="all") {
     datatable <- datatable[agegroups, ]
     pyr       <- pyr[agegroups, ]
  }

## If needed; Standardize data and Collapse agegroups
  if (!is.null(standpop)) {
     datainc <- (datatable/pyr)*100000
     if (sum(is.na(datainc)) > 0) {
        datainc[is.na(datainc)] <- 0
     } 	
    res <- apply(datainc*standpop, 2, sum) 
  } else {
    if (!byage) {
       datatable <- apply(datatable, 2, sum)
       pyr       <- apply(pyr, 2, sum) 
    }
    if (incidence) {
       res <- (datatable/pyr)*100000
       if (sum(is.na(res)) > 0) {
          res[is.na(res)] <- 0
       } 	
    }  else {
          res <- datatable
    } 
  }

## Select data:
  if (excludeobs) {
     if (is.matrix(res)) {
        predstart <- dim(res)[2] - acproj.object$nopred + 1
        res <- res[ , predstart:(predstart+acproj.object$nopred-1)]
     }  else {
          predstart <- length(res) - acproj.object$nopred + 1
          res <- res[predstart:(predstart+acproj.object$nopred-1)]
     }
  } 
  
## Return data: 
  return(res)
}


## Summary Annual Projection Results:
acproj.getproj <- function(cdat, pdat, startp, acproj.object, standpop=NULL) {
  r0 <- acproj.getpred(acproj.object, incidence=T)
  outasp <- asrpy(r0, cdat, pdat, startp=startp, nagg=5)
  if (is.null(standpop)) {
    return(outasp)
  } else {
    outann <- asry(outasp, pdat, standpop=standpop)
    return(outann)
  }
}



## Summaring methods used for the projection:
summary.acproj <- function(acproj.object, printpred=F, printcall=F, digits=0) {
  method <- "Age-Cohort Model"
  if (class(acproj.object)!="acproj") {
    stop("Variable \"acproj.object\" must be of type \"acproj\"")	
  } 
  # Setting internal variables:
  obsto <- names(acproj.object$predictions)[dim(acproj.object$predictions)[2]-acproj.object$nopred]

  if (!is.null(acproj.object$gofpvalue)) {
    gofpvalue <- round(acproj.object$gofpvalue,4)
  } else { gofpvalue <- NA }
      
  # Print information about object: 
  if (printpred) {
    cat("Observed and predicted values:")
    cat("(observations up to",obsto,")\n")
    print(round(as.matrix(acproj.object$predictions),digits=digits))
    cat("\n")
  }
  cat("\nPrediction done with:\n")
    
  moptions <- matrix(NA,7,2)
  moptions[,1] <- c("Method:",
                    "Number of periods predicted (nopred):",
                    "Trend used in new cohort estimation (cuttrend):",
                    "Number of periods used in estimate (noperiod):",
                    "Distribution function of regression:",
                    "P-value for goodness of fit:",
                    "First age group estimated (startestage):")        
  moptions[,2] <- c(method,
                    acproj.object$nopred, 
                    paste(acproj.object$cuttrend, collapse=" , "),
                    acproj.object$noperiod, 
                    acproj.object$distribution,
                    gofpvalue,
                    acproj.object$startestage)

  maxl <- max(nchar(moptions[,1]))
    
  for (i in 1:dim(moptions)[1]) {
    spaces <- rep (" ", maxl-nchar(moptions[i,1])+2)
    cat(moptions[i,1], spaces, moptions[i,2], "\n", sep="")	
  }    
  if (printcall) {
    cat("\n  Call: ")
    dput(attr(acproj.object, "Call"))
  }
  invisible(acproj.object)
}

## Obtaining modeling info and parameter estimates:
glm.acproj <- function(acproj.object) {
     summary(acproj.object$glm)
}

## Ploting prejected standardized rates:
plot.acproj<- function(cdat, pdat, startp, acproj.object, standpop, startplot=1, 
                         xlab="Calendar Year", ylab="Rates",
                main="",labels=NULL,ylim=NULL,lty=c(1,3),col=c(1,1),new=T,...) {
  if (class(acproj.object)!="acproj") {
    stop("Variable \"acproj.object\" must be of type \"acproj\"")	
  } 
  
  # Reading & formating data:
  indat <- acproj.getproj(cdat, pdat, startp=startp, acproj.object, standpop=standpop)
  indata <- indat[, 1]
  indata <- indata[startplot:length(indata)]

  # Seting internal variables:
  obsy <- dim(cdat)[2]
  nopredy <- length(indata) - obsy 
  if (is.null(labels)) {
    labels <- row.names(indat)
  }

  # Create plots:
  maxx <- length(indata)
  if (new) { 	
    maxy <- max(indata)*(21/20)	
    if (is.null(ylim)) {
      ylim <- c(0,maxy)	
    }
    plot(c(1,maxx),ylim,type="n",ylab=ylab,xlab=xlab,axes=F,...)
    axis(2)
    axis(1,at=1:maxx,labels=labels)
    box()
    title(main)
  }
  lines(1:(maxx-nopredy),indata[1:(maxx-nopredy)],type="o",pch=20,lty=lty[1],col=col[1],...)
  lines((maxx-nopredy):maxx,indata[(maxx-nopredy):maxx],lty=lty[2],col=col[2],...)

  # Returning object as invisible
  invisible(acproj.object)
}



############################################################################## 
#    HYBDPROJ:  R functions for projection of cancer incidence/mortality     #
# using the modified Hybrid methods                                          #
# Modified Hybrid by adding choice of age-model, cut-trend parameter,        #
# and power 5 link function                                                  #
# Written by: Dr. Zhenguo QIU                                                #
# QA by: Zhichang JIANG, 2011                                                #
# License:	GNU version 1                                                    #
############################################################################## 

## Hybrid method control program:
hybdproj <- function(cdat, pdat, projfor="incidence", nagg=NULL, ncase=NULL, 
                     cuttrd=0.04, shortp=0, linkfunc="power5", pD=0.05, pGOF=0.05) {
#cdat: 19 (age groups) x N (years) historical cancer data, N >= 15
#pdat: 19 (age groups) x (N + M) (years) observed and projected population, M <= 25
#nagg: number of years for data aggregation (by years), default: 1, annual data.
#ncases: minimum number of cancer cases/deaths for splitting data.
#cuttrd: degenerating percent of trends per year after 5 years (shortp=0)
#        or the first projection year.
#linkfunc: link function, default is power5, can be log, sqrt and identity

# Define number of cases for data splitting:
    if (is.null(ncase)) {
       if (projfor=="incidence") {
          ncase <- 1
       } else {
          ncase <- 3/5
       } 
    }

# Define number of years for data aggregation:
    if (is.null(nagg)) {
      masr <- mean(obasr(cdat, pdat)[,1])
      if (projfor=="incidence") {   
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
  
## Data aggregation:
  aggdata <- datagg(cdat, pdat, nagg)
  cases <- aggdata$cases
  pyr <- aggdata$pyr

## Setting default and checking data:
  # minimum number of periods in observed data
  percases <- dim(cases)[2]   
  if (percases < 5) { 
     stop("Minimum number of period is 5 (5*nagg years) in \"cases\"") 
  }
   
  # maximum number of periods for projection base:
  noperiod <- percases 
  if (percases > 25) {
    noperiod <- 25 # maximum points is 25, cut-off ancient data:
    cases <- cases[ , -c(1:(percases-25))]
    pyr <- pyr[ , -c(1:(percases-25))]
  }
   
## Perform estimation and prediction:
  est  <- hybdproj.estimate(cases=cases, pyr=pyr, nagg=nagg, ncase=ncase, 
                            linkfunc=linkfunc, pD=pD, pGOF=pGOF)
  pred <- hybdproj.prediction(hybdproj.estimate.object=est, cuttrd=cuttrd, shortp=shortp)
  return(pred)
}


## Fitting models for data B:
hybdproj.estimate <- function(cases, pyr, nagg, ncase, linkfunc="power5", pD=0.05, pGOF=0.05) {
## Checking data
  if ( dim(cases)[2]>dim(pyr)[2]) {
    stop("\"pyr\" must include information about all periods in \"cases\"")	
  }
  
  if (dim(pyr)[2]==dim(cases)[2]) {
    stop("\"pyr\" must include information on future rates")	
  }

  if ((dim(pyr)[2]-dim(cases)[2]) > (30/nagg)) {
    stop("Package can not project more than 30 years")	
  }

  if (dim(cases)[2]< 5) {
    stop("\"noperiod\" must be 5 or larger")   	
  }  

## Splitting data by age groups: A for average method, B for regression.
   nage <- 1:19
   mcase <- apply(cases, 1, mean) # calculate the mean cases
   caseage <- as.data.frame(cbind(nage, cases))
   pyrage <- as.data.frame(cbind(nage, pyr))
   casesA <- caseage[mcase < (nagg*ncase), ]
   pyrA <- pyrage[mcase < (nagg*ncase), ]
   casesB <- caseage[mcase >= (nagg*ncase), ]
   pyrB <- pyrage[mcase >= (nagg*ncase), ]

## Setting internal variables:
   dnoperiods <- dim(casesB)[2] - 1
   dnoagegr   <- dim(casesB)[1]
   agpreg <- casesB[ , 1]
   agpave <- casesA[ , 1]
#
   if (dnoagegr < 2) {
      stop("\"dnoagegr\" must be 2 or larger")   	
   }  

## Transform data format:
  ageno    <- rep(agpreg, dnoperiods)
  periodno <- sort(rep(1:dnoperiods,dnoagegr))
  y        <- c(as.matrix(pyrB[,2:(dnoperiods+1)]))
  apdata   <- data.frame(Age=ageno, Period=periodno, 
                       Cases=c(as.matrix(casesB[,-1])), y=y)

## Select cutting year:
  if (dim(casesB)[2] == 5) {
      apdatan <- apdata
	lastper <- 5
      projbase <- 5*nagg
  } else {
      mod1 <- glm(Cases~as.factor(Age)+Period+offset(log(y))-1, 
                        data=apdata, family=poisson) 
      mod2 <- glm(Cases~as.factor(Age)+as.factor(Age)*Period+offset(log(y))-1, 
                        data=apdata, family=poisson)   
      pdiff <- 1 - pchisq((mod1$deviance-mod2$deviance), (mod1$df.residual-mod2$df.residual))
      if (is.null(pdiff)) {
           mod <- "common"
      } else { 
        if (pdiff < pGOF) {
             mod <- "age-specific"
        } else {
             mod <- "common"
        }
      }
   
  # likelihood searching:
      trydat <- apdata
      likeh <- rep(0, (dnoperiods-5))
      for (i in 1: (dnoperiods-5)) {
         trydat$Period[trydat$Period <= i] <- 0
         if (mod=="common") {   
           likeh[i] <- glm(Cases~as.factor(Age)+Period+offset(log(y))-1, 
                                 data=trydat, family=poisson)$deviance/(-2)
         } else {
           likeh[i] <- glm(Cases~as.factor(Age)+as.factor(Age)*Period+offset(log(y))-1, 
                                 data=trydat, family=poisson)$deviance/(-2)           
         }
      }
      cuty <- which(likeh==max(likeh))   
      apdatan <- apdata[apdata$Period >= cuty, ] 
      if (cuty > 1) {
        apdatan$Period <- apdatan$Period - cuty + 1
      }   
	lastper <- length(cuty:dnoperiods)
      projbase <- lastper*nagg
  }

## finalize model using the projection base:
   mode1 <- glm(Cases~as.factor(Age)+Period+offset(log(y))-1, data=apdatan, family=poisson) 
   mode2 <- glm(Cases~as.factor(Age)+as.factor(Age)*Period+offset(log(y))-1, 
                      data=apdatan, family=poisson)   
   pdiffn <- 1 - pchisq((mod1$deviance-mod2$deviance), (mod1$df.residual-mod2$df.residual))
   pd1 <- summary(mode1)$coef["Period", 4]
   pd2 <- 1 - pchisq(mode2$deviance, mode2$df.residual)

   if (is.null(pdiffn)) {
      if (pd1 > pD) {
         fmodel <- "average"
      } else {
         fmodel <- "common-trend"
      } 
   } else {
     if (pdiffn > pGOF) {
        if (pd1 > pD) {
           fmodel <- "average"
        } else {
           fmodel <- "common-trend"
        } 
     } else {
        if (pd2 > pD) {
           fmodel <- "age-specific"
        } else {
           fmodel <- "nba-specific"
        }
     }
    }

## Estimation:
  if (fmodel=="common-trend") {
    if (linkfunc=="power5") {
      # Creation of power5 link for poisson family:
      y <- apdatan$y
      power5link <- poisson()
      power5link$link <- "0.2 root link Poisson family"
      power5link$linkfun <- function(mu)  { (mu/y)^0.2 }
      power5link$linkinv <- function(eta) { pmax(.Machine$double.eps, y*eta^5) }
      power5link$mu.eta <- function(eta)  { pmax(.Machine$double.eps, 5*y*eta^4) }
      #
      res.glm <- glm(Cases ~ as.factor(Age)+Period - 1,
                             data=apdatan, family=power5link)
    } else  if (linkfunc=="log") {
      res.glm <- glm(Cases ~ as.factor(Age)+Period+ offset(log(y))-1,
                             data=apdatan, family=poisson(link=log))  	
    } else  if (linkfunc=="sqrt") {
      res.glm <- glm(Cases/y ~ as.factor(Age)+Period -1,
                               data=apdatan, family=poisson(link=sqrt))  
    } else  if (linkfunc=="identity") {
      res.glm <- glm(Cases/y ~ as.factor(Age)+Period -1,
                               data=apdatan, family=poisson(link=identity))  
    } else {
      stop("Unknown \"linkfunc\"")	
    }
    pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df.residual)
  } else if (fmodel=="age-specific") {
    if (linkfunc=="power5") {
      # Creation of power5 link for poisson family:
      y <- apdatan$y
      power5link <- poisson()
      power5link$link <- "0.2 root link Poisson family"
      power5link$linkfun <- function(mu)  { (mu/y)^0.2 }
      power5link$linkinv <- function(eta) { pmax(.Machine$double.eps, y*eta^5) }
      power5link$mu.eta <- function(eta)  { pmax(.Machine$double.eps, 5*y*eta^4) }
      #
      res.glm <- glm(Cases ~ as.factor(Age)+as.factor(Age)*Period - 1 - Period,
                             data=apdatan, family=power5link)
    } else  if (linkfunc=="log") {
      res.glm <- glm(Cases ~ as.factor(Age)+as.factor(Age)*Period+ offset(log(y))-1-Period,
                             data=apdatan, family=poisson(link=log))  	
    } else  if (linkfunc=="sqrt") {
      res.glm <- glm(Cases/y ~ as.factor(Age)+as.factor(Age)*Period - 1 - Period,
                               data=apdatan, family=poisson(link=sqrt))  
    } else  if (linkfunc=="identity") {
      res.glm <- glm(Cases/y ~ as.factor(Age)+as.factor(Age)*Period - 1 - Period,
                               data=apdatan, family=poisson(link=identity))  
    } else {
      stop("Unknown \"linkfunc\"")	
    }
    pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df.residual)
  } else if (fmodel=="nba-specific") {
  options(warn=-1)
    if (linkfunc=="power5") {
    # Make power5 link function for negative binomial family:
      y <- apdatan$y
      glmnb <- glm.nb(Cases ~ as.factor(Age)+as.factor(Age)*Period - 1-Period
                            + offset(log(y)), data=apdatan, link=log)
      theta <- as.numeric(theta.md(apdatan$Cases, fitted(glmnb), dfr = df.residual(glmnb)))
      nbpower5link <- negative.binomial(theta)
      nbpower5link$link <- "0.2 root link negative.binomial(theta) family"
      nbpower5link$linkfun <- function(mu)  { (mu/y)^0.2 }
      nbpower5link$linkinv <- function(eta) { pmax(.Machine$double.eps, y*eta^5) }
      nbpower5link$mu.eta <- function(eta)  { pmax(.Machine$double.eps, 5*y*eta^4) }
      #
      res.glm <- glm(Cases ~ as.factor(Age)+as.factor(Age)*Period - 1 - Period,
                             data=apdatan, family=nbpower5link)
    } else  if (linkfunc=="log") {
      res.glm <- glm.nb(Cases ~ as.factor(Age)+as.factor(Age)*Period+ offset(log(y))-1-Period,
                                data=apdatan, link=log)  	
    } else  if (linkfunc=="sqrt") {
      res.glm <- glm.nb(Cases/y ~ as.factor(Age)+as.factor(Age)*Period - 1 - Period,
                                  data=apdatan, link=sqrt)  
    } else  if (linkfunc=="identity") {
      res.glm <- glm.nb(Cases/y ~ as.factor(Age)+as.factor(Age)*Period - 1 - Period,
                                  data=apdatan, link=identity) 
    } else {
      stop("Unknown \"linkfunc\"")	
    }
    pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df.residual)
  options(warn=0)
  } else { ## Using age-model -- average method: 
    if (linkfunc=="power5") {
      # Creation of power5 link for poisson family:
      y <- apdatan$y
      power5link <- poisson()
      power5link$link <- "0.2 root link Poisson family"
      power5link$linkfun <- function(mu)  { (mu/y)^0.2 }
      power5link$linkinv <- function(eta) { pmax(.Machine$double.eps, y*eta^5) }
      power5link$mu.eta <- function(eta)  { pmax(.Machine$double.eps, 5*y*eta^4) }
      #
      res.glm <- glm(Cases ~ as.factor(Age) - 1,
                             data=apdatan, family=power5link)
    } else  if (linkfunc=="log") {
      res.glm <- glm(Cases ~ as.factor(Age) + offset(log(y))-1,
                             data=apdatan, family=poisson(link=log))  	
    } else  if (linkfunc=="sqrt") {
      res.glm <- glm(Cases/y ~ as.factor(Age) -1,
                               data=apdatan, family=poisson(link=sqrt))  
    } else  if (linkfunc=="identity") {
      res.glm <- glm(Cases/y ~ as.factor(Age) -1,
                               data=apdatan, family=poisson(link=identity))  
    } else {
      stop("Unknown \"linkfunc\"")	
    }
    pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df.residual)

  }

## Set class and return results
  res <- list(glm=res.glm, cases=cases, pyr=pyr, agrpave=agpave, lastper=lastper, cuty=cuty,
              noperiod=lastper, noyearagg=nagg, nocaseagp=ncase, linkfunc=linkfunc, 
              agrpmod=agpreg, projbase=projbase, finalmod=fmodel, gofpvalue=pvalue)
  class(res) <- "hybdproj.estimate"
  attr(res, "Call") <- sys.call()
  return(res)
}


## Project age-specific numbers:
hybdproj.prediction <- function(hybdproj.estimate.object, cuttrd=0.04, shortp=0) {
## running conditions:  
  if (class(hybdproj.estimate.object)!="hybdproj.estimate") {
    stop("Variable \"hybdproj.estimate.object\" must be of type \"hybdproj.estimate\"")	
  } 
## Setting local variables:
  cases     <- hybdproj.estimate.object$cases
  pyr	      <- hybdproj.estimate.object$pyr
  nagg      <- hybdproj.estimate.object$noyearagg
  ncase     <- hybdproj.estimate.object$nocaseagp
  projbase  <- hybdproj.estimate.object$projbase
  agpreg    <- hybdproj.estimate.object$agrpmod
  agpave    <- hybdproj.estimate.object$agrpave
  fmodel    <- hybdproj.estimate.object$finalmod
  noperiod  <- hybdproj.estimate.object$noperiod
  lastper   <- hybdproj.estimate.object$lastper
#  linkfunc  <- hybdproj.estimate.object$linkfunc
  nototper  <- dim(pyr)[2]
  noobsper  <- dim(cases)[2]
  nonewpred <- nototper - noobsper
## Define degenerating parameters:
  cuttrend <- rep(shortp, nonewpred)
  if (nagg==1) {
      if (nonewpred <= 5) {
          cuttrend <- rep(shortp, nonewpred) # first 5 years
      } else { 
          for (i in 6:nonewpred) { 
             cuttrend[i] <- shortp + (i-5)*cuttrd
          }
      }  
      cuttrend[cuttrend > 1] <- 1  
  }
  if (nagg==2) {
      if (nonewpred <= 3) {
          cuttrend <- rep(shortp, nonewpred) # first 6 years
      } else { 
          for (i in 4:nonewpred) { 
             cuttrend[i] <- shortp + (i-3)*2*cuttrd
          }
      }  
      cuttrend[cuttrend > 1] <- 1  
  }
  if (nagg==3) {
      if (nonewpred <= 2) {
          cuttrend <- rep(shortp, nonewpred) # first 6 years
      } else { 
          for (i in 3:nonewpred) { 
             cuttrend[i] <- shortp + (i-2)*3*cuttrd
          }
      }  
      cuttrend[cuttrend > 1] <- 1  
  }
  if (nagg==4) {
      if (nonewpred == 1) {
          cuttrend <- shortp # first 4 years
      } else { 
          for (i in 2:nonewpred) { 
             cuttrend[i] <- shortp + (i-1)*4*cuttrd
          }
      }  
      cuttrend[cuttrend > 1] <- 1  
  }
  if (nagg==5) {
      if (nonewpred == 1) {
          cuttrend <- shortp # first 5 years
      } else { 
          for (i in 2:nonewpred) { 
             cuttrend[i] <- shortp + (i-1)*5*cuttrd
          }
      }  
      cuttrend[cuttrend > 1] <- 1  
  }

## Future period values:
  driftmp   <- cumsum(1-cuttrend)
## define years or period names:  
  if (is.data.frame(pyr)) {
    years <- names(pyr)
  } else {
    if (is.null(dimnames(pyr))) {
      years <- paste("Periode",1:nototper)	
    } else {
      years <- dimnames(pyr)[[2]]
    }	
  }  
## Making data object:
  datatable <- matrix(NA, 19, nototper)
  # fill in observed cases:
  datatable[ , 1:(nototper-nonewpred)] <- as.matrix(cases)
  datatable <- data.frame(datatable)
  row.names(datatable) <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34",
                            "35-39","40-44","45-49","50-54","55-59","60-64",
                            "65-69","70-74","75-79","80-84","85-89","90+")
  names(datatable) <- years 
## calculate average rates:
  obsrat <- data.frame(matrix(0, dim(cases)[1], 5))  
  for (age in 1:(dim(cases)[1]) ) {
       obsrat[age, ] <- as.vector(cases[age,(noobsper-4):noobsper])/as.vector(
                                                pyr[age,(noobsper-4):noobsper])
  }
  obsrate <- apply(obsrat, 1, mean)

## extract parameter estimation and project age-specific numbers:
  if (fmodel=="age-specific" || fmodel=="nba-specific") {
     coef <- hybdproj.estimate.object$glm$coef
     m.eff <- cbind(coef[1:(length(coef)/2)],coef[(length(coef)/2+1):length(coef)],agpreg)
     row.names(m.eff) <- NULL
     avef <- cbind(rep(NA,length(agpave)),rep(0,length(agpave)),agpave)
     mcoef <- data.frame(rbind(avef, m.eff))
     acoef <- mcoef[with(mcoef, order(mcoef[ , 3])), ]
     row.names(acoef) <- NULL
     colnames(acoef) <- c("a.eff", "p.eff", "agrp")
  # Project age-specific numbers:
    for (age in 1:19) {
      if (is.na(acoef$a.eff[acoef$agrp==age])) {
         rate <- rep(obsrate[age], length(driftmp))
      } else {
        if (hybdproj.estimate.object$linkfunc=="power5") { 
           rate <- (acoef$a.eff[age]+acoef$p.eff[age]*(lastper+driftmp))^5
        } else  if (hybdproj.estimate.object$linkfunc=="log") {
           rate <- exp(acoef$a.eff[age]+acoef$p.eff[age]*(lastper+driftmp))
        } else  if (hybdproj.estimate.object$linkfunc=="sqrt") {
           rate <- (acoef$a.eff[age]+acoef$p.eff[age]*(lastper+driftmp))^2     
        } else { # identity link:
           rate <- (acoef$a.eff[age]+acoef$p.eff[age]*(lastper+driftmp))
        }
      }
      datatable[age,(noobsper+1):nototper] <- rate*pyr[age,(noobsper+1):nototper]
    }
  } else if (fmodel=="common-trend") {
     coef <- hybdproj.estimate.object$glm$coef
     m.eff <- cbind(coef[1:(length(coef)-1)],rep(coef[length(coef)],length(agpreg)),agpreg)
     row.names(m.eff) <- NULL
     avef <- cbind(rep(NA,length(agpave)),rep(0,length(agpave)),agpave)
     mcoef <- data.frame(rbind(avef, m.eff))
     acoef <- mcoef[with(mcoef, order(mcoef[ , 3])), ]
     row.names(acoef) <- NULL
     colnames(acoef) <- c("a.eff", "p.eff", "agrp")
  # Project age-specific numbers:
    for (age in 1:19) {
      if (is.na(acoef$a.eff[acoef$agrp==age])) {
         rate <- rep(obsrate[age], length(driftmp))
      } else {
        if (hybdproj.estimate.object$linkfunc=="power5") { 
           rate <- (acoef$a.eff[age]+acoef$p.eff[age]*(lastper+driftmp))^5
        } else  if (hybdproj.estimate.object$linkfunc=="log") {
           rate <- exp(acoef$a.eff[age]+acoef$p.eff[age]*(lastper+driftmp))
        } else  if (hybdproj.estimate.object$linkfunc=="sqrt") {
           rate <- (acoef$a.eff[age]+acoef$p.eff[age]*(lastper+driftmp))^2     
        } else { # identity link:
           rate <- (acoef$a.eff[age]+acoef$p.eff[age]*(lastper+driftmp))
        }
      }
      datatable[age,(noobsper+1):nototper] <- rate*pyr[age,(noobsper+1):nototper]
    }
  } else { # Average method
     coef <- hybdproj.estimate.object$glm$coef
     m.eff <- cbind(coef[1:length(coef)],rep(0,length(agpreg)),agpreg)
     row.names(m.eff) <- NULL
     avef <- cbind(rep(NA,length(agpave)),rep(0,length(agpave)),agpave)
     mcoef <- data.frame(rbind(avef, m.eff))
     acoef <- mcoef[with(mcoef, order(mcoef[ , 3])), ]
     row.names(acoef) <- NULL
     colnames(acoef) <- c("a.eff", "p.eff", "agrp")
  # Project age-specific numbers:
    for (age in 1:19) {
      if (is.na(acoef$a.eff[acoef$agrp==age])) {
         rate <- rep(obsrate[age], length(driftmp))
      } else {
        if (hybdproj.estimate.object$linkfunc=="power5") { 
           rate <- (acoef$a.eff[age]+acoef$p.eff[age]*(lastper+driftmp))^5
        } else  if (hybdproj.estimate.object$linkfunc=="log") {
           rate <- exp(acoef$a.eff[age]+acoef$p.eff[age]*(lastper+driftmp))
        } else  if (hybdproj.estimate.object$linkfunc=="sqrt") {
           rate <- (acoef$a.eff[age]+acoef$p.eff[age]*(lastper+driftmp))^2     
        } else { # identity link:
           rate <- (acoef$a.eff[age]+acoef$p.eff[age]*(lastper+driftmp))
        }
      }
      datatable[age,(noobsper+1):nototper] <- rate*pyr[age,(noobsper+1):nototper]
    }
  }

## Structure and return results:

   res <- list(predictions=datatable, pyr=pyr, cuttrd=cuttrd, shortp=shortp, cuttrend=cuttrend, 
              nopred=nonewpred, noperiod=hybdproj.estimate.object$noperiod,
              lastperiod=lastper, noobsper=noobsper, nototper=nototper,
              noyearagg=hybdproj.estimate.object$noyearagg,
              nocaseagp=hybdproj.estimate.object$nocaseagp,
              agrpave=hybdproj.estimate.object$agrpave,
              agrpmod=hybdproj.estimate.object$agrpmod, 
              linkfunc=hybdproj.estimate.object$linkfunc, 
              projbase=hybdproj.estimate.object$projbase,
              finalmod=hybdproj.estimate.object$finalmod,
              gofpvalue=hybdproj.estimate.object$gofpvalue, 
              glm=hybdproj.estimate.object$glm)

  class(res) <- "hybdproj"
  attr(res,"Call") <- sys.call()
  return(res)
}


## Summary Projection Results (ASR and total number):
hybdproj.getpred <- function(hybdproj.object,
                    incidence=T, standpop=NULL, excludeobs=F, byage, agegroups="all") {  
## Setting defaults:
  if (missing(byage)) {
     byage <- ifelse(is.null(standpop), T, F)
  }
  
## Checking input:
  if (class(hybdproj.object)!="hybdproj") {
     stop("Variable \"hybdproj.object\" must be of type \"hybdproj\"")	
  } 
  
  if ((!is.null(standpop)) && (!incidence)) {
     stop("\"standpop\" should only be used with incidence predictions (incidence=T)") 	
  }
 
  if (!is.null(standpop)) {
    if (round(sum(standpop),5)!=1) {
       stop("\"standpop\" must be of sum 1") 	
    }
    if ((length(standpop)!=length(agegroups)) && (agegroups[1]!="all")) {
       stop("\"standpop\" must be the same length as \"agegroups\"") 	
    }
    if (byage) {
       stop("\"standpop\" is only valid for \"byage=T\"") 
    }
  }
 
## Seting local data:
  datatable <- hybdproj.object$predictions
  pyr       <- data.frame(hybdproj.object$pyr)

## Secting agegroups:
  if (agegroups[1]!="all") {
     datatable <- datatable[agegroups, ]
     pyr       <- pyr[agegroups, ]
  }

## If needed; Standardize data and Collapse agegroups
  if (!is.null(standpop)) {
     datainc <- (datatable/pyr)*100000
     if (sum(is.na(datainc)) > 0) {
        datainc[is.na(datainc)] <- 0
     } 	
    res <- apply(datainc*standpop, 2, sum) 
  } else {
    if (!byage) {
       datatable <- apply(datatable, 2, sum)
       pyr       <- apply(pyr, 2, sum) 
    }
    if (incidence) {
       res <- (datatable/pyr)*100000
       if (sum(is.na(res)) > 0) {
          res[is.na(res)] <- 0
       } 	
    }  else {
          res <- datatable
    } 
  }

## Select data:
  if (excludeobs) {
     if (is.matrix(res)) {
        predstart <- dim(res)[2] - hybdproj.object$nopred + 1
        res <- res[ , predstart:(predstart+hybdproj.object$nopred-1)]
     }  else {
          predstart <- length(res) - hybdproj.object$nopred + 1
          res <- res[predstart:(predstart+hybdproj.object$nopred-1)]
     }
  } 
  
## Return data: 
  return(res)
}


## Summary Annual Projection Results:
hybdproj.getproj <- function(cdat, pdat, startp, hybdproj.object, standpop=NULL, Ave5=NULL, sum5=NULL) {
  finalmod <- hybdproj.object$finalmod
  r0 <- hybdproj.getpred(hybdproj.object, incidence=T)
  nagg <- hybdproj.object$noyearagg
  outasp <- asrpy(r0, cdat, pdat, startp=startp, nagg=nagg)
  if (finalmod != "average") {
      if (is.null(standpop)) {
           return(outasp)
      } else {
           outann <- asry(outasp, pdat, standpop=standpop)
           return(outann)
      } 
  } else {
      if (is.null(Ave5)) {
          if (is.null(standpop)) {
               return(outasp)
          } else {
               outann <- asry(outasp, pdat, standpop=standpop)
               return(outann)
          } 
      } else {
          mod <- ave5proj(cdat, pdat, startp, sum5=sum5)
          outasp <- mod$agsproj
          if (is.null(standpop)) {
              return(outasp)
          } else {
              outann <- asry(outasp, pdat, standpop=standpop)
              return(outann)
          } 
      }
  }
}


## summary results:
summary.hybdproj <- function(hybdproj.object, printpred=F, printcall=F, digits=0) {
  method <- "Hybrid approach"
  if (class(hybdproj.object)!="hybdproj") {
    stop("Variable \"hybdproj.object\" must be of type \"hybdproj\"")	
  }  
  # Setting internal variables:
  nototper <- hybdproj.object$nototper
  noobsper <- hybdproj.object$noobsper
  nopred <- hybdproj.object$nopred
  noypred <- (hybdproj.object$nopred)*(hybdproj.object$noyearagg)
  noperiod <- hybdproj.object$noperiod
  obsto <- names(hybdproj.object$predictions)[dim(hybdproj.object$predictions)[2]-nopred]
  predcases <- hybdproj.object$predictions[ , (noobsper+1): (nototper)]   
  if (!is.null(hybdproj.object$gofpvalue)) {
    gofpvalue <- round(hybdproj.object$gofpvalue,4)
  } else { gofpvalue <- NA }      
  # Print information about object: 
  if (printpred) {
    cat("Predicted number of cases or deaths:")
    cat("(observations up to",obsto,")\n")
    print(round(as.matrix(predcases), digits=digits))
    cat("\n")
  }
  cat("\nPrediction done with:\n")    
  moptions <- matrix(NA,12,2)
  moptions[,1] <- c("Method:",
                    "Number of prediction years:",
                    "First period cutting trend:",
                    "Degenerating trend per year:",
                    "Projection base (years):",
                    "Aggregating years (nagg):",
                    "Age-cases per year (ncase):", 
                    "Model for regression:",
                    "Link function for GLM:",
                    "P-value for goodness of fit:",
                    "Age group for regression:",
                    "Age group for average method:")        
  moptions[,2] <- c(method,
                    noypred, 
                    hybdproj.object$shortp,
                    hybdproj.object$cuttrd,
                    hybdproj.object$projbase, 
                    hybdproj.object$noyearagg,
                    hybdproj.object$nocaseagp,
                    hybdproj.object$finalmod,
                    hybdproj.object$linkfunc,
                    gofpvalue,
                    paste(hybdproj.object$agrpmod, collapse=","),
                    paste(hybdproj.object$agrpave, collapse=",") )
  maxl <- max(nchar(moptions[,1]))    
  for (i in 1:dim(moptions)[1]) {
    spaces <- rep (" ", maxl-nchar(moptions[i,1])+2)
    cat(moptions[i,1], spaces, moptions[i,2], "\n", sep="")	
  }
  if (printcall) {
    cat("\n  Call: ")
    dput(attr(hybdproj.object, "Call"))
  }
  invisible(hybdproj.object)
}

glm.hybdproj <- function(hybdproj.object) {
     summary(hybdproj.object$glm)
}

plot.hybdproj<- function(cdat, pdat, startp, hybdproj.object, standpop, startplot=1, 
                         xlab="Calendar Year", ylab="Rates",
                main="",labels=NULL,ylim=NULL,lty=c(1,3),col=c(1,1),new=T,...) {
  if (class(hybdproj.object)!="hybdproj") {
    stop("Variable \"hybdproj.object\" must be of type \"hybdproj\"")	
  } 
  
  # Reading & formating data:
  indat <- hybdproj.getproj(cdat, pdat, startp=startp, hybdproj.object, standpop=standpop)
  indata <- indat[, 1]
  indata <- indata[startplot:length(indata)]

  # Seting internal variables:
  obsy <- dim(cdat)[2]
  nopredy <- length(indata) - obsy 
  if (is.null(labels)) {
    labels <- row.names(indat)
  }

  # Create plots:
  maxx <- length(indata)
  if (new) { 	
    maxy <- max(indata)*(21/20)	
    if (is.null(ylim)) {
      ylim <- c(0,maxy)	
    }
    plot(c(1,maxx),ylim,type="n",ylab=ylab,xlab=xlab,axes=F,...)
    axis(2)
    axis(1,at=1:maxx,labels=labels)
    box()
    title(main)
  }
  lines(1:(maxx-nopredy),indata[1:(maxx-nopredy)],type="o",pch=20,lty=lty[1],col=col[1],...)
  lines((maxx-nopredy):maxx,indata[(maxx-nopredy):maxx],lty=lty[2],col=col[2],...)

  # Returning object as invisible
  invisible(hybdproj.object)
}



############################################################################## 
#    Average5:  R functions for projection of cancer incidence/mortality     #
# using the average methods based on the recent 5 years data                 #
# (i) average numbers and population sizes then calculate rates (default), or#
# (ii) average the calculated yearly rates (sum5=T)                          #
# Written by: Dr. Zhenguo QIU                                                #
# License:	GNU version 1                                                    #
############################################################################## 

ave5proj <- function(cdat, pdat, startp, sum5=NULL){
    nc <- dim(cdat)[2]
    np <- dim(pdat)[2]
    noypred <- np - nc
  # fill in observed rates:
    obasr <- matrix(NA, 19, nc)
    for (i in 1:nc) {
        obasr[, i] <- 100000*cdat[, i]/pdat[, i]
    }
    datatab <- matrix(NA, 19, np)
    datatab[ , 1:nc] <- as.matrix(obasr)
    datatab <- data.frame(datatab)
    row.names(datatab) <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
                            "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+")
    colnames(datatab) <- (startp-nc):(startp+(np-nc)-1) 
  # Calculate predictions in age-specific rates:
    for (age in 1:19) {
        if (is.null(sum5)) {
            obsinc <- apply(cdat[age,(nc-5):nc], 1, sum)/apply(pdat[age,(nc-5):nc], 1, sum)
        } else {
            obsinc <- cdat[age,(nc-5):nc]/pdat[age,(nc-5):nc]
        }
        if (sum(is.na(obsinc))) {
            obsinc[is.na(obsinc)] <- 0 	
        }
        if (is.null(sum5)) {
            datatab[age,(nc+1):np]<-100000*obsinc
        } else {
            datatab[age,(nc+1):np]<-100000*(obsinc[,1]+obsinc[,2]+obsinc[,3]+obsinc[,4]+obsinc[,5])/5
        }
    }
    res <- list(agsproj=datatab, startp=startp, sum5=sum5, noypred=noypred)
    class(res) <- "ave5proj"
    attr(res,"Call") <- sys.call()
    return(res)
}

ave5proj.getproj <- function(pdat, ave5proj.object, standpop=NULL) {
  outasp <- ave5proj.object$agsproj
  if (is.null(standpop)) {
    return(outasp)
  } else { 
    annproj <- asry(outasp, pdat, standpop=standpop)
    return(annproj)
  }
}


## summary methods:
summary.ave5proj <- function(ave5proj.object, printcall=F, digits=0) {
  if (class(ave5proj.object)!="ave5proj") {
    stop("Variable \"ave5proj.object\" must be of type \"ave5proj\"")	
  }  
  method <- "Five-Year Average"
  if (is.null(ave5proj.object$sum5)) {
      sum5 <- "5-year period"
  } else {
      sum5 <- "average yearly-rates"
  }
  noypred <- ave5proj.object$noypred
  
  cat("\nPrediction done with:\n")    
  moptions <- matrix(NA,3,2)
  moptions[,1] <- c("Method:",
                    "Age-Specific Rate by:",
                    "Number of prediction years:" )     
  moptions[,2] <- c(method,
                    sum5,
                    noypred) 

  maxl <- max(nchar(moptions[,1]))    
  for (i in 1:dim(moptions)[1]) {
    spaces <- rep (" ", maxl-nchar(moptions[i,1])+2)
    cat(moptions[i,1], spaces, moptions[i,2], "\n", sep="")	
  }
  if (printcall) {
    cat("\n  Call: ")
    dput(attr(ave5proj.object, "Call"))
  }
  invisible(ave5proj.object)
}



###############################################################
## Complimentary Functions:                                   #
## Written by Dr. Zhenguo QIU                                 #
###############################################################

## Data aggregation function: 
##aggregating cases and population by nagg=x years
datagg <- function(cdat, pdat, nagg) {
  # aggregating cancer data:
   ny <- dim(cdat)[2]  # of years for observed data
   nper <- floor(ny/nagg) # of period if aggregated by nagg    
   if (ny==nper*nagg) {
     cdat0 <- cdat
   } else { # truncate the beginning years less than nagg 
     cdat0 <- cdat[ , -c(1:(ny-nper*nagg))] 
   }
   pern <- sort(rep(1:nper, nagg)) # period vector
   cases <- as.data.frame(t(aggregate(t(cdat0), list(Period=pern), sum)))[-1, ]
   colnames(cases) <- 1:nper
   rownames(cases) <- 1:19

  # aggregating population data:
   pdat1 <- pdat[ , 1:ny] # observed population
   if (ny==nper*nagg) {
     pdat0 <- pdat1
   } else { # truncate the beginning years less than nagg 
     pdat0 <- pdat1[ , -c(1:(ny-nper*nagg))] 
   }
   pyr1 <- as.data.frame(t(aggregate(t(pdat0), list(Period=pern), sum)))[-1, ]

  # aggregating projection population data:
   nypop <- dim(pdat)[2]
   pdat2 <- pdat[ , c((ny+1):nypop)] # projected population
   nyp <- nypop - ny
   nperp <- floor(nyp/nagg) # of period if aggregated by nagg    
   if (nyp==nperp*nagg) {
     pdatn <- pdat2
   } else { # truncate the end years great than nagg 
     pdatn <- pdat2[ , -c((nperp*nagg+1):nyp)] 
   }
   pernp <- sort(rep(1:nperp, nagg)) # period vector
   pyr2 <- as.data.frame(t(aggregate(t(pdatn), list(Period=pernp), sum)))[-1, ]

  # combine population data:
   pyr <- as.data.frame(cbind(pyr1, pyr2))
   colnames(pyr) <- 1:(nper+nperp)
   rownames(pyr) <- 1:19
   return(list(cases=cases, pyr=pyr))
}


## Linear interpolation:
## Convert to annual age-specific rates
## by linear interpolation for each two-points segment
## if rate is projected age-specific rates in period:

asrpy <- function(rate, cdat, pdat, startp, nagg) {
# nagg: number of years used for aggregation: 1, 2, ..., 5
# r0: 19*m matrix, age-specific rate of m periods
    nt <- dim(cdat)[2]
    nperd <- floor(nt/nagg) # of observed periods if aggregated by nagg
    np <- dim(pdat)[2]
    if (nperd <= 25) {
       nypop <- dim(pdat)[2]  # of years for observed and projected data
       ny <- dim(cdat)[2]  # of years for observed data
       nyp <- nypop - ny  # of years for projection
       pdat0 <- pdat
    } else {
       ny <- 25
       nypop <- dim(pdat)[2] - (nt - 25)
       nyp <- nypop - ny
       pdat0 <- pdat[, -c(1:(nt - 25))]
    }

    nper <- floor(ny/nagg) # of observed periods if aggregated by nagg 
    nperp <- floor(nyp/nagg) # of projected period if aggregated by nagg
    
    nycp <- nperp*nagg # of projection years be aggregated
    ry <- nyp - nycp # of rest projection years not aggregated

    pdat1 <- pdat0[ , 1:ny] # yearly observed population
    pdat2 <- pdat0[ , (ny+1):nypop] # yearly projected population

    mt <- dim(rate)[2] # of periods from output
    r0 <- as.matrix(rate[, (nper+1):mt]) # of projected periods 
    m <- dim(r0)[2] 
# producing the end periods rates:
    rc1 <- rate[, nper]
    if (m > 1) {
       rc2 <- 2*r0[, m] - r0[, (m-1)]
    } else {
       rc2 <- 2*r0[, m] - rc1
    }
    rc3 <- 2*rc2 - r0[, m]
    r1 <- cbind(rc1, r0, rc2)
# producing annual age-specific rates:
  rr <- matrix(NA, 19, nyp)
  if (nagg==1) {
    rr <- r0
  } else if (nagg==2) {
    for (i in 2:(m+1)) {
         rr[,(i-2)*2+1] <- (1/4)*r1[,i-1] + (3/4)*r1[,i]
         rr[,(i-2)*2+2] <- (1/4)*r1[,i+1] + (3/4)*r1[,i]
    }
    if (nyp > nycp) {
        rr[, nyp] <- (1/4)*r1[,(m+1)] + (3/4)*r1[,(m+2)]
    }
  } else if (nagg==3) {
     for (i in 2:(m+1)) {
         rr[,(i-2)*3+1] <- (1/3)*r1[,i-1] + (2/3)*r1[,i]
         rr[,(i-2)*3+2] <- (0/3)*r1[,i-1] + (3/3)*r1[,i]
         rr[,(i-2)*3+3] <- (1/3)*r1[,i+1] + (2/3)*r1[,i]
    }
    if (nyp == (nycp+1)) {
         rr[, nyp] <- (1/3)*r1[,(m+1)] + (2/3)*r1[,(m+2)]
    } 
    if (nyp == (nycp+2)) {
         rr[, (nyp-1)] <- (1/3)*r1[,(m+1)] + (2/3)*r1[,(m+2)]
         rr[, nyp] <- rc2
    }   
  } else if (nagg==4) {
     for (i in 2:(m+1)) {
         rr[,(i-2)*4+1] <- (3/8)*r1[,i-1] + (5/8)*r1[,i]
         rr[,(i-2)*4+2] <- (1/8)*r1[,i-1] + (7/8)*r1[,i]
         rr[,(i-2)*4+3] <- (1/8)*r1[,i+1] + (7/8)*r1[,i]
         rr[,(i-2)*4+4] <- (3/8)*r1[,i+1] + (5/8)*r1[,i]
    }
    if (nyp == (nycp+1)) {
         rr[, nyp] <- (3/8)*r1[,(m+1)] + (5/8)*r1[,(m+2)]
    } 
    if (nyp == (nycp+2)) {
         rr[, (nyp-1)] <- (3/8)*r1[,(m+1)] + (5/8)*r1[,(m+2)]
         rr[, nyp] <- (1/8)*r1[,(m+1)] + (7/8)*r1[,(m+2)]
    }   
    if (nyp == (nycp+3)) {
         rr[, (nyp-2)] <- (3/8)*r1[,(m+1)] + (5/8)*r1[,(m+2)]
         rr[, (nyp-1)] <- (1/8)*r1[,(m+1)] + (7/8)*r1[,(m+2)]
         rr[, nyp] <- (7/8)*r1[,(m+2)] + (1/8)*rc3
    }   
  } else if (nagg==5) {
    for (i in 2:(m+1)) {
         rr[,(i-2)*5+1] <- (2/5)*r1[,i-1] + (3/5)*r1[,i]
         rr[,(i-2)*5+2] <- (1/5)*r1[,i-1] + (4/5)*r1[,i]
         rr[,(i-2)*5+3] <- (0/5)*r1[,i-1] + (5/5)*r1[,i]
         rr[,(i-2)*5+4] <- (1/5)*r1[,i+1] + (4/5)*r1[,i]
         rr[,(i-2)*5+5] <- (2/5)*r1[,i+1] + (3/5)*r1[,i]
    }
    if (nyp == (nycp+1)) {
         rr[, nyp] <- (2/5)*r1[,(m+1)] + (3/5)*r1[,(m+2)]
    } 
    if (nyp == (nycp+2)) {
         rr[, (nyp-1)] <- (2/5)*r1[,(m+1)] + (3/5)*r1[,(m+2)]
         rr[, nyp] <- (1/5)*r1[,(m+1)] + (4/5)*r1[,(m+2)]
    }   
    if (nyp == (nycp+3)) {
         rr[, (nyp-2)] <- (2/5)*r1[,(m+1)] + (3/5)*r1[,(m+2)]
         rr[, (nyp-1)] <- (1/5)*r1[,(m+1)] + (4/5)*r1[,(m+2)]
         rr[, nyp] <- rc2
    }   
    if (nyp == (nycp+4)) {
         rr[, (nyp-3)] <- (2/5)*r1[,(m+1)] + (3/5)*r1[,(m+2)]
         rr[, (nyp-2)] <- (1/5)*r1[,(m+1)] + (4/5)*r1[,(m+2)]
         rr[, (nyp-1)] <- rc2
         rr[, nyp] <- (1/5)*rc3 + (4/5)*rc2
    }   
  } else {
    stop("Years of aggregation \"nagg\" must be integer between 1 and 5")
  }

# fill in observed age-specific rates per 100,000:
    obasr <- matrix(NA, 19, nt)
    for (i in 1:nt) {
        obasr[, i] <- 100000*cdat[, i]/pdat[, i]
    }
    datatab <- matrix(NA, 19, np)
    datatab[ , 1:nt] <- as.matrix(obasr)
    datatab <- data.frame(datatab)
    row.names(datatab) <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34",
                            "35-39","40-44","45-49","50-54","55-59","60-64",
                            "65-69","70-74","75-79","80-84","85-89","90+")
    colnames(datatab) <- (startp-nt):(startp+nyp-1)
# fill in projected age-specific rates per 100,000:
    datatab[ , (nt+1):np] <- rr
    return(datatab)
}


## Convert to annual ASRs and Counts:
#asrny <- function(r0, cdat, pdat, startp, nagg, standpop=ca91) {
#    rr <- asrpy(r0, cdat, pdat, startp=startp, nagg)
#    c1 <- rr * pdat / 100000
#    a1 <- rr * standpop
#    asr <- round(apply(a1, 2, sum), 6)
#    case <- round(apply(c1, 2, sum), 0)
#    for (i in 1: length(asr)) {if (asr[i] < 0) asr[i] <- 0}
#    for (i in 1: length(case)) {if (case[i] < 0) case[i] <- 0}
#    return(cbind(asr, case))
#}

## Or
asry <- function(rr, pdat, standpop=ca11) {
    c1 <- rr * pdat / 100000
    a1 <- rr * standpop
    asr <- round(apply(a1, 2, sum), 6)
    case <- round(apply(c1, 2, sum), 0)
    for (i in 1: length(asr)) {if (asr[i] < 0) asr[i] <- 0}
    for (i in 1: length(case)) {if (case[i] < 0) case[i] <- 0}
    return(cbind(asr, case))
}

## Standardizing observed rates and total numbers:
obasr <- function(cdat, pdat, standpop=ca11){
    m <- dim(cdat)[2]
    popu <- pdat[ , 1:m]
    rr <- matrix(NA, 19, m)
    for (i in 1:m) {
        rr[, i] <- 100000*cdat[, i]/popu[, i]
    }
    aspr <- rr * standpop
    asr <- apply(aspr, 2, sum)
    case <- apply(cdat, 2, sum)
    return(cbind(asr, case))
}

asrsd <- function(cdat, pdat, standpop=ca11){
## calculate the age-standardized rate with correspond standard error
## cdat, pdat: 19 age groups by N years
    m <- dim(cdat)[2]
    popu <- pdat[ , 1:m]
    rr <- matrix(NA, 19, m)
    ww <- matrix(NA, 19, m)
    for (i in 1:m) {
        rr[, i] <- cdat[, i]/popu[, i]
        ww[, i] <- cdat[, i]*(standpop/popu[, i])^2
    }
    aspr <- rr*standpop
    asr <- 100000*apply(aspr, 2, sum)
    asd <- 100000*sqrt(apply(ww, 2, sum)) 
    case <- apply(cdat, 2, sum)
    return(cbind(asr, asd))
}


chper <- function(aspr, pdat, byear, cyear, starty=NULL) {
## Calculating the percentage changes due to risk, population growth and aging.
## pdat: population data over historical and projection years,
## aspr: data frame for annual age-specific rates,
## byear: baseline calendar year, cyear: comparison calendar year,
## starty: the first calendar year in historical data.

# Check data:
    if ( dim(aspr)[2]!= dim(pdat)[2]) {
      stop("\"aspr\" and \"pdat\" must have data for same calendar years")	
    }
# Identify the first calendar year:
    if (is.null(starty)) {
         caly <- as.numeric(names(aspr))
         starty <- caly[1]
    }
# Identify the numbers of columns for baseline and comparison:
    n1 <- byear - starty + 1
    n2 <- cyear - starty + 1
    asr12 <- aspr[, c(n1, n2)]
    pat12 <- pdat[, c(n1, n2)]
# Total numbers in baseline and comparison:
    N1 <- sum(pat12[,1]*asr12[,1]/100000)
    N2 <- sum(pat12[,2]*asr12[,2]/100000)
# Crude rate in baseline year:
    R1 <- N1/sum(pat12[,1])
# Total number in comparison if no change in risk:
    N3 <- sum(pat12[,2]*asr12[,1]/100000)
# Total number in comparison if no change by aging:
    N4 <- R1 * sum(pat12[,2])

# Overall change:
    C.all <- 100 * (N2 - N1)/N1
# Change due to population growth and aging:
    C.pop <- 100 * (N3 - N1)/N1
# Change due to risk:
    C.rsk <- C.all - C.pop
# Change due to population growth:
    C.gwh <- 100 * (N4 - N1)/N1
# Change due to population aging:
    C.age <- C.pop - C.gwh

# Output:
    out <- data.frame(matrix(c(N1, N2, C.all, C.rsk, C.gwh, C.age), ncol=6))
    colnames(out) <- c("ref.case", "comp.case", "overall", "risk", "p.growth", "p.aging")
    return(out)
}


# Projection plot contains both ASRs and total numbers:
projplot <- function(site.asr, sex=NULL, ma=2, mr=1.02, starty=1986, startp=2011) {
    case <- site.asr[,2]
    asr <- site.asr[,1]
    n <- dim(site.asr)[1] #total number of years
    n1 <- startp - starty #n1 in observations
    n2 <- n - n1 #n2 in projections
    maxn <- max(case)
    maxr <- max(asr)
    if (is.null(sex)) {
        mycol <- c("red4", "red2", "red4", "red2")
    } else if (sex=="M") {
        mycol <- c("darkblue", "blue", "darkblue", "blue2")
    } else {
        mycol <- c("hotpink", "lightpink", "hotpink", "lightpink2")
    }
    par(mar=c(4,4,2,4),mgp=c(3,1,0), cex=0.7)  
    barplot(c(case[1:n1],rep(0,n2)),col=mycol[1],space=1,ylim=c(0,ma*maxn),
#            cex.axis=1,xlab="", las=1, xaxs="r", border = NA, axis.lty = 1)
            cex.axis=1,xlab="", las=1, xaxs="r", border = NA)
    abline(h=0)
    par(mgp=c(3,1,0))
    barplot(c(rep(0,n1),case[(n1+1):n]),add=T,col=mycol[2],space=1,
             cex.axis=1,las=1,xaxs="r",border = NA)
    par(new=T, mar=c(4,4,2,4),mgp=c(3,1,0), cex=1)  
    plot(c(starty:(starty+n-1)),c(asr[1:n1],rep("",n2)),type="o",pch=20, ylim=c(0,maxr),
         xaxt="n",yaxt="n",bty="n",lwd=3,ylab="",col=mycol[3],xlab="",cex.axis=1)
    lines(c(starty:(starty+n-1)),c(rep("",(n1-1)),asr[n1:n]),type="l",pch=20, ylim=c(0,maxr),
          xaxt="n",yaxt="n",bty="n",lwd=3,lty=2,col=mycol[4], cex.axis=1)
    axis(4, at=seq(0,(mr*maxr),ceiling(maxr/6)),las=1, cex.axis=1)
#    axis(4, las=1, cex.axis=1)
    mtext(side=4,"", line=2.5, cex=1)
}


## Standard populations:
# ca91 <- c(0.069465, 0.069454, 0.068034, 0.068495, 0.075016, 0.089944, 0.092400, 0.083388, 
#           0.076063, 0.059536, 0.047649, 0.044041, 0.042326, 0.038570, 0.029660, 0.022127, 
#           0.013595, 0.010237)
# 
# ca96 <- c(0.066235, 0.067985, 0.067716, 0.067841, 0.067761, 0.072914, 0.087030, 0.088510, 
#           0.080055, 0.071847, 0.055812, 0.044869, 0.040705, 0.037858, 0.032589, 0.023232, 
#           0.015424, 0.011617)

#ca11 <- c(0.055297, 0.052717, 0.055853, 0.065194, 0.068555, 0.069006, 0.067786, 0.066188,
#           0.069474, 0.079199, 0.078365, 0.068518, 0.059705, 0.044636, 0.033597, 0.026769,
#           0.020416, 0.018725)
#Updated June 2019
ca11 <- c(0.055325, 0.052736, 0.055857, 0.065121, 0.068506, 0.068967, 0.067743, 0.066176,
          0.069477, 0.079209, 0.078366, 0.068519, 0.059696, 0.044613, 0.033573, 0.026761,
          0.020406, 0.012420, 0.006529)

# Segi world population
# wdsd <- c(.12, .1, .09, .09, .08, .08, .06, .06, .06, .06, .05, .04, .04, .03, .02, .01, 
#           .005, .005)

# New WHO world population
wdsd <- c(0.08855, 0.08687, 0.08597, 0.08474, 0.08222, 0.07928, 0.07605, 0.07145, 0.0659, 0.06038,  
          0.05371, 0.04547, 0.03723, 0.02955, 0.0221,  0.01515, 0.00905, 0.0044,  0.00193)

