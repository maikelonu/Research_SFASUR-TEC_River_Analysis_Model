# Streamflow and Flood Analysis Using R (SFAUR-TEC)
# Instituto Tecnologico de Costa Rica (www.tec.ac.cr)
# Maikel Mendez-M (mamendez@itcr.ac.cr);(maikel.mendez@gmail.com)
# Luis Alexander Calvo-V (lcalvo@itcr.ac.cr);(lualcava.sa@gmail.com)
# This script is structured in R (www.r-project.org)
# General purpose: Generate various graphical and numerical products 
# for streamflow and flood analysis at river catchments
# Custom functions: FrequencyAnalysis() and BootstrapCI() created by
# Dave Hutchinson (mtb_dave[at]yahoo[dot]ca)
# Script conceptualization: based on SAAS (Streamflow Analysis and Assessment Software)
# by Robert A. Metcalfe (http://people.trentu.ca/~rmetcalfe/SAAS.html)
# Input files: 
# Output files:

# Workspace is cleared
rm(list = ls())

# working directory is defined
# setwd("B:\\R_ITC\\SAAS_diario")
setwd("B:\\R_ITC\\SAAS_diario")

# CRAN libraries are loaded
require(dplyr)
require(EcoHydRology)
require(evd)
require(extRemes)
require(ggplot2)
require(ggthemes)
require(hydroTSM)
require(lmom)
require(lmomco)
require(lubridate)
require(pastecs)
require(reshape)
require(reshape2)
require(scales)
require(tidyr)
require(xts)
require(zoo)

# ////////////////////////
# BLOCK: Custom Functions
# ////////////////////////

# Custom Function: FrequencyAnalysis
# Fits a given extreme value distribution to an extreme value series
# @param series A vector representing an extreme value series (e.g., annual maximum flood)
# @param distribution A three-character name of the extreme value distribution (see ?dist.list())
# @param nep A vector of non-exceedance probabilities
# @return A list object containing: (1) distribution information and (2) output
# (quantile estimates at various non-exceedance probabilities)
# @export
# @import lmomco

FrequencyAnalysis <- function( series, distribution, nep = nonexceeds() ) {
  
  distribution <- tolower(distribution)
  transformed <- FALSE
  
  # add log Pearson Type 3 to list of distributions supported
  # by lmomco package
  base.dist <- c('lp3', dist.list())
  
  if( any(distribution %in% base.dist) ) {
    
    # log transform series 
    if( distribution == 'lp3' ) {
      series <- log10(series)
      transformed <- TRUE
      distribution <- 'pe3'
    }
    
    # compute L-moments
    samLmom <- lmom.ub(series)
    
    # estimate distribution parameters
    distPar <- lmom2par(samLmom, type = distribution)
    
    # compute quantiles for nonexceedances
    quant <- par2qua(f = nep, para = distPar)
    
    if( distribution == 'pe3' & transformed ) {
      distribution <- 'lp3'
      quant <- 10^quant
    }
    
    # return result as list object
    return(
      list(
        distribution = list(
          name = distribution,
          logTransformed = transformed,
          parameters = distPar),
        output = data.frame(nep = nep, rp = prob2T(nep), estimate = quant) 
      ) )
    
  } else {
    stop(
      sprintf('Distribution \'%s\' not recognized!', distribution))
  }
}

# Custom Function: BootstrapCI
# Conducts bootstrap to randomly sample an extreme value series 'n' times for a 
# specified distribution to estimate confidence interval for each given 
# non-exceedance probability.
# @param fitted.model Fitted distribution (see ?frequencyAnalysis)
# @param series A vector representing an extreme value series (e.g., annual maximum flood)
# @param distribution A three-character name of the extreme value distribution (see ?dist.list())
# @param n.resamples An integer representing number of re-samples to conduct
# @param nep A vector of non-exceedance probabilities
# @param ci The confidence interval 
# @export
# @import lmomco
# @return A list containing a data frame of confidence bounds for quantile estimates for each 
# non-exceedance probability, a matrix containing estimated distribution parameters for each resample,
# and a matrix of quantile estimates for each resample

BootstrapCI <- function(series, distribution, n.resamples=1E3, nep=nonexceeds(), ci=0.90) {
  
  # compute frequency analysis
  fa <- FrequencyAnalysis(series=series, distribution=distribution, nep=nep)
  
  # extract fitted model parameters and flag as to whether the 
  # distribution is based on log transformed data
  base.params <- fa$distribution$parameters
  isTransformed <- fa$distribution$logTransformed
  
  # create output matrices to store parameter sets and quantile estimates
  param.sets <- matrix(NA, nrow = n.resamples, ncol = length(base.params$para))
  quantile.estimates <- matrix(NA, nrow = n.resamples, ncol = length(nep), 
                               dimnames = list(NULL, nep) ) 
  
  # begin bootstrapping procedure
  for(i in 1:n.resamples) {
    
    valid.moments <- FALSE
    j <- 0
    
    # allow up to 20 re-tries to re-sample 
    while(!valid.moments & j < 20) {  
      
      # sample 'n' random variates from base distribution
      data <- rlmomco(n=length(series), base.params)
      
      # compute sample l-moments
      sample.moms = lmom.ub(data)
      
      valid.moments <- are.lmom.valid(sample.moms)
      j <- j + 1
    }
    
    # error handling
    if(!valid.moments) {
      stop("Bootstrapping failed to sample valid l-moments")
    } else {
      # estimate distribution parameters
      dist.par <- lmom2par(sample.moms, base.params$type)
      
      # store the distribution parameters
      param.sets[i,] <- dist.par$para
      
      # estimate quantiles at NEP
      estimated <- qlmomco(nep, dist.par)
      
      # convert quantile estimates to real values if
      # distribution was transformed
      if(isTransformed) estimated <- 10^estimated
      
      # store the quantiles at the desired AEP values
      quantile.estimates[i,] <- estimated
    } 
    
  }
  
  # now calculate confidence limits for quantiles
  p <- c((1-ci)/2, (1+ci)/2)
  ci <- sapply(colnames(quantile.estimates), 
               FUN=function(x){
                 quantile(quantile.estimates[,x], probs=p, na.rm=TRUE)})
  
  # now return list object containing output
  return(
    list(
      ci = data.frame(
        nonexceed_prob = nep,
        lower = as.vector(ci[1,]),
        true = fa$output$estimate,
        upper = as.vector(ci[2,]) ),
      parameters = param.sets,
      quantiles = quantile.estimates)
  )
  
}

boundary.LM <- function(max.X,min.X,max.Y,min.Y) {
# Defines boundaries of labels-position for a LM model-plot using ggplot2
#
# Args:  
#  max.X: horizontal-axis maximum extension
#  min.X: horizontal-axis minumum extension
#  max.Y: vertical-axis maximum extension
#  min.Y: vertical-axis minumum extension
#
# Returns:
#  a vector containing label X and Y position in ggplot2-plot
  
delta.X <- (max.X - min.X)
  pos.X <- (min.X + (delta.X)*0.05)  
delta.Y <- (max.Y - min.Y)
  pos.Y <- (min.Y + (delta.Y)*1.10)
return(c(pos.X,pos.Y))
}

eq.PAR <- function(lm.model) {
# Creates a label parameter-summary for LM models in ggplot2
#
# Args:  
#  lm.model: a lm model previously defined
#
# Returns:
#  a character-pasted vector to be displayed in ggplot2

paste("Adj R2 = ",round(summary(lm.model)$adj.r.squared,4),
      "; Intercept =",round(lm.model$coef[[1]],4),
      "; Slope =",round(lm.model$coef[[2]],4),
      "; p-value =",round(summary(lm.model)$coef[2,4],4))
}

# /////////////////////////////////////////////////
# BLOCK: Creating and organizing input data.frames
# /////////////////////////////////////////////////

# Observations data.frame is loaded
df.obs <- read.table("inputflow.txt",header=T,sep="\t",quote="")

# "DATE" class character is converted to class-date and added as a new column named "DATE2"
        temp <- df.obs$DATE
df.obs$DATE2 <- as.Date(temp, format = "%d/%m/%Y")

# A date-class query is requested (TRUE or FALSE)
is.Date(df.obs$DATE2)

# lubridate Library functions are applied to df.obs to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH, WEEK, DAY, DAY_MONTH   
     df.obs$YEAR <- year(df.obs$DATE2) # Years component of a date-time
  df.obs$YEAR_CH <- as.character(year(df.obs$DATE2)) # Years component of a date-time as character
    df.obs$MONTH <- month(df.obs$DATE2, label = FALSE) # Months component of a date-time
 df.obs$MONTH_CH <- month(df.obs$DATE2, label = TRUE) # Months component of a date-time as character
     df.obs$WEEK <- week(df.obs$DATE2) # Weeks component of a date-time
      df.obs$DAY <- yday(df.obs$DATE2) # Days component of a date-time
df.obs$DAY_MONTH <- days_in_month(df.obs$DATE2) # Number of days in the month of a date-time

# Descriptive statistics are extracted for df.obs data.frame for "FLOW" column only
df.obs.desc <- round((as.data.frame(stat.desc(df.obs$FLOW))),3)

# colnames in df.obs.desc data.frame are renamed
colnames(df.obs.desc) <- c("FLOW")

# Total length of df.obs data.frame is requested
obs.length <- length(df.obs$FLOW)

# A subset data.frame containing only "YEAR","DAY" and "FLOW" is created      
df.obs.sub <- df.obs[,c("YEAR","DAY","FLOW")]

# A pivot data.frame is created organizing data by "YEAR" and "FLOW"
df.pivot <- df.obs.sub %>% spread(YEAR,FLOW)

# Apply function is applied to df.pivot to generate mean, median and sd at row level
  df.pivot$MEAN <- round(apply(df.pivot [ , 2:18],1,mean),3)
df.pivot$MEDIAN <- round(apply(df.pivot [ , 2:18],1,median),3)
    df.pivot$SD <- round(apply(df.pivot [ , 2:18],1,sd),3)

# NA values are omitted from df.pivot (leap years)
df.pivot <- na.omit(df.pivot)

# df.pivot data.frame is ordered (sorted) by "DAY" value  
df.pivot <- df.pivot[order(df.pivot[ ,1]) , ]

# A prefix "YEAR" is added to df.pivot column names
colnames(df.pivot) <- paste("YEAR", colnames(df.pivot), sep = "_")

# Descriptive statistics are extracted from df.pivot data.frame
df.pivot.desc <- round(stat.desc (df.pivot[ , 2:21]),3)

# A temporal date-class is re-incorporated in df.pivot data.frame
df.pivot$tempDATE <- as.Date(df.pivot$YEAR_DAY - 1, origin = "2015-01-01")

# lubridate Library functions are applied to create new columns to df.pivot data.frame
# contaning: MONTH, MONTH_CH, and DAY_YEAR, DAY_MONTH   
   df.pivot$MONTH <- month(df.pivot$tempDATE, label = FALSE) # Months component of a date-time
df.pivot$MONTH_CH <- month(df.pivot$tempDATE, label = TRUE) # Months component of a date-time as character
df.pivot$DAY_YEAR <- yday(df.pivot$tempDATE) # Days component of a date-time

# ////////////////////////////////
# BLOCK: Flood Frequency Analysis
# ////////////////////////////////

# CRAN library data.table is loaded at this point to avoid conflicts with other libraries
require(data.table)

# A subset data.frame is created base on df.obs data.frame
df.obs.sub2 <- data.table(date = as.IDate(df.obs$DATE2), df.obs[-1])

# Maximum annual flow, along with other annual statistics are extracted 
# from df.obs.sub2 into a new data.frame called df.annual.flow
df.annual.flow <- as.data.frame(df.obs.sub2[, list(mean.FLOW = mean(FLOW),
                                                   median.FLOW = median(FLOW),
                                                   min.FLOW = min(FLOW),
                                                   max.FLOW = max(FLOW)), 
                                                   by = year(date)])

# df.annual.flow data.frame is rounded to three significant digits
df.annual.flow <- round(df.annual.flow,3)

# CRAN library data.table is detached to avoid conflicts with other libraries
detach(package:data.table)

# Maximum flow dataset is extracted df.annual.flow data.frame
# and stated as input_flow (numeric vector)
input_flow <- df.annual.flow$max.FLOW

# A Distribution function is selected
# (see ?dist.list for available distributions)
# log Pearson 3 is not one of them, but is a log-transformed equivalent of pe3
# note: this script recognizes 'lp3' to stand for log Pearson Type 3
dist <- "lp3" # this script uses "gev" by default; not "lp3"

# Frequency distribution is fitted as function of "dist"
# and parameters "mu", "sigma" and "gamma" are calculated
fa <- FrequencyAnalysis(series=input_flow, distribution=dist)

# A data.frame containing estimated-deterministic fitted values, 
# non-exceedance probabilities (nep) and return period (rp) is created
df.fa.out <- fa$output

# 95% confidence intervals are estimated for chosen frequency-distribution
ci <- BootstrapCI(series=input_flow,   # flow data
                  distribution=dist,   # distribution
                  n.resamples = 2.5E4, # number of re-samples to conduct
                  ci = 0.95)           # confidence interval level

# A data.frame containing non-exceedance probabilities (nep) for "central", "lower"
# and "upper" CI is created
df.ci.out <- ci$ci

# Maximum flow values are sorted df.annual.flow data.frame
df.annual.flow <- df.annual.flow[order(df.annual.flow$max.FLOW),]

# length of max.FLOW per year is requested from df.annual.flow data.frame 
n.length <- length(df.annual.flow$max.FLOW)

# Weibull probabilities plotting-positions are calculated and added to df.annual.flow data.frame
df.annual.flow$PROB <- round(((1:n.length)/(1+n.length)),3) 

# Frequency analysisis plot is configured by defining: specific x-breaks, ceiling and floor
# range and scale logarithmic transformation
  bwpeaks <- data.frame(PROB = df.annual.flow$PROB , FLOW = df.annual.flow$max.FLOW)
  xbreaks <- c(0.002,0.01,0.1,0.3335,0.5,0.8,0.9,0.95,0.975,0.99,0.995, 0.998) # Plotting positions are determined
log.range <- log10(range(bwpeaks$FLOW))
    lower <- 10^floor(log.range[1])
    upper <- 10^ceiling(log.range[2])
      cap <- lower
  ybreaks <- NULL
while(cap < upper) {
  ybreaks <- c(ybreaks, seq(cap, cap*15, by = cap))
      cap <- cap * 10
  }

# Flood Frequency Analysisis plot is generated and saved
g.ffa <- ggplot(bwpeaks) + 
         geom_line(data=df.ci.out, aes(x=nonexceed_prob, y=true),size = 1.25, color="#cc0000", alpha = 0.50) +
         geom_line(data=df.ci.out, aes(x=nonexceed_prob, y=lower),size = 0.75, color="#333333",alpha = 0.75, lty=2) +
         geom_line(data=df.ci.out, aes(x=nonexceed_prob, y=upper),size = 0.75, color="#333333",alpha = 0.75, lty=2) +
         geom_point(aes(x=PROB, y=FLOW),size = 4.00) + 
         scale_y_continuous(trans="log10", breaks=ybreaks) +
         scale_x_continuous(trans=probability_trans(distribution="norm"),
                            breaks=xbreaks, labels=signif(prob2T(xbreaks), digits=4),
                            name="Return period [yrs]") +
         ylab(label = 'Q (m3/s)') +
         xlab(label = 'Recurrence Interval (years)') +
         ggtitle(label = 'Upper Toro River Catchment. Flood Frequency Analysisis (1994-2010). Log Pearson Type 3') +
         theme_bw(base_size = 18.0)

# Flood Frequency Analysisis plot is requested
g.ffa

# Flow events with recurrence interval = 1.5 years are identified and a new subset data.frame is created 
# Their percentage is also calculated
        limit1.5 <- df.fa.out[12,3] 
        FFATr1.5 <- sum(df.obs$FLOW > limit1.5)
 df.FFATr1.5_sel <- df.obs [df.obs$FLOW > limit1.5, ]
FFATr1.5_percent <- round(((FFATr1.5 / obs.length)*100),3)

# Flow events with recurrence interval = 10 years are identified and a new subset data.frame is created 
# Their percentage is also calculated
        limit10 <- df.fa.out[23,3] 
        FFATr10 <- sum(df.obs$FLOW > limit10)
 df.FFATr10_sel <- df.obs [df.obs$FLOW > limit10, ]
FFATr10_percent <- round(((FFATr10 / obs.length)*100),3)

# A complete-series hydrograph is generated and saved
# Data points over FFATr1.5 and FFATr10 are shown above threshold lines
g.hydro01 <- ggplot() +
             geom_point(aes(x = DATE2,y = FLOW,colour = YEAR_CH),data=df.obs,size = 1.50) +
             geom_point(aes(x = DATE2,y = FLOW),data=df.FFATr1.5_sel,shape = 21,colour = '#999900',size = 3.0) +
             geom_line(aes(x = DATE2,y = FLOW),data=df.obs,colour = '#666666',size = 0.15,alpha = 0.5126) +
             scale_y_continuous(breaks = scales::pretty_breaks(n = 8.0,min.n = 1.0),expand = c(0.05,0.50),limits = c(0,70)) +
             scale_x_date(breaks = scales::pretty_breaks(n = 6.0)) +
             geom_hline(data=df.obs,size = 0.7,alpha = 0.5,yintercept = (df.obs.desc[8,1])) +
             geom_hline(data=df.obs,size = 0.7,linetype = 2,alpha = 0.5,yintercept = (df.obs.desc[9,1]), color="black") +
             geom_hline(data=df.obs,colour = '#999900',yintercept = limit1.5) +
             geom_hline(data=df.obs,colour = '#999900',linetype = 2,yintercept = limit10) +
             ylab(label = 'Q (m3/s)') +
             xlab(label = 'Period (years)') +
             ggtitle(label = 'Upper Toro River Catchment Streamflow Analisis (1994-2010)') +
             geom_text(aes(df.obs[1,4],limit1.5,label = "     R.int > 1.5 years", vjust = -1)) +
             geom_text(aes(df.obs[1,4],limit10,label = "      R.int > 10 years", vjust = -1)) +
             geom_text(aes(df.obs[1,4],df.obs.desc[8,1],label = "median", vjust = 1.5)) +
             geom_text(aes(df.obs[1,4],df.obs.desc[9,1],label = "mean", vjust = -1)) +
             theme_bw(base_size = 18.0)

# Complete-series hydrograph is requested
g.hydro01

# Only columns containing "YEAR_" are selected from df.pivot data.frame 
i1 <- grep("YEAR_", names(df.pivot))

# A df.L data.frame is created
df.L <- (df.pivot)[i1]

# Irrelevant columns are erased from df.L data.frame 
   df.L$YEAR_DAY <- NULL
  df.L$YEAR_MEAN <- NULL
df.L$YEAR_MEDIAN <- NULL
    df.L$YEAR_SD <- NULL

# "tempDATE" is incorporated in df.L data.frame  
df.L$tempDATE <- df.pivot$tempDATE

# melt function from reshape2 CRAN-Library  is applied  
# and a new data.frame is created for graphical purposes
df.L2 <- melt(df.L,id.vars="tempDATE")

# An annual-summary hydrograph is generated. The mean and the median are highlighted
# Both df.L2 and df.pivot are called  
g.hydro02 <- ggplot() +
             geom_line(aes(x = tempDATE,y = value, group=variable),data=df.L2,size = 0.1,alpha = 0.20) +
             geom_line(aes(x = tempDATE,y = YEAR_MEAN),data=df.pivot,colour = '#ff3300',size = 0.75) +
             geom_line(aes(x = tempDATE,y = YEAR_MEDIAN),data=df.pivot,colour = '#0000cc',size = 0.75) +   
             scale_y_continuous(breaks = scales::pretty_breaks(n = 8.0,min.n = 8.0)) +
             scale_x_date(breaks = scales::date_breaks(),labels = date_format(format = '%m')) +
             ylab(label = 'Q (m3/s)') +
             xlab(label = 'Period (months)') +
             ggtitle(label = 'Upper Toro River Catchment Streamflow Analisis (1994-2010)') +
             geom_text(aes(df.pivot[1,22],5,label = "mean", vjust = -1)) +  
             geom_text(aes(df.pivot[1,22],3,label = "median", vjust = 1)) +             
             theme_bw(base_size = 18.0)

# annual-summary hydrograph is requested
g.hydro02

# An annual-summary hydrograph is generated (colored by "YEAR"). The mean and the median are highlighted
# Both df.L2 and df.pivot are called
g.hydro03 <- ggplot() +
             geom_line(aes(x = tempDATE,y = value, group=variable, colour=variable),data=df.L2,size = 0.1,alpha = 0.55) +
             geom_line(aes(x = tempDATE,y = YEAR_MEAN),data=df.pivot,colour = '#ff3300',size = 0.75) +
             geom_line(aes(x = tempDATE,y = YEAR_MEDIAN),data=df.pivot,colour = '#0000cc',size = 0.75) +   
             scale_y_continuous(breaks = scales::pretty_breaks(n = 8.0,min.n = 8.0)) +
             scale_x_date(breaks = scales::date_breaks(),labels = date_format(format = '%m')) +
             ylab(label = 'Q (m3/s)') +
             xlab(label = 'Period (months)') +
             ggtitle(label = 'Upper Toro River Catchment Streamflow Analisis (1994-2010)') +
             geom_text(aes(df.pivot[1,22],5,label = "mean", vjust = -1)) +  
             geom_text(aes(df.pivot[1,22],3,label = "median", vjust = 1)) +             
             theme_grey()

# annual-summary hydrograph is requested
g.hydro03

# ////////////////////////////
# BLOCK: Flow Duration Curves
# ////////////////////////////

# A subset df.obs.month data.frame containing only "MONTH_CH" and "FLOW" is created      
# based on df.obs data.frame
df.obs.month <- df.obs[,c("MONTH_CH","FLOW")]

# df.obs.month data.frame is ordered (sorted) by "MONTH_CH"" value for the whole period  
df.obs.month <- df.obs.month[order(df.obs.month[ ,1]) , ]

# A rep function is introduced
df.obs.month$ID <-rep (1:length(df.obs.month$FLOW),1)

# dcast function is requested to convert from long to wide format
# (unmelt) NAs are introduced by coercion !!!!!!!!!!!!!!!!!
df.obs.month <- dcast(df.obs.month, ID ~ MONTH_CH, value.var="FLOW")

# ID variable is erased from df.obs.month
df.obs.month$ID <- NULL

# Descriptive statistics are extracted for df.obs.month data.frame 
# and rounded to three significant digits
df.obs.month.desc <- round((as.data.frame(stat.desc(df.obs.month))),3)

# Rows 2 and 3 are deleted from df.obs.month.desc data.frame to avoid confusion
df.obs.month.desc <- df.obs.month.desc[-c(2, 3), ] 

# hydroTSM fdc function is applied and plotting is suppressed
df.FDC <- fdc(df.obs.month, plot=FALSE)

# df.FDC matrix is transformed into data.frame class 
df.FDC <- as.data.frame(df.FDC)

# A data.frame is created for each month
df.Jan <- data.frame(x=(df.FDC$Jan)*100, y=df.obs.month$Jan)
df.Feb <- data.frame(x=(df.FDC$Feb)*100, y=df.obs.month$Feb)
df.Mar <- data.frame(x=(df.FDC$Mar)*100, y=df.obs.month$Mar)
df.Apr <- data.frame(x=(df.FDC$Apr)*100, y=df.obs.month$Apr)
df.May <- data.frame(x=(df.FDC$May)*100, y=df.obs.month$May)
df.Jun <- data.frame(x=(df.FDC$Jun)*100, y=df.obs.month$Jun)
df.Jul <- data.frame(x=(df.FDC$Jul)*100, y=df.obs.month$Jul)
df.Aug <- data.frame(x=(df.FDC$Aug)*100, y=df.obs.month$Aug)
df.Sep <- data.frame(x=(df.FDC$Sep)*100, y=df.obs.month$Sep)               
df.Oct <- data.frame(x=(df.FDC$Oct)*100, y=df.obs.month$Oct)
df.Nov <- data.frame(x=(df.FDC$Nov)*100, y=df.obs.month$Nov)
df.Dec <- data.frame(x=(df.FDC$Dec)*100, y=df.obs.month$Dec)

# NA values are omitted from each data.frame
df.Jan <- na.omit(df.Jan)
df.Feb <- na.omit(df.Feb)
df.Mar <- na.omit(df.Mar)
df.Apr <- na.omit(df.Apr)
df.May <- na.omit(df.May)
df.Jun <- na.omit(df.Jun)
df.Jul <- na.omit(df.Jul)
df.Aug <- na.omit(df.Aug)
df.Sep <- na.omit(df.Sep)
df.Oct <- na.omit(df.Oct)
df.Nov <- na.omit(df.Nov)
df.Dec <- na.omit(df.Dec)

# data.frames are ordered (sorted) by "% Excedeed" value
df.Jan <- df.Jan[order(df.Jan[ ,1]) , ]
df.Feb <- df.Feb[order(df.Feb[ ,1]) , ]
df.Mar <- df.Mar[order(df.Mar[ ,1]) , ]
df.Apr <- df.Apr[order(df.Apr[ ,1]) , ]
df.May <- df.May[order(df.May[ ,1]) , ]
df.Jun <- df.Jun[order(df.Jun[ ,1]) , ]
df.Jul <- df.Jul[order(df.Jul[ ,1]) , ]
df.Aug <- df.Aug[order(df.Aug[ ,1]) , ]
df.Sep <- df.Sep[order(df.Sep[ ,1]) , ]
df.Oct <- df.Oct[order(df.Oct[ ,1]) , ]
df.Nov <- df.Nov[order(df.Nov[ ,1]) , ]
df.Dec <- df.Dec[order(df.Dec[ ,1]) , ]

# A "Month" character column is added to each adta.frame
df.Jan$Month <- c("Jan")
df.Feb$Month <- c("Feb")
df.Mar$Month <- c("Mar")
df.Apr$Month <- c("Apr")
df.May$Month <- c("May")
df.Jun$Month <- c("Jun")
df.Jul$Month <- c("Jul")
df.Aug$Month <- c("Aug")
df.Sep$Month <- c("Sep")
df.Oct$Month <- c("Oct")
df.Nov$Month <- c("Nov")
df.Dec$Month <- c("Dec")

# colnames in monthly data.frames are renamed
colnames(df.Jan) <- c("Perc_Exceedance", "FLOW", "MONTH")
colnames(df.Feb) <- c("Perc_Exceedance", "FLOW", "MONTH")
colnames(df.Mar) <- c("Perc_Exceedance", "FLOW", "MONTH")
colnames(df.Apr) <- c("Perc_Exceedance", "FLOW", "MONTH")
colnames(df.May) <- c("Perc_Exceedance", "FLOW", "MONTH")
colnames(df.Jun) <- c("Perc_Exceedance", "FLOW", "MONTH")
colnames(df.Jul) <- c("Perc_Exceedance", "FLOW", "MONTH")
colnames(df.Aug) <- c("Perc_Exceedance", "FLOW", "MONTH")
colnames(df.Sep) <- c("Perc_Exceedance", "FLOW", "MONTH")
colnames(df.Oct) <- c("Perc_Exceedance", "FLOW", "MONTH")
colnames(df.Nov) <- c("Perc_Exceedance", "FLOW", "MONTH")
colnames(df.Dec) <- c("Perc_Exceedance", "FLOW", "MONTH")

# data.frames are binded by row
df.rbind.FDC <- rbind(df.Jan, 
                      df.Feb,
                      df.Mar,
                      df.Apr,
                      df.May,
                      df.Jun,
                      df.Jul,
                      df.Aug,
                      df.Sep,
                      df.Oct,
                      df.Nov,
                      df.Dec)

# Continuous standard monthly Flow Duration Curves are generated
g.cont <- ggplot() +
	         geom_line(aes(x = Perc_Exceedance,y = FLOW,colour = MONTH, linetype = MONTH,group = MONTH),
	         data=df.rbind.FDC,size = 0.75) +
          scale_y_continuous(trans='log10',
                             breaks = c(c(0.1,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,30,40,50,60,70,80,90,100))) +
          scale_x_continuous(breaks = c(c(1,2,3,4,5,10,20,30,40,50,60,70,80,90,100))) + 
          ylab(label = 'Q (m3/s)') +
          xlab(label = 'Percent Exceedance') +
          ggtitle(label = 'Upper Toro River Catchment. Monthly Flow Duration Curves. Continuous Values (1994-2010)') +
          theme_bw(base_size = 18.0)

# Continuous standard monthly Flow Duration Curves are requested
g.cont

# NULL vectors are declared for each month
vt.month.output.Jan <- NULL
vt.month.output.Feb <- NULL
vt.month.output.Mar <- NULL
vt.month.output.Apr <- NULL
vt.month.output.May <- NULL
vt.month.output.Jun <- NULL
vt.month.output.Jul <- NULL
vt.month.output.Aug <- NULL
vt.month.output.Sep <- NULL
vt.month.output.Oct <- NULL
vt.month.output.Nov <- NULL
vt.month.output.Dec <- NULL

# Empty vectors are declared for each month
vt.month.output.Jan <- vector()
vt.month.output.Feb <- vector()
vt.month.output.Mar <- vector()
vt.month.output.Apr <- vector()
vt.month.output.May <- vector()
vt.month.output.Jun <- vector()
vt.month.output.Jul <- vector()
vt.month.output.Aug <- vector()
vt.month.output.Sep <- vector()
vt.month.output.Oct <- vector()
vt.month.output.Nov <- vector()
vt.month.output.Dec <- vector()

# Discrete Percentage Exceedance interval are defined
breaks.fdc=c(0.01,0.1,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99,99.9,99.99)

# A loop is executed to find closest Percentage Exceedance to discrete intervals for each month
for (i in 1:(length(breaks.fdc))) {
  
  vt.month.output.Jan [i] <- which.min(abs((df.Jan$Perc_Exceedance - breaks.fdc[i])))
  vt.month.output.Feb [i] <- which.min(abs((df.Feb$Perc_Exceedance - breaks.fdc[i])))
  vt.month.output.Mar [i] <- which.min(abs((df.Mar$Perc_Exceedance - breaks.fdc[i])))
  vt.month.output.Apr [i] <- which.min(abs((df.Apr$Perc_Exceedance - breaks.fdc[i])))
  vt.month.output.May [i] <- which.min(abs((df.May$Perc_Exceedance - breaks.fdc[i])))
  vt.month.output.Jun [i] <- which.min(abs((df.Jun$Perc_Exceedance - breaks.fdc[i])))
  vt.month.output.Jul [i] <- which.min(abs((df.Jul$Perc_Exceedance - breaks.fdc[i])))
  vt.month.output.Aug [i] <- which.min(abs((df.Aug$Perc_Exceedance - breaks.fdc[i])))
  vt.month.output.Sep [i] <- which.min(abs((df.Sep$Perc_Exceedance - breaks.fdc[i])))
  vt.month.output.Oct [i] <- which.min(abs((df.Oct$Perc_Exceedance - breaks.fdc[i])))
  vt.month.output.Nov [i] <- which.min(abs((df.Nov$Perc_Exceedance - breaks.fdc[i])))
  vt.month.output.Dec [i] <- which.min(abs((df.Dec$Perc_Exceedance - breaks.fdc[i])))
  
}

# Matching values are extracted from monthly data.frames
vt.month.output.Jan <- (df.Jan[vt.month.output.Jan,2])
vt.month.output.Feb <- (df.Feb[vt.month.output.Feb,2])
vt.month.output.Mar <- (df.Mar[vt.month.output.Mar,2])
vt.month.output.Apr <- (df.Apr[vt.month.output.Apr,2])
vt.month.output.May <- (df.May[vt.month.output.May,2])
vt.month.output.Jun <- (df.Jun[vt.month.output.Jun,2])
vt.month.output.Jul <- (df.Jul[vt.month.output.Jul,2])
vt.month.output.Aug <- (df.Aug[vt.month.output.Aug,2])
vt.month.output.Sep <- (df.Sep[vt.month.output.Sep,2])
vt.month.output.Oct <- (df.Oct[vt.month.output.Oct,2])
vt.month.output.Nov <- (df.Nov[vt.month.output.Nov,2])
vt.month.output.Dec <- (df.Dec[vt.month.output.Dec,2])

# A monthly summary data.frame is created based on discrete Percentage Exceedance intervals 
df.comp.FDC <- data.frame(breaks.fdc,
                          vt.month.output.Jan,
                          vt.month.output.Feb,
                          vt.month.output.Mar,
                          vt.month.output.Apr,
                          vt.month.output.May,
                          vt.month.output.Jun,
                          vt.month.output.Jul,
                          vt.month.output.Aug,
                          vt.month.output.Sep,
                          vt.month.output.Oct,
                          vt.month.output.Nov,
                          vt.month.output.Dec)

# Colnames in df.comp.FDC data.frames are renamed
names(df.comp.FDC) <- c("Perc_Exceedance",
                        "Jan",
                        "Feb",
                        "Mar",
                        "Apr",
                        "May",
                        "Jun",
                        "Jul",
                        "Aug",
                        "Sep",
                        "Oct",
                        "Nov",
                        "Dec")

# Monthy character labels are repeated as a function of length(breaks.fdc)
Jan_ch <- rep("Jan",length(breaks.fdc))
Feb_ch <- rep("Feb",length(breaks.fdc))
Mar_ch <- rep("Mar",length(breaks.fdc))
Apr_ch <- rep("Apr",length(breaks.fdc))
May_ch <- rep("May",length(breaks.fdc))
Jun_ch <- rep("Jun",length(breaks.fdc))
Jul_ch <- rep("Jul",length(breaks.fdc))
Aug_ch <- rep("Aug",length(breaks.fdc))
Sep_ch <- rep("Sep",length(breaks.fdc))
Oct_ch <- rep("Oct",length(breaks.fdc))
Nov_ch <- rep("Nov",length(breaks.fdc))
Dic_ch <- rep("Dic",length(breaks.fdc))

# A vector is repeated which length is equal to breaks.fdc is repeated 12 times,
# one for each month
V.C1 <- rep(1:length(breaks.fdc), 12)

# Discrete Percentage Exceedance intervals are repeated 12 times,
# one for each month
V.C2 <- rep(breaks.fdc,12)

# A summation vector is created based on monthly vectors
V.C3 <- c(vt.month.output.Jan,
          vt.month.output.Feb,
          vt.month.output.Mar,
          vt.month.output.Apr,
          vt.month.output.May,
          vt.month.output.Jun,
          vt.month.output.Jul,
          vt.month.output.Aug,
          vt.month.output.Sep,
          vt.month.output.Oct,
          vt.month.output.Nov,
          vt.month.output.Dec)

# A summation vector is created based on monthy character-labeled vectors
V.C4 <- c(Jan_ch,
          Feb_ch,
          Mar_ch,
          Apr_ch,
          May_ch,
          Jun_ch,
          Jul_ch,
          Aug_ch,
          Sep_ch,
          Oct_ch,
          Nov_ch,
          Dic_ch)

# A summary data.frame is created based on all monthly vectors for graphical purposes
df.discrete <- data.frame(V.C1, V.C2, V.C3, V.C4)

# colnames in df.discrete data.frame are renamed
colnames(df.discrete) <- c("Sequence","Perc_Exceedance","FLOW","Month")

# Discrete standard monthly Flow Duration Curves are generated
g.disc <- ggplot() +
	        geom_line(aes(x = Perc_Exceedance,y = FLOW,colour = Month,linetype = Month,group = Month),
	                  data=df.discrete,size = 0.75) +
          scale_y_continuous(trans='log10',
                             breaks = c(c(0.1,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,30,40,50,60,70,80,90,100))) +
          scale_x_continuous(breaks = c(c(1,2,3,4,5,10,20,30,40,50,60,70,80,90,100))) + 
          ylab(label = 'Q (m3/s)') +
          xlab(label = 'Percent Exceedance') +
          ggtitle(label = 'Upper Toro River Catchment. Monthly Flow Duration Curves. Discrete Values (1994-2010)') +
          theme_grey()

# Discrete standard monthly Flow Duration Curves are requested
g.disc

# /////////////////////////////////////////
# BLOCK: Temporal and Graphical Indicators
# /////////////////////////////////////////

# A flow monthly boxplot is generated
boxplot.Month <- ggplot() +
  geom_point(aes(x = MONTH_CH,y = FLOW),data=df.obs,size = 0.5) +
  geom_boxplot(aes(y = FLOW,x = MONTH_CH,colour = MONTH_CH),
               data=df.obs,size = 0.75,alpha = 0.75,outlier.colour = '#ff0000',outlier.size = 1.5) +
  scale_y_continuous(trans='log10',
                     breaks = c(c(0.1,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,30,40,50,60,70,80,90,100))) +
  ylab(label = 'Q (m3/s)') +
  xlab(label = 'Month') +
  ggtitle(label = 'Upper Toro River Catchment. Monthly Boxplot (1994-2010)') +
  theme_bw(base_size = 18.0)

# A flow monthly boxplot is requested
boxplot.Month

# A Log-Scaled flow monthly boxplot-violin.plot is generated
boxplot.Month.violin.log <- ggplot() +
  geom_point(aes(x = MONTH_CH,y = FLOW),data=df.obs,
             size = 0.95,position = position_jitter(width = 0.15)) +
  geom_violin(aes(x = MONTH_CH,y = FLOW,colour = MONTH_CH),data=df.obs,size = 0.85,alpha = 0.75) +
  geom_boxplot(aes(y = FLOW,x = MONTH_CH),data=df.obs,
               size = 0.15,alpha = 0.0,outlier.size = 0.15) +
  scale_y_continuous(trans='log10',
                     breaks = c(c(0.5,2,4,6,8,10,12,14,16,18,20,30,40,50,60,70,80,90,100))) +
  ylab(label = 'Q (m3/s)_log scaled') +
  xlab(label = 'Month') +
  ggtitle(label = 'Upper Toro River Catchment. Monthly Boxplot-Violin (1994-2010)') +
  theme_grey()

# A Log-Scaled flow monthly boxplot-violin.plot is requested
boxplot.Month.violin.log

# A flow monthly boxplot-violin.plot is generated
boxplot.Month.violin <- ggplot() +
  geom_point(aes(x = MONTH_CH,y = FLOW),data=df.obs,
             size = 0.95,position = position_jitter(width = 0.15)) +
  geom_violin(aes(x = MONTH_CH,y = FLOW,colour = MONTH_CH),data=df.obs,size = 0.85,alpha = 0.75) +
  geom_boxplot(aes(y = FLOW,x = MONTH_CH),data=df.obs,size = 0.15,alpha = 0.0,outlier.size = 0.15) +
  scale_y_continuous(breaks = c(c(0.5,2,4,6,8,10,12,14,16,18,20,30,40,50,60,70,80,90,100))) +
  ylab(label = 'Q (m3/s)') +
  xlab(label = 'Month') +
  ggtitle(label = 'Upper Toro River Catchment. Monthly Boxplot-Violin (1994-2010)') +
  theme_grey()

# A flow monthly boxplot-violin.plot is requested
boxplot.Month.violin

# A flow yearly boxplot is generated
boxplot.Year <- ggplot() +
  geom_point(aes(x = YEAR_CH,y = FLOW),data=df.obs,size = 0.5) +
  geom_boxplot(aes(y = FLOW,x = YEAR_CH,colour = YEAR_CH),data=df.obs,
               size = 0.75,alpha = 0.75,outlier.colour = '#ff0000',outlier.size = 1.5) +
  scale_y_continuous(trans='log10',
                     breaks = c(c(0.1,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,30,40,50,60,70,80,90,100))) +
  ylab(label = 'Q (m3/s)_log scaled') +
  xlab(label = 'Year') +
  ggtitle(label = 'Upper Toro River Catchment. Yearly Boxplot (1994-2010)') +
  theme_bw(base_size = 18.0)

# A flow yearly boxplot-violin.plot is requested
boxplot.Year

# ///////////////////////////////////////
# BLOCK: Baseflow Separation and Lowflow
# ///////////////////////////////////////

# BaseflowSeparation function from EcoHydRology library is used to get an approximation
# of baseflow using a 3-pass filter and rounded to 3 significant digit
df.bfs <- round((BaseflowSeparation(df.obs$FLOW, passes=3)),3)

# colnames in df.bfs data.frame are renamed
colnames(df.bfs) <- c("BASEFLOW","QUICKFLOW")

# df.obs data.frame is cbinded to df.bfs data.frame
df.bfs.union <- cbind(df.obs,df.bfs)

# CRAN library data.table is loaded at this point to avoid conflicts with other libraries
require(data.table)

# A subset data.frame is created base on df.bfs.union data.frame
df.bfs.union.sub2 <- data.table(date = as.IDate(df.bfs.union$DATE2), df.bfs.union[-1])

# Baseflow annual statistics are extracted 
# from df.bfs.union.sub2 into a new data.frame called df.bfs.union.annual
df.bfs.union.annual <- as.data.frame(df.bfs.union.sub2[, list(mean.BASEFLOW = mean(BASEFLOW),
                                                   median.BASEFLOW = median(BASEFLOW),
                                                   min.BASEFLOW = min(BASEFLOW),
                                                   max.BASEFLOW = max(BASEFLOW)), 
                                                   by = year(date)])

# df.bfs.union.annual data.frame is rounded to three significant digits
df.bfs.union.annual <- round(df.bfs.union.annual,3)

# CRAN library data.table is detached to avoid conflicts with other libraries
detach(package:data.table)

# A subset data.frame containing only "MONTH_CH" and "BASEFLOW" is created      
df.bfs.union.temp <- df.bfs.union[,c("MONTH_CH","BASEFLOW")]


#//////////////////////////////////////////


# NULL vectors for baseflow are declared for each month
Jan.base <- NULL
Jan.base <- c()
Feb.base <- NULL
Feb.base <- c()
Mar.base <- NULL
Mar.base <- c()
Apr.base <- NULL
Apr.base <- c()
May.base <- NULL
May.base <- c()
Jun.base <- NULL
Jun.base <- c()
Jul.base <- NULL
Jul.base <- c()
Aug.base <- NULL
Aug.base <- c()
Sep.base <- NULL
Sep.base <- c()
Oct.base <- NULL
Oct.base <- c()
Nov.base <- NULL
Nov.base <- c()
Dec.base <- NULL
Dec.base <- c()

# A baseflow counter-loop is declared
counter01 <- length(df.bfs.union.temp$MONTH_CH)

# A disaggregation baseflow loop is executed to separate baseflow by month
for(i in 1:counter01) {

  if(df.bfs.union.temp$MONTH_CH[i] == "Jan" ) {
    Jan.base[i] <- df.bfs.union.temp$BASEFLOW[i] }
  else if(df.bfs.union.temp$MONTH_CH[i] == "Feb" ) {
         Feb.base[i] <- df.bfs.union.temp$BASEFLOW[i] }
  else if(df.bfs.union.temp$MONTH_CH[i] == "Mar" ) {
         Mar.base[i] <- df.bfs.union.temp$BASEFLOW[i] }
  else if(df.bfs.union.temp$MONTH_CH[i] == "Apr" ) {
         Apr.base[i] <- df.bfs.union.temp$BASEFLOW[i] }
  else if(df.bfs.union.temp$MONTH_CH[i] == "May" ) {
         May.base[i] <- df.bfs.union.temp$BASEFLOW[i] }
  else if(df.bfs.union.temp$MONTH_CH[i] == "Jun" ) {
         Jun.base[i] <- df.bfs.union.temp$BASEFLOW[i] }
  else if(df.bfs.union.temp$MONTH_CH[i] == "Jul" ) {
         Jul.base[i] <- df.bfs.union.temp$BASEFLOW[i] }
  else if(df.bfs.union.temp$MONTH_CH[i] == "Aug" ) {
         Aug.base[i] <- df.bfs.union.temp$BASEFLOW[i] }
  else if(df.bfs.union.temp$MONTH_CH[i] == "Sep" ) {
         Sep.base[i] <- df.bfs.union.temp$BASEFLOW[i] }
  else if(df.bfs.union.temp$MONTH_CH[i] == "Oct" ) {
         Oct.base[i] <- df.bfs.union.temp$BASEFLOW[i] }
  else if(df.bfs.union.temp$MONTH_CH[i] == "Nov" ) {
         Nov.base[i] <- df.bfs.union.temp$BASEFLOW[i] }
  else if(df.bfs.union.temp$MONTH_CH[i] == "Dec" ) {
         Dec.base[i] <- df.bfs.union.temp$BASEFLOW[i] }

}

# NA values are removed from monthly baseflow vectors
Jan.base <- Jan.base[!is.na(Jan.base)]
Feb.base <- Feb.base[!is.na(Feb.base)]
Mar.base <- Mar.base[!is.na(Mar.base)]
Apr.base <- Apr.base[!is.na(Apr.base)]
May.base <- May.base[!is.na(May.base)]
Jun.base <- Jun.base[!is.na(Jun.base)]
Jul.base <- Jul.base[!is.na(Jul.base)]
Aug.base <- Aug.base[!is.na(Aug.base)]
Sep.base <- Sep.base[!is.na(Sep.base)]
Oct.base <- Oct.base[!is.na(Oct.base)]
Nov.base <- Nov.base[!is.na(Nov.base)]
Dec.base <- Dec.base[!is.na(Dec.base)]

# Monthly baseflow vectors are converted to data.frames and rounded to 3 significant digits

df.Jan.base <- as.data.frame(round((stat.desc(Jan.base)),3))
df.Feb.base <- as.data.frame(round((stat.desc(Feb.base)),3))
df.Mar.base <- as.data.frame(round((stat.desc(Mar.base)),3))
df.Apr.base <- as.data.frame(round((stat.desc(Apr.base)),3))
df.May.base <- as.data.frame(round((stat.desc(May.base)),3))
df.Jun.base <- as.data.frame(round((stat.desc(Jun.base)),3))
df.Jul.base <- as.data.frame(round((stat.desc(Jul.base)),3))
df.Aug.base <- as.data.frame(round((stat.desc(Aug.base)),3))
df.Sep.base <- as.data.frame(round((stat.desc(Sep.base)),3))
df.Oct.base <- as.data.frame(round((stat.desc(Oct.base)),3))
df.Nov.base <- as.data.frame(round((stat.desc(Nov.base)),3))
df.Dec.base <- as.data.frame(round((stat.desc(Dec.base)),3))

# Monthly baseflow data.frames are binded by column
df.rbind.base <- cbind(df.Jan.base, 
                       df.Feb.base,
                       df.Mar.base,
                       df.Apr.base,
                       df.May.base,
                       df.Jun.base,
                       df.Jul.base,
                       df.Aug.base,
                       df.Sep.base,
                       df.Oct.base,
                       df.Nov.base,
                       df.Dec.base)

# colnames in df.rbind.base data.frame are renamed
colnames(df.rbind.base) <- c("Jan_BASEFLOW",
                             "Feb_BASEFLOW",
                             "Mar_BASEFLOW",
                             "Apr_BASEFLOW",
                             "May_BASEFLOW",
                             "Jun_BASEFLOW",
                             "Jul_BASEFLOW",
                             "Aug_BASEFLOW",
                             "Sep_BASEFLOW",
                             "Oct_BASEFLOW",
                             "Nov_BASEFLOW",
                             "Dec_BASEFLOW")

# A Baseflow monthly boxplot is generated
boxplot.baseflow.Month <- ggplot() +
                          geom_boxplot(aes(y = BASEFLOW,x = MONTH_CH, colour = MONTH_CH),data=df.bfs.union,
                                       size = 0.75,alpha = 0.75,outlier.colour = '#ff0000',outlier.size = 1.5) +
                          scale_y_continuous(breaks = c(c(1,2,3,4,5,6,7,8,9,10,15,20,25,30,50,100))) +
                          ylab(label = 'Q (m3/s)') +
                          xlab(label = 'Month') +
                          ggtitle(label = 'Upper Toro River Catchment. Monthy Baseflow Boxplot (1994-2010)') +
                          theme_grey()

# A Baseflow monthly boxplot is requested
boxplot.baseflow.Month

# A Quickflow-Baseflow ratio is calculated and added to df.bfs.union data.frame
df.bfs.union$RATIO <- round(((df.bfs.union$QUICKFLOW / df.bfs.union$FLOW)),3)

# Descriptive statistics are extracted for df.bfs.union data.frame for "RATIO" column only
df.bfs.desc.RATIO <- round((as.data.frame(stat.desc(df.bfs.union$RATIO))),3)

# A Baseflow monthly boxplot (as fraction of total flow) is generated
boxplot.baseflow.RATIO <- ggplot() +
                          geom_boxplot(aes(y = RATIO,x = MONTH_CH, colour = MONTH_CH),data=df.bfs.union,
                                       size = 0.75,alpha = 0.75,outlier.colour = '#ff0000',outlier.size = 1.5) +
                          scale_y_continuous(breaks = c(c(0.025,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))) +
                          #geom_hline(data=df.bfs.union,size = 0.95,alpha = 0.5,yintercept = (df.bfs.desc.RATIO[8,1])) +
                          #geom_hline(data=df.bfs.union,size = 0.95,linetype = 2,
                                     #alpha = 0.5,yintercept = (df.bfs.desc.RATIO[9,1]), color="black") +
                          geom_text(aes(df.bfs.union[1,8],df.bfs.desc.RATIO[8,1],label = "median", vjust = -0.5)) +
                          (aes(df.bfs.union[1,8],df.bfs.desc.RATIO[9,1],label = "mean", vjust = -1)) +
                          ylab(label = 'Baseflow Index (BFI)') +
                          xlab(label = 'Month') +
                          ggtitle(label = 'Upper Toro River Catchment. Monthy Baseflow Boxplot (as fraction of total flow)  (1994-2010)') +
                          theme_bw(base_size = 18.0)

# A Baseflow monthly boxplot (as fraction of total flow) is requested
boxplot.baseflow.RATIO

# A complete-series hydrograph is generated which includes Baseflow
# Data points over FFATr1.5 and FFATr10 are shown above threshold lines
g.hydro04 <- ggplot() +
  #geom_point(aes(x = DATE2,y = FLOW,colour = YEAR_CH),data=df.obs,size = 1.50) +
  geom_point(aes(x = DATE2,y = FLOW),data=df.FFATr1.5_sel,shape = 21,colour = '#999900',size = 3.0) +
  geom_line(aes(x = DATE2,y = FLOW),data=df.obs,colour = '#666666',size = 0.15,alpha = 0.75) +
  geom_line(aes(x = DATE2,y = BASEFLOW),data=df.bfs.union,colour = '#0000cc',size = 0.75) +
  scale_y_continuous(trans='log10',breaks = c(c(0.1,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12,14,16,18,20,30,40,50,60,70,80,90,100))) +
  scale_x_date(breaks = scales::pretty_breaks(n = 15.0)) +
  geom_hline(data=df.obs,size = 0.95,alpha = 0.5,yintercept = (df.obs.desc[8,1])) +
  geom_hline(data=df.obs,size = 0.95,linetype = 2,alpha = 0.5,yintercept = (df.obs.desc[9,1]), color="black") +
  geom_hline(data=df.obs,,size = 0.95,colour = '#999900',yintercept = limit1.5) +
  geom_hline(data=df.obs,,size = 0.95,colour = '#999900',linetype = 2,yintercept = limit10) +
  ylab(label = 'Q (m3/s)') +
  xlab(label = 'Period (years)') +
  ggtitle(label = 'Upper Toro River Catchment Streamflow Analisis (1994-2010)') +
  geom_text(aes(df.obs[1,4],limit1.5,label = "     R.int > 1.5 years", vjust = -1)) +
  geom_text(aes(df.obs[1,4],limit10,label = "      R.int > 10 years", vjust = -1)) +
  geom_text(aes(df.obs[1,4],df.obs.desc[8,1],label = "median", vjust = -0.5)) +
  geom_text(aes(df.obs[1,4],df.obs.desc[9,1],label = "mean", vjust = -1)) +
  geom_text(aes(df.obs[1,4],df.obs.desc[9,1],label = "Baseflow", vjust = 6)) +
  theme_bw(base_size = 18.0)

# Baseflow annual-summary hydrograph is requested
g.hydro04

# //////////////////////////////
# BLOCK: Flow Assessment Models
# //////////////////////////////

# Maximum flow values are sorted in df.annual.flow data.frame
df.annual.flow.model <- df.annual.flow[order(df.annual.flow$year),]

# Character columns are created in df.annual.flow.model data.frame for graphical purposes
 df.annual.flow.model$data <- c("observations")
df.annual.flow.model$model <- c("CI_95")

# Character columns are created in df.bfs.union.annual data.frame for graphical purposes
df.bfs.union.annual$data <- c("observations")
df.bfs.union.annual$model <- c("CI_95")

# A LM model for the Annual-MedianFlow is created
lm.01 <- lm(median.FLOW ~ year, data=df.annual.flow.model)
sm.lm.01 <- summary(lm.01)

# A LM model for the Annual-MeanFlow is created
lm.02 <- lm(mean.FLOW ~ year, data=df.annual.flow.model)
sm.lm.02 <- summary(lm.02)

# A LM model for the Annual-MaxFlow is created
lm.03 <- lm(max.FLOW ~ year, data=df.annual.flow.model)
sm.lm.03 <- summary(lm.03)

# A LM model for the Annual-MedianBaseFlow is created
lm.04 <- lm(median.BASEFLOW ~ year, data=df.bfs.union.annual)
sm.lm.04 <- summary(lm.04)

# A LM model for the Annual-MeanBaseFlow is created
lm.05 <- lm(mean.BASEFLOW ~ year, data=df.bfs.union.annual)
sm.lm.05 <- summary(lm.05)

# A LM model for the Annual-MaxBaseFlow is created
lm.06 <- lm(max.BASEFLOW ~ year, data=df.bfs.union.annual)
sm.lm.06 <- summary(lm.06)

# A LM model-plot for the Annual-Median is created 
g.lm01 <- ggplot(data=df.annual.flow.model) +
          geom_smooth(aes(x = year,y = median.FLOW,colour = model),
                      fill = '#cccc00',size = 1.25,alpha = 0.25,method = lm) +
          geom_line(aes(x = year,y = median.FLOW),size = 0.45,linetype = 2,alpha = 0.55) +
          geom_point(aes(x = year,y = median.FLOW,fill = data),colour = '#333333',size = 4.0,alpha = 0.99) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0),
	                            expand = c(0.05,((max(df.annual.flow.model$median.FLOW))*0.075))) +
	         scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
          theme_bw(base_size = 18.0) +
          ylab(label = 'Q (m3/s)') +
          xlab(label = 'Period (years)') +
          ggtitle(label = 'Upper Toro River Catchment. Annual-Median Linear-Model (1994-2010)')

# boundary.LM function is called
post.lm01 <- boundary.LM(max(df.annual.flow.model$year),
                         min(df.annual.flow.model$year),
                         max(df.annual.flow.model$median.FLOW),
                         min(df.annual.flow.model$median.FLOW))

# Boundaries for LM model-plot for the Annual-Max labels are incorporated
g.lm01 <- g.lm01 + (geom_text(data=df.annual.flow.model,x = post.lm01[1],y = post.lm01[2],
                    label = eq.PAR(lm.01), # eq.PAR function is called
                    colour = 'black',hjust = 0.0,vjust = 0.50,alpha = 0.05,size=5.0, parse = FALSE))

# A LM model-plot for the Annual-Median is requested
g.lm01 

# A LM model-plot for the Annual-mean is created 
g.lm02 <- ggplot(data=df.annual.flow.model) +
          geom_smooth(aes(x = year,y = mean.FLOW,colour = model),
                      fill = '#cccc00',size = 1.25,alpha = 0.25,method = lm) +
          geom_line(aes(x = year,y = mean.FLOW),size = 0.45,linetype = 2,alpha = 0.55) +
          geom_point(aes(x = year,y = mean.FLOW,fill = data),colour = '#333333',size = 4.0,alpha = 0.99) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0),
                             expand = c(0.05,((max(df.annual.flow.model$mean.FLOW))*0.075))) +
          scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
          theme_bw(base_size = 18.0) +
          ylab(label = 'Q (m3/s)') +
          xlab(label = 'Period (years)') +
          ggtitle(label = 'Upper Toro River Catchment. Annual-Mean Linear-Model (1994-2010)')

# boundary.LM function is called
post.lm02 <- boundary.LM(max(df.annual.flow.model$year),
                         min(df.annual.flow.model$year),
                         max(df.annual.flow.model$mean.FLOW),
                         min(df.annual.flow.model$mean.FLOW))

# Boundaries for LM model-plot for the Annual-Max labels are incorporated
g.lm02 <- g.lm02 + (geom_text(data=df.annual.flow.model,x = post.lm02[1],y = post.lm02[2],
                    label = eq.PAR(lm.02), # eq.PAR function is called
                    colour = 'black',hjust = 0.0,vjust = 0.50,alpha = 0.05,size=5.0, parse = FALSE))

# A LM model-plot for the Annual-Mean is requested
g.lm02 

# A LM model-plot for the Annual-Max is created 
g.lm03 <- ggplot(data=df.annual.flow.model) +
          geom_smooth(aes(x = year,y = max.FLOW,colour = model),
                      fill = '#cccc00',size = 1.25,alpha = 0.25,method = lm) +
          geom_line(aes(x = year,y = max.FLOW),size = 0.45,linetype = 2,alpha = 0.55) +
          geom_point(aes(x = year,y = max.FLOW,fill = data),colour = '#333333',size = 4.0,alpha = 0.99) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0),
                             expand = c(0.05,((max(df.annual.flow.model$max.FLOW))*0.075))) +
          scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
          theme_bw(base_size = 18.0) +
          ylab(label = 'Q (m3/s)') +
          xlab(label = 'Period (years)') +
          ggtitle(label = 'Upper Toro River Catchment. Annual-Max Linear-Model (1994-2010)')
      
# boundary.LM function is called
post.lm03 <- boundary.LM(max(df.annual.flow.model$year),
                         min(df.annual.flow.model$year),
                         max(df.annual.flow.model$max.FLOW),
                         min(df.annual.flow.model$max.FLOW))

# Boundaries for LM model-plot for the Annual-Max labels are incorporated
g.lm03 <- g.lm03 + (geom_text(data=df.annual.flow.model,x = post.lm03[1],y = post.lm03[2],
                    label = eq.PAR(lm.03), # eq.PAR function is called
                    colour = 'black',hjust = 0.0,vjust = 0.50,alpha = 0.05,size=5.0, parse = FALSE))

# A LM model-plot for the Annual-Max is requested
g.lm03 

# A LM model-plot for the Annual-median.BASEFLOW is created 
g.lm04 <- ggplot(data=df.bfs.union.annual) +
          geom_smooth(aes(x = year,y = median.BASEFLOW,colour = model),
                      fill = 'green',size = 1.25,alpha = 0.25,method = lm) +
          geom_line(aes(x = year,y = median.BASEFLOW),size = 0.45,linetype = 2,alpha = 0.55) +
          geom_point(aes(x = year,y = median.BASEFLOW,fill = data),colour = '#333333',size = 4.0,alpha = 0.99) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0),
                            expand = c(0.05,((max(df.bfs.union.annual$median.BASEFLOW))*0.075))) +
          scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
          theme_bw(base_size = 18.0) +
          ylab(label = 'Q (m3/s)') +
          xlab(label = 'Period (years)') +
          ggtitle(label = 'Upper Toro River Catchment. Annual-Median BaseFlow Linear-Model (1994-2010)')

# boundary.LM function is called
post.lm04 <- boundary.LM(max(df.bfs.union.annual$year),
                         min(df.bfs.union.annual$year),
                         max(df.bfs.union.annual$median.BASEFLOW),
                         min(df.bfs.union.annual$median.BASEFLOW))

# Boundaries for LM model-plot for the Annual-median.BASEFLOW labels are incorporated
g.lm04 <- g.lm04 + (geom_text(data=df.bfs.union.annual,x = post.lm04[1],y = post.lm04[2],
                              label = eq.PAR(lm.04), # eq.PAR function is called
                              colour = 'black',hjust = 0.0,vjust = 0.50,alpha = 0.05,size=5.0, parse = FALSE))

# A LM model-plot for the Annual-median.BASEFLOW is requested
g.lm04 

# A LM model-plot for the Annual-mean.BASEFLOW is created 
g.lm05 <- ggplot(data=df.bfs.union.annual) +
          geom_smooth(aes(x = year,y = mean.BASEFLOW,colour = model),
                      fill = 'green',size = 1.25,alpha = 0.25,method = lm) +
          geom_line(aes(x = year,y = mean.BASEFLOW),size = 0.45,linetype = 2,alpha = 0.55) +
          geom_point(aes(x = year,y = mean.BASEFLOW,fill = data),colour = '#333333',size = 4.0,alpha = 0.99) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0),
                             expand = c(0.05,((max(df.bfs.union.annual$mean.BASEFLOW))*0.075))) +
          scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
          theme_bw(base_size = 18.0) +
          ylab(label = 'Q (m3/s)') +
          xlab(label = 'Period (years)') +
          ggtitle(label = 'Upper Toro River Catchment. Annual-Mean BaseFlow Linear-Model (1994-2010)')

# boundary.LM function is called
post.lm05 <- boundary.LM(max(df.bfs.union.annual$year),
                         min(df.bfs.union.annual$year),
                         max(df.bfs.union.annual$mean.BASEFLOW),
                         min(df.bfs.union.annual$mean.BASEFLOW))

# Boundaries for LM model-plot for the Annual-mean.BASEFLOW labels are incorporated
g.lm05 <- g.lm05 + (geom_text(data=df.bfs.union.annual,x = post.lm05[1],y = post.lm05[2],
                              label = eq.PAR(lm.05), # eq.PAR function is called
                              colour = 'black',hjust = 0.0,vjust = 0.50,alpha = 0.05,size=5.0, parse = FALSE))

# A LM model-plot for the Annual-mean.BASEFLOW is requested
g.lm05 

# A LM model-plot for the Annual-max.BASEFLOW is created 
g.lm06 <- ggplot(data=df.bfs.union.annual) +
          geom_smooth(aes(x = year,y = max.BASEFLOW,colour = model),
                      fill = 'green',size = 1.25,alpha = 0.25,method = lm) +
          geom_line(aes(x = year,y = max.BASEFLOW),size = 0.45,linetype = 2,alpha = 0.55) +
          geom_point(aes(x = year,y = max.BASEFLOW,fill = data),colour = '#333333',size = 4.0,alpha = 0.99) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0),
                             expand = c(0.05,((max(df.bfs.union.annual$max.BASEFLOW))*0.075))) +
          scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
          theme_bw(base_size = 18.0) +
          ylab(label = 'Q (m3/s)') +
          xlab(expression(atop("Period (years)", paste("data.frame-source: df.bfs.union.annual")))) +
          ggtitle(label = 'Upper Toro River Catchment. Annual-Max BaseFlow Linear-Model (1994-2010)')
    
# boundary.LM function is called
post.lm06 <- boundary.LM(max(df.bfs.union.annual$year),
                         min(df.bfs.union.annual$year),
                         max(df.bfs.union.annual$max.BASEFLOW),
                         min(df.bfs.union.annual$max.BASEFLOW))

# Boundaries for LM model-plot for the Annual-max.BASEFLOW labels are incorporated
g.lm06 <- g.lm06 + (geom_text(data=df.bfs.union.annual,x = post.lm06[1],y = post.lm06[2],
                              label = eq.PAR(lm.06), # eq.PAR function is called
                              colour = 'black',hjust = 0.0,vjust = 0.50,alpha = 0.05,size=5.0, parse = FALSE))
 
# A LM model-plot for the Annual-max.BASEFLOW is requested
g.lm06 

# /////////////////////////////////////////////
# BLOCK: Export and display and of data.frames
# /////////////////////////////////////////////

# Observed daily flow data.frame
write.csv(df.obs, file = "df.obs.csv")  

# Observed daily flow descriptive statistics data.frame
write.csv(df.obs.desc, file = "df.obs.desc.csv") 

# Observed daily flow  disaggregated by year data.frame					
write.csv(df.pivot, file = "df.pivot.csv") 

# Observed daily flow  descriptive statistics disaggregated by year data.frame							
write.csv(df.pivot.desc, file = "df.pivot.desc.csv") 

# Annual descriptive statistics data.frame				
write.csv(df.annual.flow, file = "df.annual.flow.csv") 

# Estimated-deterministic fitted values, non-exceedance probabilities (nep), 
# return period (rp) and estimated FLOW data.frame
write.csv(df.fa.out, file = "df.fa.out.csv") 

# Non-exceedance probabilities (nep) for "central", "lower" and "upper" CI data.frame
write.csv(df.ci.out, file = "df.ci.out.csv")

# Observed daily flow descriptive statistics disaggregated by month data.frame								
write.csv(df.obs.month.desc, file = "df.obs.month.desc.csv")

# Summarized monthly flows by discrete Percentage Exceedance intervals data.frame								
write.csv(df.comp.FDC, file = "df.comp.FDC.csv")

# Summarized monthly flows by continuous Percentage Exceedance intervals data.frame							
write.csv(df.Jan, file = "df.Jan.csv")
write.csv(df.Feb, file = "df.Feb.csv")
write.csv(df.Mar, file = "df.Mar.csv")
write.csv(df.Apr, file = "df.Apr.csv")
write.csv(df.May, file = "df.May.csv")
write.csv(df.Jun, file = "df.Jun.csv")
write.csv(df.Jul, file = "df.Jul.csv")
write.csv(df.Aug, file = "df.Aug.csv")
write.csv(df.Sep, file = "df.Sep.csv")
write.csv(df.Oct, file = "df.Oct.csv")
write.csv(df.Nov, file = "df.Nov.csv")
write.csv(df.Dec, file = "df.Dec.csv")

# Daily baseflow-quickflow separation	data.frame		
write.csv(df.bfs, file = "df.bfs.csv")

# Observed flow daily including baseflow and quickflow data.frame					
write.csv(df.bfs.union, file = "df.bfs.union.csv")

# Annual baseflow descriptive statistics data.frame		
write.csv(df.bfs.union.annual, file = "df.bfs.union.annual.csv")

# names(df.rbind.base) # baseflow daily descriptive statistics disaggregated by month							
write.csv(df.rbind.base, file = "df.rbind.base.csv")

# LM-Models summaries are requested, saved and exported as *.TXT 
capture.output(sm.lm.01, file = "lm01_median.txt")
capture.output(sm.lm.02, file = "lm02_mean.txt")
capture.output(sm.lm.03, file = "lm03_max.txt")
capture.output(sm.lm.04, file = "lm04_median_baseflow.txt")
capture.output(sm.lm.05, file = "lm05_mean_baseflow.txt")
capture.output(sm.lm.06, file = "lm06_max_baseflow.txt")

# Relevant data.frames are displayed
#View(df.obs)
#View(df.obs.desc)
#View(df.pivot)
#View(df.pivot.desc)
#View(df.annual.flow)
#View(df.fa.out)
#View(df.ci.out)
#View(df.obs.month.desc)
#View(df.comp.FDC)
#View(df.bfs)
#View(df.bfs.union)
#View(df.bfs.union.annual)
#View(df.rbind.base)
