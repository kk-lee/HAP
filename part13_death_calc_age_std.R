# load packages####
library(tableone)
library(caret)
library(reshape2)
library(tidyr)
library(dplyr)
library(binom)
library(survival)
library(gridExtra)
library(grid)
library(survMisc)
library(scales)
library(ggplot)
library(tidyr)
library(metafor)
library(xlsx)
library(xlsxjars)
library(rJava)
library(readxl)
library(stringr)
library(reshape)
library(data.table)
library(plyr)
library(RColorBrewer)
library(rworldmap)
library(classInt)
library(gdata)
library(forestplot)
library(Hmisc)
library(data.table)
library(tidyverse)

# load beta pharm function ####
beta.parms.from.quantiles <- function(q, p=c(0.025,0.975),
                                      precision=0.001, derivative.epsilon=1e-3, start.with.normal.approx=T, start=c(1, 1), plot=F)
{
  # Version 1.3 (February 2017)
  #
  # Function developed by 
  # Lawrence Joseph and pbelisle
  # Division of Clinical Epidemiology
  # Montreal General Hospital
  # Montreal, Qc, Can
  #
  # patrick.belisle@rimuhc.ca
  # http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/BetaParmsFromQuantiles.html
  #
  # Please refer to our webpage for details on each argument.
  
  f <- function(x, theta){dbeta(x, shape1=theta[1], shape2=theta[2])}
  F.inv <- function(x, theta){qbeta(x, shape1=theta[1], shape2=theta[2])}
  f.cum <- function(x, theta){pbeta(x, shape1=theta[1], shape2=theta[2])}
  f.mode <- function(theta){a <- theta[1]; b <- theta[2]; mode <- ifelse(a>1, (a-1)/(a+b-2), NA); mode}
  theta.from.moments <- function(m, v){a <- m*m*(1-m)/v-m; b <- a*(1/m-1); c(a, b)}
  plot.xlim <- c(0, 1)
  
  dens.label <- 'dbeta'
  parms.names <- c('a', 'b')
  
  if (length(p) != 2) stop("Vector of probabilities p must be of length 2.")
  if (length(q) != 2) stop("Vector of quantiles q must be of length 2.")
  p <- sort(p); q <- sort(q)
  
  #_____________________________________________________________________________________________________
  
  print.area.text <- function(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim, M=30, M0=50)
  {
    par.usr <- par('usr')
    par.din <- par('din')
    
    p.string <- as.character(round(c(0,1) + c(1,-1)*p.check, digits=4))
    str.width <- strwidth(p.string, cex=cex)
    str.height <- strheight("0", cex=cex)
    
    J <- matrix(1, nrow=M0, ncol=1)
    
    x.units.1in <- diff(par.usr[c(1,2)])/par.din[1]
    y.units.1in <- diff(par.usr[c(3,4)])/par.din[2]
    aspect.ratio <- y.units.1in/x.units.1in
    
    # --- left area  -----------------------------------------------------------
    
    scatter.xlim <- c(max(plot.xlim[1], par.usr[1]), q[1])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=0, to=p[1], length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-1]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the right of the mode, if any
    w <- which(x>mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[1]+str.width[1]) <= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y, fromLast=T))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[1], adj=c(1,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[1], mean(par.usr[c(3,4)]), labels=p.string[1], col='gray', cex=cex, srt=90, adj=c(1,0))
    }
    
    # --- right area  ----------------------------------------------------------
    
    scatter.xlim <- c(q[2], plot.xlim[2])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=p[2], to=f.cum(plot.xlim[2], theta), length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-length(h)]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the left of the mode, if any
    w <- which(x<mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[2]-str.width[2]) >= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[2], adj=c(0,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[2], mean(par.usr[c(3,4)]), labels=p.string[2], col='gray', cex=cex, srt=-90, adj=c(1,0))
    }
  }
  
  # ......................................................................................................................................
  
  Newton.Raphson <- function(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start)
  {
    Hessian <- matrix(NA, 2, 2)
    
    if (start.with.normal.approx)
    {
      # Probably not a very good universal choice, but proved good in most cases in practice
      m <-  diff(q)/diff(p)*(0.5-p[1]) + q[1]
      v <- (diff(q)/diff(qnorm(p)))^2
      theta <- theta.from.moments(m, v)
    }
    else theta <- start
    
    
    change <- precision + 1
    niter <- 0
    # Newton-Raphson multivariate algorithm
    while (max(abs(change)) > precision)
    {
      Hessian[,1] <- (f.cum(q, theta) - f.cum(q, theta - c(derivative.epsilon, 0))) / derivative.epsilon
      Hessian[,2] <- (f.cum(q, theta) - f.cum(q, theta - c(0, derivative.epsilon))) / derivative.epsilon
      
      f <- f.cum(q, theta) - p
      change <- solve(Hessian) %*% f
      last.theta <- theta
      theta <- last.theta - change
      
      # If we step out of limits, reduce change
      
      if (any(theta<0))
      {
        w <- which(theta<0)
        k <- min(last.theta[w]/change[w])
        theta <- last.theta - k/2*change
      }
      
      niter <- niter + 1
    }
    
    list(theta=as.vector(theta), niter=niter, last.change=as.vector(change))
  }
  
  # ...............................................................................................................
  
  plot.density <- function(p, q, f, f.cum, F.inv, mode, theta, plot.xlim, dens.label, parms.names, cex)
  {
    if (length(plot.xlim) == 0)
    {
      plot.xlim <- F.inv(c(0, 1), theta)
      
      if (is.infinite(plot.xlim[1]))
      {
        tmp <- min(c(0.001, p[1]/10))
        plot.xlim[1] <- F.inv(tmp, theta)
      }  
      
      if (is.infinite(plot.xlim[2]))
      {
        tmp <- max(c(0.999, 1 - (1-p[2])/10))
        plot.xlim[2] <- F.inv(tmp, theta)
      }
    }
    plot.xlim <- sort(plot.xlim)
    
    
    x <- seq(from=min(plot.xlim), to=max(plot.xlim), length=1000)
    h <- f(x, theta)
    x0 <- x; f0 <- h
    ylab <- paste(c(dens.label, '(x, ', parms.names[1], ' = ', round(theta[1], digits=5), ', ', parms.names[2], ' = ', round(theta[2], digits=5), ')'), collapse='')
    plot(x, h, type='l', ylab=ylab)
    
    # fill in area on the left side of the distribution
    x <- seq(from=plot.xlim[1], to=q[1], length=1000)
    y <- f(x, theta)
    x <- c(x, q[1], plot.xlim[1]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # fill in area on the right side of the distribution
    x <- seq(from=max(plot.xlim), to=q[2], length=1000)
    y <- f(x, theta)
    x <- c(x, q[2], plot.xlim[2]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # draw distrn again
    points(x0, f0, type='l')
    h <- f(q, theta)
    points(rep(q[1], 2), c(0, h[1]), type='l', col='orange')
    points(rep(q[2], 2), c(0, h[2]), type='l', col='orange')
    # place text on both ends areas
    print.area.text(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim)  
    
    xaxp <- par("xaxp")
    x.ticks <- seq(from=xaxp[1], to=xaxp[2], length=xaxp[3]+1)
    q2print <- as.double(setdiff(as.character(q), as.character(x.ticks)))
    
    mtext(q2print, side=1, col='orange', at=q2print, cex=0.6, line=2.1)
    points(q, rep(par('usr')[3]+0.15*par('cxy')[2], 2), pch=17, col='orange')
  }
  
  #________________________________________________________________________________________________________________
  
  
  parms <- Newton.Raphson(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start=start)
  p.check <- f.cum(q, parms$theta)
  
  if (plot) plot.density(p, q, f, f.cum, F.inv, f.mode(parms$theta), parms$theta, plot.xlim, dens.label, parms.names, 0.8)
  
  list(a=parms$theta[1], b=parms$theta[2], last.change=parms$last.change, niter=parms$niter, q=q, p=p, p.check=p.check)
}

###=========================================================
### Part 2: estimating death using death and prevalence data
###=========================================================
##=================================
## Part 2a - load prevalence files
##=================================

fuel <- readRDS("data_derived/fuel.rds")

##=================================
## Part 2b - load death files
##=================================

death_std <- readRDS("data_raw/ihme_age_std.rds")

death_std <- death_std %>% 
  filter(measure_name=="Deaths" & 
           sex_name=="Both" &
           year %in% c(2000:2017))

death <- death_std %>% 
  rename.vars(from=c("location_name", "cause_name","val", "upper", "lower"),
              to=c("country", "subgroup","death.ce", "death.ul", "death.ll")) %>% 
  select(country, year, subgroup,death.ce, death.ll,death.ul) %>% 
  mutate_at(vars(country, subgroup), funs(as.factor(.)))

death <- death %>% 
  mutate(country= recode(country, 
                         "The Bahamas" = "Bahamas",
                         "Bolivia" = "Bolivia (Plurinational State of)",
                         "Brunei" = "Brunei Darussalam",
                         "Cape Verde" = "Cabo Verde",
                         "Czech Republic" = "Czechia",
                         "North Korea" = "Democratic People's Republic of Korea",
                         "The Gambia" = "Gambia",
                         "Iran" = "Iran (Islamic Republic of)",
                         "Laos" = "Lao People's Democratic Republic",
                         "Federated States of Micronesia" = "Micronesia (Federated States of)",
                         "South Korea" = "Republic of Korea",
                         "Syria" = "Syrian Arab Republic",
                         "Macedonia" = "The former Yugoslav Republic of Macedonia",
                         "Tanzania" = "United Republic of Tanzania",
                         "United States" = "United States of America",
                         "Venezuela" = "Venezuela (Bolivarian Republic of)",
                         "Vietnam" = "Viet Nam",
                         "Moldova" = "Republic of Moldova")) %>% 
  mutate(subgroup= recode(subgroup,
                          "Ischemic heart disease" = "ihd",
                          "Lower respiratory infections" = "lri",
                          "Upper respiratory infections" = "uri",
                          "Stroke" = "stroke",
                          "Tracheal, bronchus, and lung cancer" = "lungca",
                          "Tuberculosis" = "tb",
                          "Asthma" = "asthma",
                          "Chronic obstructive pulmonary disease" = "copd"))%>%
  mutate(subgroup = case_when(
    subgroup %in% c("asthma", "copd", "lri", "uri", "lungca", "tb") ~ "resp death",
    subgroup %in% c("stroke", "ihd") ~ "cv death",
    TRUE ~ as.character(subgroup)))


# combine all resp and cv deaths

death <- death %>% 
  group_by(country, year, subgroup) %>% 
  summarise_at(.vars= c("death.ce", "death.ll", "death.ul"),
               funs(sum))  %>% 
  filter(year %in% c(2000:2017))%>% 
  ungroup()

##=========================================
## Part 2c - merge prevelance data to death
##=========================================

final <- merge(x=fuel,y=death,by.x = c("country","year"), by.y=c("country","year")) 


## filter out country/years with zero prevalence 
prev_zero <- final %>% 
  filter(fuel.ce==0)

## calculate death for prevalence above zero
final <- final %>% 
  filter(fuel.ce>0)

# checks
# dat.fuel <- data.frame(fuel$country)
# dat.death <- data.frame(death$country)
# mismatch <- anti_join(dat.fuel, dat.death, by= c("fuel.country" = "death.country"))
# mismatch2 <- anti_join(dat.death, dat.fuel, by= c("death.country"="fuel.country"))

## derive alpha and beta using beta pharm function for prevalence data 

alpha <- function(x){
  alpha <- as.vector(beta.parms.from.quantiles(c(as.numeric(x[4]),as.numeric(x[5])))[1])
}

beta <- function(x){
  beta <- as.vector(beta.parms.from.quantiles(c(as.numeric(x[4]),as.numeric(x[5])))[2])
}

final$alpha <- as.vector(unlist(apply(final, 1,alpha)))
final$beta <- as.vector(unlist(apply(final, 1,beta)))


# check beta results
final$uniq_id <- seq(1000, 1000+nrow(final)-1, 1)

check_res <- by(final, final$uniq_id, function(x) qbeta(shape1 = x$alpha, shape2= x$beta, 
                                                        p = c(0.025, 0.5, 0.975)))
check_res <- do.call(rbind, check_res)
colnames(check_res) <- c("est", "lci", "uci")
check_res <- cbind(final, check_res)
plot(check_res$fuel.ll, check_res$lci)
rm(check_res)
final$uniq_id <- NULL

a <- as.list(NULL)
num_iter <- 10000

## derive central estimate and 95% uncertainty term for death using log distribution

final[,c("death.ce","death.ll","death.ul")] <- lapply(final[,c("death.ce","death.ll","death.ul")],function(x){log(x,base = exp(1))})

final$nm_se <- (as.numeric(final$death.ul)-as.numeric(final$death.ll))/3.92
final$nm_mean <- as.numeric(final$death.ce)

# check if log distribution appropriate

final2 <- final
final2[, c('fuel.ce', 'fuel.ll', 'fuel.ul')] <- lapply(final2[, c('fuel.ce', 'fuel.ll', 'fuel.ul')],log) 
final2$upper_diff <- final2$fuel.ul - final2$fuel.ce
final2$lower_diff <- final2$fuel.ce - final2$fuel.ll
final2$diff_diff <- final2$upper_diff - final2$lower_diff
hist(final2$diff_diff, breaks = 100)

##=========================================
## Part 2d - merge rr parameters into final dataset
##=========================================

rr.sim.para.death <- read.csv("data_derived.v2/rr.sim.para.death.csv")
rr.sim.para.death <- rr.sim.para.death %>% 
  select(endpoints, rr.beta, rr.se) %>% 
  rename.vars(from=c("endpoints"),
              to=c("subgroup")) 

final <- inner_join(final, rr.sim.para.death)

###=======================================================
### Part 2e - calculate paf and corresponding death due to IAP 
###=======================================================
num_iter <- 10000

newvars <- data.frame(id = 1:nrow(final), 
                      deathce = NA,
                      deathll = NA,
                      deathul = NA,
                      pafce= NA,
                      pafll = NA,
                      paful = NA)
newvars$id <- NULL


final <- cbind(final, newvars)


for(i in 1:nrow(final)){
  x <- rbeta(num_iter, shape1 = final$alpha[i], shape2 = final$beta[i])
  y <- rnorm(num_iter,final$nm_mean[i],final$nm_se[i])
  y <- exp(y)
  z <- exp(rnorm(num_iter, as.vector(final$rr.beta[i]), final$rr.se[i]))
  p <- (x*(z-1))/(1+x*(z-1))
  death <- p*y
  
  final$deathce[i] <- quantile(death,0.5)
  final$deathll[i] <- quantile(death,0.025)
  final$deathul[i] <- quantile(death,0.975)
  final$pafce[i] <- quantile(p,0.5)
  final$pafll[i] <- quantile(p,0.025)
  final$paful[i] <- quantile(p,0.975)
}


# in thousands

final$deathce <-final$deathce/1000
final$deathll <-final$deathll/1000
final$deathul <-final$deathul/1000

#### add iso3 code

iso <- readRDS("data_derived/iso.rds")

final <- left_join(final,iso) %>% 
  select(country, iso3, region, year, subgroup, deathce, deathll, deathul)

### add prevalence zero country/years

final_prev_zero <- left_join(prev_zero, iso) %>% 
  select(country, iso3, region, year, subgroup) %>% 
  mutate(deathce=0, deathll=0, deathul=0)

final <- rbind(final, final_prev_zero)

final <- final %>% 
  filter(year==2017) %>% 
  group_by(country, iso3, region) %>% 
  summarise_at(.vars= c("deathce", "deathll", "deathul"),
               funs(sum)) 
  
final <- final %>% 
  mutate_at(vars(deathce, deathll, deathul), funs(.*1000)) %>% 
  mutate(death=paste0(round(deathce, digits=0)," [", round(deathll, digits=0), " - ", round(deathul, digits=0), "]")) %>% 
  select(country, death) 
  
#### write csv
write.csv(final, "results.v2/country_death_agestd.table.csv")


### country table
burden <- read.csv("results.v2/country_burden.table.csv") %>% select(-X)
burden_age <- read.csv("results.v2/country_burden_agestd.table.csv") %>% select(-X) %>% 
  dplyr::rename(burden.age=burden)
death <- read.csv("results.v2/country_death.table.csv")%>% select(-X)
death_age <- read.csv("results.v2/country_death_agestd.table.csv")%>% select(-c(X,iso3)) %>% 
  dplyr::rename(death.age=death)

country_table <- left_join(burden, burden_age) %>% 
  left_join(.,death) %>% 
  left_join(.,death_age)

write.csv(country_table, "results.v2/country_burden.death_table.csv")







