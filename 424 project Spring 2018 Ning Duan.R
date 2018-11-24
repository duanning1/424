#
# 1. S&P 500 index (vfinx)
# 2. Alaska Air (alk)
# 3. Costco (cost)
# 4. F5 (ffiv)
# 5. Nordstrom (jwn)

options(digits=3, width=70)
# load packages
library(IntroCompFinR)
library("PerformanceAnalytics")
library("tseries")
library("zoo")
library("boot")
library("corrplot")
savePath="Downloads/"

resetPar = par()

# You will also need to download portfolio_noshorts.r to help compute portfolio calculations
# source portfolio functions
loadPath = "~//Downloads/"
source(file=paste(loadPath, "portfolio_noshorts.r", sep=""))


#
# load data from Yahoo!
#

asset.names = c("alk","cost","ffiv","jwn")

start.date = "2013-04-01"
end.date = "2018-04-30"

vfinx.prices = get.hist.quote(instrument="vfinx", start=start.date,
                              end=end.date, quote="AdjClose",
                              provider="yahoo", origin="1970-01-01",
                              compression="m", retclass="zoo") 
alk.prices = get.hist.quote(instrument="alk", start=start.date,
                             end=end.date, quote="AdjClose",
                             provider="yahoo", origin="1970-01-01",
                             compression="m", retclass="zoo")
cost.prices = get.hist.quote(instrument="cost", start=start.date,
                             end=end.date, quote="AdjClose",
                             provider="yahoo", origin="1970-01-01",
                             compression="m", retclass="zoo")
ffiv.prices = get.hist.quote(instrument="ffiv", start=start.date,
                            end=end.date, quote="AdjClose",
                            provider="yahoo", origin="1970-01-01",
                            compression="m", retclass="zoo")
jwn.prices = get.hist.quote(instrument="jwn", start=start.date,
                            end=end.date, quote="AdjClose",
                            provider="yahoo", origin="1970-01-01",
                            compression="m", retclass="zoo")

# change time indices to class yearmon, which is most appropriate for monthly data
index(alk.prices) = as.yearmon(index(alk.prices))
index(cost.prices) = as.yearmon(index(cost.prices))
index(ffiv.prices) = as.yearmon(index(ffiv.prices))
index(jwn.prices) = as.yearmon(index(jwn.prices))

projectPrices.z = merge(alk.prices,cost.prices,ffiv.prices,
                        jwn.prices)

colnames(projectPrices.z) = asset.names

# create data.frame for downloading
projectPrices.df = as.data.frame(coredata(projectPrices.z))
rownames(projectPrices.df) = as.character(index(projectPrices.z))

#
# compute cc and simple returns
#
projectReturns.z = diff(log(projectPrices.z)) 
# create data.frame for downloading
projectReturns.df = coredata(projectReturns.z)
rownames(projectReturns.df) = as.character(index(projectReturns.z))
projectReturnsSimple.z = exp(projectReturns.z) - 1

#
# plot data
#
my.panel <- function(...) {
  lines(...)
  abline(h=0)
}
plot(projectPrices.z, main="Monthly Prices", col="blue", lwd=2)
plot(projectReturnsSimple.z, main="Monthly Simple Returns", panel=my.panel, col="blue", lwd=2)
plot(projectReturns.z, main="Monthly CC Returns", panel=my.panel, col="blue", lwd=2)

# plot growth of $1 over the five years using PerformanceAnalytics function
# chart.CumReturns
chart.CumReturns(projectReturnsSimple.z, wealth.index=TRUE, legend.loc="topleft", 
                 lwd=2, main="Growth of $1") 

#
# Create matrix of return data and compute pairwise scatterplots
#
ret.mat = coredata(projectReturns.z)


# Note: expand graphics window in Rstudio to see legend

# show a 4 panel plot for vfinx returns - notice the use of drop=FALSE.
# This preserves the column name
fourPanelPlot(projectReturns.z[, "alk", drop=FALSE])
fourPanelPlot(projectReturns.z[, "cost", drop=FALSE])
fourPanelPlot(projectReturns.z[, "ffiv", drop=FALSE])
fourPanelPlot(projectReturns.z[, "jwn", drop=FALSE])

#
# boxplots of returns
#
boxplot(ret.mat, main="Vanguard Returns", col="cornflowerblue")

#
# compute descriptive statistics
#
muhat.vals = colMeans(projectReturns.z)
sd.vals = apply(projectReturns.z, 2, sd)
skew.vals = apply(projectReturns.z, 2, skewness)
ekurt.vals = apply(projectReturns.z, 2, kurtosis)
cov.mat = var(projectReturns.z)
cor.mat = cov2cor(cov.mat)

covhat.vals = cov.mat[lower.tri(cov.mat)]
rhohat.vals = cor.mat[lower.tri(cor.mat)]
names(covhat.vals) <- names(rhohat.vals) <- 
  c("alk,cost", "alk,ffiv", "alk,jwn", 
    "cost,ffiv", "cost,jwn", 
    "ffiv,jwn")
rhohat.vals

# empirical quantiles for VaR calculations
q.vals = apply(projectReturns.z, 2, quantile, prob=c(0.01,0.05))

# display results in a table
stats.mat = rbind(muhat.vals, 
                  sd.vals,
                  skew.vals,
                  ekurt.vals,
                  q.vals)
rownames(stats.mat) = c("Mean", "Std Dev", "Skewness", 
                        "Excess Kurtosis", "1% Quantile", 
                        "5% Quantile")
# print statistics
stats.mat
cov.mat
cor.mat
corrplot.mixed(cor.mat, lower="number", upper="ellipse")

# 7
# compute estimated standard error for mean
nobs = nrow(ret.mat)
nobs
se.muhat = sd.vals/sqrt(nobs)
se.muhat
# show estimates with SE values underneath
rbind(muhat.vals,se.muhat)

# compute approx 95% confidence intervals
mu.lower = muhat.vals - 2*se.muhat
mu.upper = muhat.vals + 2*se.muhat
cbind(mu.lower,mu.upper)

cbind(muhat.vals, se.muhat, mu.lower,mu.upper)

# compute estimated standard errors for variance and sd
se.sigmahat = sd.vals/sqrt(2*nobs)
se.sigmahat

rbind(sd.vals,se.sigmahat)

# compute approx 95% confidence intervals
sigma.lower = sd.vals - 2*se.sigmahat
sigma.upper = sd.vals + 2*se.sigmahat
cbind(sigma.lower,sigma.upper)

cbind(sd.vals,se.sigmahat, sigma.lower,sigma.upper)
#Bootstrap
bootstraped.boot = function(x, idx) {
  muhat = mean(x[idx])
  sigmahat = sd(x[idx])
}
bootstraped.boot.mean = function(x, idx) {
  muhat = mean(x[idx])
}
bootstraped.boot.sd = function(x, idx) {
  sigmahat = sd(x[idx])
}
alk.boot.1000.mean = boot(ret.mat[, "alk"], 
                     statistic=bootstraped.boot.mean, R=1000)
alk.boot.1000.mean
alk.boot.10000.mean = boot(ret.mat[, "alk"], 
                          statistic=bootstraped.boot.mean, R=10000)
alk.boot.10000.mean
boot.ci(alk.boot.1000.mean, conf = 0.95, type = c("norm","perc"))
boot.ci(alk.boot.10000.mean, conf = 0.95, type = c("norm","perc"))

alk.boot.1000.sd = boot(ret.mat[, "alk"], 
                          statistic=bootstraped.boot.sd, R=1000)
alk.boot.1000.sd
alk.boot.10000.sd = boot(ret.mat[, "alk"], 
                      statistic=bootstraped.boot.sd, R=10000)
alk.boot.10000.sd
boot.ci(alk.boot.1000.sd, conf = 0.95, type = c("norm","perc"))
boot.ci(alk.boot.10000.sd, conf = 0.95, type = c("norm","perc"))

alk.boot.1000
alk.boot.10000 = boot(ret.mat[, "alk"], 
                       statistic=bootstraped.boot, R=10000)
alk.boot.10000
boot.ci(alk.boot.1000, conf = 0.95, type = c("norm","perc"))
boot.ci(alk.boot.10000, conf = 0.95, type = c("norm","perc"))

plot(alk.boot.1000)


# SE and 95% confidence values for estimated correlations
n.obs = nrow(projectReturns.z)
SE.rhohat = (1 - rhohat.vals^2)/sqrt(n.obs)
rho.lower = rhohat.vals - 2*SE.rhohat
rho.upper = rhohat.vals + 2*SE.rhohat

rhohat.mat = cbind(rhohat.vals, SE.rhohat, rho.lower, rho.upper)
colnames(rhohat.mat) = c("rho.hat", "SE", "LCL (0.95)", "UCL (0.95)")
rhohat.mat


#
# VaR analysis
#

# function to compute normal and empirical VaR for a matrix of cc returns
## .05
Value.at.Risk = function(x, p=0.05, w=100000, 
                         method=c("normal", "empirical"),
                         return.type=c("cc", "simple")) {
  method=method[1]
  return.type=return.type[1]
  x = as.matrix(x)
  if (method == "normal") {
    q = apply(x, 2, mean) + apply(x, 2, sd)*qnorm(p)
  } else {    
    q = apply(x, 2, quantile, p)
  }
  if (return.type == "simple") {
    VaR = q*w
  } else {
    VaR = (exp(q) - 1)*w
  }
  VaR
}

# compute 5% and 1% normal VaR for all assets
VaR.normal.05 = Value.at.Risk(ret.mat, p=0.05, 
                              method="normal",
                              return.type="cc")
VaR.normal.01 = Value.at.Risk(ret.mat, p=0.01)

# empirical VaR
VaR.empirical.05 = Value.at.Risk(ret.mat, p=0.05, 
                                 method="empirical",
                                 return.type="cc")
VaR.empirical.01 = Value.at.Risk(ret.mat, p=0.01, 
                                 method="empirical",
                                 return.type="cc")

VaR.normal.mat = cbind(VaR.normal.05, VaR.normal.01)
VaR.normal.mat
VaR.empirical.mat = cbind(VaR.empirical.05, VaR.empirical.01)
VaR.empirical.mat

# function to bootstrap VaR for a single asset

## .05
ValueAtRisk.boot = function(x, idx, p=0.05, w=100000,
                            method=c("normal", "empirical"),
                            return.type=c("cc", "simple")) {
  method = method[1]
  return.type = return.type[1]
  if (method == "normal") {
    q = mean(x[idx]) + sd(x[idx])*qnorm(p)
  } else {
    q = quantile(x[idx], p)
  }
  if (return.type == "cc") {
    VaR = (exp(q) - 1)*w
  } else {
    VaR = q*w
  }
  VaR
}
## .01
ValueAtRisk.boot.01 = function(x, idx, p=0.01, w=100000,
                            method=c("normal", "empirical"),
                            return.type=c("cc", "simple")) {
  method = method[1]
  return.type = return.type[1]
  if (method == "normal") {
    q = mean(x[idx]) + sd(x[idx])*qnorm(p)
  } else {
    q = quantile(x[idx], p)
  }
  if (return.type == "cc") {
    VaR = (exp(q) - 1)*w
  } else {
    VaR = q*w
  }
  VaR
}

## .05
ValueAtRisk.boot.emp = function(x, idx, p=0.05, w=100000,
                            method=c("normal", "empirical"),
                            return.type=c("cc", "simple")) {
  method = method[2]
  return.type = return.type[1]
  if (method == "normal") {
    q = mean(x[idx]) + sd(x[idx])*qnorm(p)
  } else {
    q = quantile(x[idx], p)
  }
  if (return.type == "cc") {
    VaR = (exp(q) - 1)*w
  } else {
    VaR = q*w
  }
  VaR
}
## .01
ValueAtRisk.boot.01.emp = function(x, idx, p=0.01, w=100000,
                               method=c("normal", "empirical"),
                               return.type=c("cc", "simple")) {
  method = method[2]
  return.type = return.type[1]
  if (method == "normal") {
    q = mean(x[idx]) + sd(x[idx])*qnorm(p)
  } else {
    q = quantile(x[idx], p)
  }
  if (return.type == "cc") {
    VaR = (exp(q) - 1)*w
  } else {
    VaR = q*w
  }
  VaR
}

computeSEconfintVaR = function(x, p=0.05, w=100000,
                               method=c("normal", "empirical"),
                               return.type=c("cc", "simple")) {
  VaR.boot = boot(x, statistic=ValueAtRisk.boot, R=1000)
  VaR.hat = VaR.boot$t0
  SE.VaR = sd(VaR.boot$t)
  CI.VaR = boot.ci(VaR.boot, conf = 0.95, type="norm")$normal
  CI.VaR = CI.VaR[-1]
  ans = c(VaR.hat, SE.VaR, CI.VaR)
  names(ans) = c("VaR.05", "SE", "LCL (0.95)", "UCL (0.95)")
  return(ans)
}
computeSEconfintVaR.01 = function(x, p=0.01, w=100000,
                               method=c("normal", "empirical"),
                               return.type=c("cc", "simple")) {
  VaR.boot = boot(x, statistic=ValueAtRisk.boot.01, R=1000)
  VaR.hat = VaR.boot$t0
  SE.VaR = sd(VaR.boot$t)
  CI.VaR = boot.ci(VaR.boot, conf = 0.95, type="norm")$normal
  CI.VaR = CI.VaR[-1]
  ans = c(VaR.hat, SE.VaR, CI.VaR)
  names(ans) = c("VaR.01", "SE", "LCL (0.95)", "UCL (0.95)")
  return(ans)
}

computeSEconfintVaR.emp = function(x, p=0.05, w=100000,
                               method=c("normal", "empirical"),
                               return.type=c("cc", "simple")) {
  VaR.boot = boot(x, statistic=ValueAtRisk.boot.emp, R=1000)
  VaR.hat = VaR.boot$t0
  SE.VaR = sd(VaR.boot$t)
  CI.VaR = boot.ci(VaR.boot, conf = 0.95, type="norm")$normal
  CI.VaR = CI.VaR[-1]
  ans = c(VaR.hat, SE.VaR, CI.VaR)
  names(ans) = c("VaR.05", "SE", "LCL (0.95)", "UCL (0.95)")
  return(ans)
}
computeSEconfintVaR.01.emp = function(x, p=0.01, w=100000,
                                  method=c("normal", "empirical"),
                                  return.type=c("cc", "simple")) {
  VaR.boot = boot(x, statistic=ValueAtRisk.boot.01.emp, R=1000)
  VaR.hat = VaR.boot$t0
  SE.VaR = sd(VaR.boot$t)
  CI.VaR = boot.ci(VaR.boot, conf = 0.95, type="norm")$normal
  CI.VaR = CI.VaR[-1]
  ans = c(VaR.hat, SE.VaR, CI.VaR)
  names(ans) = c("VaR.01", "SE", "LCL (0.95)", "UCL (0.95)")
  return(ans)
}


# bootstrap SE and 95% CI for assets
VaR.boot.alk = computeSEconfintVaR(ret.mat[, "alk", drop=FALSE])
VaR.boot.cost = computeSEconfintVaR(ret.mat[, "cost", drop=FALSE])
VaR.boot.ffiv = computeSEconfintVaR(ret.mat[, "ffiv", drop=FALSE])
VaR.boot.jwn = computeSEconfintVaR(ret.mat[, "jwn", drop=FALSE])

VaR.boot.mat = rbind(VaR.boot.alk,
                     VaR.boot.cost,
                     VaR.boot.ffiv,
                     VaR.boot.jwn)
rownames(VaR.boot.mat) = colnames(projectReturns.z)
VaR.boot.mat

# bootstrap SE and 99% CI for assets
VaR.boot.alk.01 = computeSEconfintVaR.01(ret.mat[, "alk", drop=FALSE])
VaR.boot.cost.01 = computeSEconfintVaR.01(ret.mat[, "cost", drop=FALSE])
VaR.boot.ffiv.01 = computeSEconfintVaR.01(ret.mat[, "ffiv", drop=FALSE])
VaR.boot.jwn.01 = computeSEconfintVaR.01(ret.mat[, "jwn", drop=FALSE])

VaR.boot.mat.01 = rbind(VaR.boot.alk.01,
                     VaR.boot.cost.01,
                     VaR.boot.ffiv.01,
                     VaR.boot.jwn.01)
rownames(VaR.boot.mat.01) = colnames(projectReturns.z)
VaR.boot.mat.01

# bootstrap SE and 95% CI for assets
VaR.boot.alk.emp = computeSEconfintVaR.emp(ret.mat[, "alk", drop=FALSE])
VaR.boot.cost.emp = computeSEconfintVaR.emp(ret.mat[, "cost", drop=FALSE])
VaR.boot.ffiv.emp = computeSEconfintVaR.emp(ret.mat[, "ffiv", drop=FALSE])
VaR.boot.jwn.emp = computeSEconfintVaR.emp(ret.mat[, "jwn", drop=FALSE])

VaR.boot.mat.emp = rbind(VaR.boot.alk.emp,
                     VaR.boot.cost.emp,
                     VaR.boot.ffiv.emp,
                     VaR.boot.jwn.emp)
rownames(VaR.boot.mat.emp) = colnames(projectReturns.z)
VaR.boot.mat.emp

# bootstrap SE and 99% CI for assets
VaR.boot.alk.01.emp = computeSEconfintVaR.01.emp(ret.mat[, "alk", drop=FALSE])
VaR.boot.cost.01.emp = computeSEconfintVaR.01.emp(ret.mat[, "cost", drop=FALSE])
VaR.boot.ffiv.01.emp = computeSEconfintVaR.01.emp(ret.mat[, "ffiv", drop=FALSE])
VaR.boot.jwn.01.emp = computeSEconfintVaR.01.emp(ret.mat[, "jwn", drop=FALSE])

VaR.boot.mat.01.emp = rbind(VaR.boot.alk.01.emp,
                        VaR.boot.cost.01.emp,
                        VaR.boot.ffiv.01.emp,
                        VaR.boot.jwn.01.emp)
rownames(VaR.boot.mat.01.emp) = colnames(projectReturns.z)
VaR.boot.mat.01.emp

# annualized VaR numbers
Value.at.Risk.annual = function(x, p=0.05, w=100000, 
                                method=c("normal", "empirical"),
                                return.type=c("cc", "simple")) {
  method=method[1]
  return.type=return.type[1]
  x = as.matrix(x)
  if (method == "normal") {
    q = 12*apply(x, 2, mean) + sqrt(12)*apply(x, 2, sd)*qnorm(p)
  } else {    
    # compute standardized returns and empirical quantile 
    zhat = apply(x, 2, scale)
    qzhat = apply(zhat, 2, quantile, p)
    q = 12*apply(x, 2, mean) + sqrt(12)*apply(x, 2, sd)*qzhat
  }
  if (return.type == "simple") {
    VaR = q*w
  } else {
    VaR = (exp(q) - 1)*w
  }
  VaR
}

VaR.normal.05.annual = Value.at.Risk.annual(ret.mat, p=0.05, 
                                            method="normal",
                                            return.type="cc")
VaR.normal.01.annual = Value.at.Risk.annual(ret.mat, p=0.01, 
                                            method="normal",
                                            return.type="cc")

# empirical VaR
VaR.empirical.05.annual = Value.at.Risk.annual(ret.mat, p=0.05, 
                                               method="empirical",
                                               return.type="cc")
VaR.empirical.01.annual = Value.at.Risk.annual(ret.mat, p=0.01, 
                                               method="empirical",
                                               return.type="cc")

VaR.normal.annual.mat = cbind(VaR.normal.05.annual, VaR.normal.01.annual)
VaR.normal.annual.mat
VaR.empirical.annual.mat = cbind(VaR.empirical.05.annual, VaR.empirical.01.annual)
VaR.empirical.annual.mat


#
# Portfolio theory allowing for short sales
#
rf = 0.005/12

SharpeRatios = (muhat.vals - rf)/sd.vals
SharpeRatios = as.matrix(SharpeRatios, ncol=1)
colnames(SharpeRatios) = "Sharpe"
SharpeRatios

# compute bootstrap standard error and 95% ci for sharpe ratio
# function to bootstrap VaR
sharpeRatio.boot = function(x, idx, risk.free) {
  muhat = mean(x[idx])
  sigmahat = sd(x[idx])
  sharpeRatio = (muhat - risk.free)/sigmahat
  sharpeRatio
}

sharpe.alk.boot = boot(ret.mat[, "alk"], 
                        statistic=sharpeRatio.boot, R=999, risk.free=rf)
sharpe.alk.boot
boot.ci(sharpe.alk.boot, conf = 0.95, type = c("norm","perc"))
plot(sharpe.alk.boot)

computeSEconfintSharpe = function(x, risk.free) {
  Sharpe.boot = boot(x, statistic=sharpeRatio.boot, R=999, risk.free=risk.free)
  Sharpe.hat = Sharpe.boot$t0
  SE.Sharpe = sd(Sharpe.boot$t)
  CI.Sharpe = boot.ci(Sharpe.boot, conf = 0.95, type="norm")$normal
  CI.Sharpe = CI.Sharpe[-1]
  ans = c(Sharpe.hat, SE.Sharpe, CI.Sharpe)
  names(ans) = c("Sharpe", "SE", "LCL (0.95)", "UCL (0.95)")
  return(ans)
}

# bootstrap SE and 95% CI for assets Sharpe Ratio
Sharpe.boot.alk = computeSEconfintSharpe(ret.mat[, "alk", drop=FALSE], risk.free=rf)
Sharpe.boot.cost = computeSEconfintSharpe(ret.mat[, "cost", drop=FALSE], risk.free=rf)
Sharpe.boot.ffiv = computeSEconfintSharpe(ret.mat[, "ffiv", drop=FALSE], risk.free=rf)
Sharpe.boot.jwn = computeSEconfintSharpe(ret.mat[, "jwn", drop=FALSE], risk.free=rf)

Sharpe.boot.mat = rbind(Sharpe.boot.alk,
                        Sharpe.boot.cost,
                        Sharpe.boot.ffiv,
                        Sharpe.boot.jwn)

rownames(Sharpe.boot.mat) = colnames(projectReturns.z)
Sharpe.boot.mat                   



# plot return-risk tradeoff and compute Sharpe ratios
#
## risk free rate

plot(sd.vals, muhat.vals, xlim=c(0, 0.15), ylim=c(0, 0.050),
     ylab="Expected Return", xlab="Standard Deviation",
     cex=2, pch=16, col="cornflowerblue")
points(0, rf, cex=2, pch=16, col="green")
points(sd.vals.vfinx, muhat.val.vfinx, cex=2, pch=16, col="orange")
text(sd.vals.vfinx, muhat.val.vfinx, labels="S&P500", pos=3)
text(0, rf, labels="Rf", pos=3)
text(sd.vals, muhat.vals, labels=colnames(projectReturns.z),
     pos=3)


## compute global minimum variance portfolio
gmin.port <- globalMin.portfolio(muhat.vals, cov.mat)
attributes(gmin.port)
print(gmin.port)
summary(gmin.port, risk.free=rf)
plot(gmin.port, col="cornflowerblue")

# annualized mean and sd
12*gmin.port$er
sqrt(12)*gmin.port$sd

# 1% and 5% VaR
w0=100000
q.05 = (gmin.port$er + gmin.port$sd*qnorm(0.05))
VaR.gmin.05 = (exp(q.05) - 1)*w0
q.01 = (gmin.port$er + gmin.port$sd*qnorm(0.01))
VaR.gmin.01 = (exp(q.01) - 1)*w0
VaR.gmin.05
VaR.gmin.01

## compute global minimum variance portfolio with no short sales
gmin.port.ns <- globalMin.portfolio(muhat.vals, cov.mat, shorts=FALSE)
attributes(gmin.port.ns)
print(gmin.port.ns)
summary(gmin.port.ns, risk.free=rf)
plot(gmin.port.ns, col="cornflowerblue")

# compare weights
wts.mat = cbind(gmin.port$weights,
                gmin.port.ns$weights)
colnames(wts.mat) = c("gmin", "gmin.ns")
wts.mat

# annualized mean and sd
12*gmin.port.ns$er
sqrt(12)*gmin.port.ns$sd

# 1% and 5% VaR
w0=100000
q.05 = (gmin.port.ns$er + gmin.port.ns$sd*qnorm(0.05))
VaR.gmin.ns.05 = (exp(q.05) - 1)*w0
q.01 = (gmin.port.ns$er + gmin.port.ns$sd*qnorm(0.01))
VaR.gmin.ns.01 = (exp(q.01) - 1)*w0
VaR.gmin.ns.05
VaR.gmin.ns.01

VaR.gmin.mat = cbind(c(gmin.port$er*12, gmin.port$sd*sqrt(12), 
                       VaR.gmin.05, VaR.gmin.01),
                     c(gmin.port.ns$er*12, gmin.port.ns$sd*sqrt(12), 
                       VaR.gmin.ns.05, VaR.gmin.ns.01))
rownames(VaR.gmin.mat) = c("Annual Mean", "Annual SD", "VaR.05", "VaR.01")
colnames(VaR.gmin.mat) = c("Shorts", "No Shorts")
VaR.gmin.mat

## efficient portfolio with target return equal to max returns
target.return <- max(muhat.vals)
e.port.max<- efficient.portfolio(muhat.vals, cov.mat, target.return)
e.port.max
summary(e.port.max, risk.free=rf)
plot(e.port.max, col="cornflowerblue")

## compute tangency portfolio
tan.port <- tangency.portfolio(muhat.vals, cov.mat, rf)
tan.port
summary(tan.port, risk.free=rf)
var.tan.port <- tan.port$sd^2
var.tan.port
plot(tan.port, col="cornflowerblue")

# annualized mean and sd
12*tan.port$er
sqrt(12)*tan.port$sd

## compute tangency portfolio with no short sales
tan.port.ns <- tangency.portfolio(muhat.vals, cov.mat, rf, shorts=FALSE)
tan.port.ns
summary(tan.port.ns, risk.free=rf)
var.tan.port.ns <- tan.port.ns$sd^2
var.tan.port.ns
## no question
plot(tan.port.ns, col="cornflowerblue")

# annualized mean and sd
12*tan.port.ns$er
sqrt(12)*tan.port.ns$sd

# compare weights
wts.mat = cbind(tan.port$weights,
                tan.port.ns$weights)
colnames(wts.mat) = c("tan", "tan.ns")

stats.tan.mat = cbind(c(tan.port$er*12, tan.port$sd*sqrt(12)), 
                      c(tan.port.ns$er*12, tan.port.ns$sd*sqrt(12)))
rownames(stats.tan.mat) = c("Annual Mean", "Annual SD")
colnames(stats.tan.mat) = c("Shorts", "No Shorts")
stats.tan.mat


# plot portfolio weights
par(mfrow=c(2,1))
plot(gmin.port, col="cornflowerblue")
plot(tan.port, col="cornflowerblue")
par(mfrow=c(1,1))

par(mfrow=c(2,1))
plot(gmin.port.ns, col="cornflowerblue")
plot(tan.port.ns, col="cornflowerblue")
par(mfrow=c(1,1))
## 23 (c)
par(mfrow=c(2,1))
plot(tan.port, col="cornflowerblue")
plot(tan.port.ns, col="cornflowerblue")
par(mfrow=c(1,1))
### and 19 (a)
par(mfrow=c(2,1))
plot(gmin.port, col="cornflowerblue")
plot(gmin.port.ns, col="cornflowerblue")
par(mfrow=c(1,1))

## compute efficient frontier
ef <- efficient.frontier(muhat.vals, cov.mat, alpha.min=-1, 
                         alpha.max=1.5, nport=20)

plot(ef, plot.assets=TRUE, col="blue", lwd=2)
points(gmin.port$sd, gmin.port$er, col="orange", lwd=5, pch=19)
text(gmin.port$sd, gmin.port$er, labels="global min", pos=4)

## plot efficient frontier allowing shortsales
plot(ef, plot.assets=TRUE, col="blue", lwd=2)
points(gmin.port$sd, gmin.port$er, col="orange", lwd=5, pch=19)
text(gmin.port$sd, gmin.port$er, labels="global min", pos=4)
points(tan.port$sd, tan.port$er, col="red", lwd=5, pch=19)
text(tan.port$sd, tan.port$er, labels="tangency", pos=4)
sr.tan = (tan.port$er - rf)/tan.port$sd
abline(a=rf, b=sr.tan, col="green", lwd=2)
abline(v=0)
points(0, rf, col="green", lwd=5, pch=19)
text(0, rf, labels="rf", pos=4)

## compute efficient frontier not allowing short sales
ef.ns <- efficient.frontier(muhat.vals, cov.mat, alpha.min=0, 
                            alpha.max=1, nport=20, shorts=FALSE)

## plot efficient frontier not allowing shortsales
plot(ef.ns, plot.assets=TRUE, col="blue", lwd=2)
points(gmin.port.ns$sd, gmin.port.ns$er, col="orange", lwd=5, pch=19)
text(gmin.port.ns$sd, gmin.port.ns$er, labels="global min", pos=4)
points(tan.port.ns$sd, tan.port.ns$er, col="red", lwd=5, pch=19)
text(tan.port.ns$sd, tan.port.ns$er, labels="tangency", pos=4)
sr.tan.ns = (tan.port.ns$er - rf)/tan.port.ns$sd
abline(a=rf, b=sr.tan.ns, col="green", lwd=2)
abline(v=0)

# show short sale and no short sale frontiers together
plot(ef, plot.assets=TRUE, col="blue", lwd=2)
points(ef.ns$sd, ef.ns$er, type="b", col="red", lwd=2)

# cost in er for portfolio with sd = 0.02 (0.4%)

# target at sd=0.05
ef$er["port 6"] - ef.ns$er["port 14"]
ef$sd["port 6"]
ef.ns$sd["port 14"]






### Comparing to the S&P 500
index(vfinx.prices) = as.yearmon(index(vfinx.prices))
tan.prices = 0.138*alk.prices+0.765*cost.prices+0.258*ffiv.prices+(-0.161)*jwn.prices
projectPrices.z.24 = merge(vfinx.prices,tan.prices)
 
colnames(projectPrices.z.24) = c("vfinx","tan")
 
# create data.frame for downloading
projectPrices.df.24 = as.data.frame(coredata(projectPrices.z.24))
rownames(projectPrices.df.24) = as.character(index(projectPrices.z.24))

#
# compute cc and simple returns
#
projectReturns.z.24 = diff(log(projectPrices.z.24))
# create data.frame for downloading
projectReturns.df.24 = coredata(projectReturns.z.24)
rownames(projectReturns.df.24) = as.character(index(projectReturns.z.24))
projectReturnsSimple.z.24 = exp(projectReturns.z.24) - 1


# plot growth of $1 over the five years using PerformanceAnalytics function
# chart.CumReturns
chart.CumReturns(projectReturnsSimple.z.24, wealth.index=TRUE, legend.loc="topleft",
                 lwd=2, main="Growth of $1")

muhat.val.vfinx = colMeans(diff(log(vfinx.prices)))
muhat.val.vfinx
sd.vals.vfinx = apply(diff(log(vfinx.prices)), 2, sd)
sd.vals.vfinx
sharpe.vfinx = (muhat.val.vfinx - rf)/sd.vals.vfinx
sharpe.vfinx

x.t = sd.vals.vfinx/tan.port$sd
x.t
x.f = 1 - x.t
x.f

mu.e = rf + x.t*(tan.port$er - rf)
mu.e
sd.e = x.t*tan.port$sd
sd.e
sharpe.e = (mu.e - rf)/sd.e
sharpe.e
stats.mat.24 = rbind(mu.e, sd.e, sharpe.e)
rownames(stats.mat.24) = c("Mean", "Std Dev", "Sharpe's Ratio")
colnames(stats.mat.24)=c("efficient")
stats.mat.24

cbind(stats.mat.24, rbind(muhat.val.vfinx, sd.vals.vfinx, sharpe.vfinx))

# Weigh your tangency portfolio
mytan.port=c(0.138, 0.765, 0.258, -0.161)
tan.port.returns=projectReturnsSimple.z*mytan.port
tan.port.returns2 <- tan.port.returns
tan.port.returns2$port <- rowSums(tan.port.returns[,c("alk", "cost", "ffiv","jwn")])
tan.port.returns2$rf <- rf
port_return=rowSums(tan.port.returns[,c("alk", "cost", "ffiv","jwn")])
port_return
port_return_xts=xts(port_return,index(projectReturnsSimple.z))
port_return_xts
eff.returns.simple <- x.t*port_return_xts+x.f*rf
eff.returns.simple2 <- x.t*tan.port.returns2$port+x.f*tan.port.returns2$rf

projectReturns.z.25 = diff(log(vfinx.prices))
# create data.frame for downloading
projectReturns.df.25 = coredata(projectReturns.z.25)
rownames(projectReturns.df.25) = as.character(index(projectReturns.z.25))
projectReturnsSimple.z.25 = exp(projectReturns.z.25) - 1

plot1 <- cbind(eff.returns.simple2,projectReturnsSimple.z.25)
colnames(plot1)=c("Efficient portfolio","S&P 500")
chart.CumReturns(plot1, wealth.index=TRUE, legend.loc="topleft",
                 lwd=2, main="Growth of $1")
