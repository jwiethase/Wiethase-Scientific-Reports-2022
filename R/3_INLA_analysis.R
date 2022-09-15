library(raster)
library(INLA)

### Load data
### ----------------------------------------------------------------------------
deg200 <- stack("data/BG200m_deg.tif")
load("data/BG200m_raw_spdf.RData")
annRain <- read.csv("data/annRain.csv")

### ----------------------------------------------------------------------------
time <- 1:nlayers(deg200)

# Compute the overall slope of degradation, accounting for rainfall in a given year
slowfun <- function(y) {
  if(all(is.na(y))){
    NA
  } else {
    m = lm(y ~ time + annRain$mean); summary(m)$coefficients[2]   # get time
  }
}

deg200_rate <- calc(deg200, slowfun)

### INLA model predicting rate from various covariates:
### ----------------------------------------------------------------------------

## Add the rate to the raw_spdf:
deg200_rate_spdf <- as(deg200_rate, "SpatialPointsDataFrame")
deg200_rate_spdf <- deg200_rate_spdf[!is.na(deg200_rate_spdf@data$layer), ]

xy_raw  <- apply(round(coordinates(raw_spdf), 6), 1, paste, collapse = "_")
xy_rate <- apply(round(coordinates(deg200_rate_spdf), 6), 1, paste, collapse = "_")

raw_spdf <- raw_spdf[order(xy_raw), ]
deg200_rate_spdf <- deg200_rate_spdf[xy_rate %in% xy_raw, ]
deg200_rate_spdf <- deg200_rate_spdf[order(xy_rate[xy_rate %in% xy_raw]), ]
raw_spdf$deg_rate <- deg200_rate_spdf$layer * 10  # Small decimal numbers will use more memory to get to integer. Also easier to read, and index is arbitrary anyways

## Transform and scale covariates:
raw_spdf$lpop <- log(raw_spdf$pop + 1)  # log(x+1) for right skewed data containing zeros
raw_spdf$ltlu <- log(raw_spdf$tlu + 1)  
raw_spdf$pop_s <- c(scale(raw_spdf$lpop))
raw_spdf$rain_s <- c(scale(raw_spdf$rain))
raw_spdf$tlu_s <- c(scale(raw_spdf$ltlu))

## Make a mesh:
max.edge = 0.05  # ~6km at equator. Max. triangle length 
mesh <- inla.mesh.2d(
      loc=coordinates(raw_spdf),
      offset = c(0.05, 0.5),    
      max.edge=c(max.edge, max.edge*4),  
      cutoff=max.edge/2)
# plot(mesh)

# and the spde:
A = inla.spde.make.A(mesh=mesh, loc=data.matrix(coordinates(raw_spdf))) 

# Covariates:
Xcov = data.frame(intercept=1, 
                  pop_s = raw_spdf$pop_s,
                  rain_s = raw_spdf$rain_s,
                  tlu_s = raw_spdf$tlu_s,
                  CA1 = (raw_spdf$CA == 1) * 1,  # CCRO. logical to numeric. In GLM, R will do this automatically
                  CA2 = (raw_spdf$CA == 2) * 1,  # NP
                  CA3 = (raw_spdf$CA == 3) * 1)  # WMA
Xcov = as.matrix(Xcov)
colnames(Xcov)
stck <- inla.stack(tag='est',     # stack is INLA data structure
                   data=list(deg_rate=raw_spdf$deg_rate),
                   effects=list(
                         s=1:mesh$n,
                         Xcov=Xcov),
                   A = list(A, 1)
)

make_model = function(prior.median.sd, prior.median.range, var){
  spde = inla.spde2.pcmatern(mesh, 
                             prior.range = c(prior.median.range,var),  
                             prior.sigma = c(prior.median.sd, var), 
                             constr = TRUE)
  
  formula = deg_rate ~ -1 + Xcov + f(s, model=spde)  
  prior.median.gaus.sd = 0.2
  family = 'gaussian'
  control.family = list(hyper = list(prec = list(
        prior = "pc.prec", fixed = FALSE, param = c(prior.median.gaus.sd,0.5))))
  
  ## make some linear combinations, for effect plots:
  pop_lc <- inla.make.lincombs(Xcov1 = rep(1, 100),
                               Xcov2 = seq(min(raw_spdf$pop_s), max(raw_spdf$pop_s), len = 100))
  names(pop_lc) <- paste0("pop", 1:100)
  rain_lc <- inla.make.lincombs(Xcov1 = rep(1, 100),
                                Xcov3 = seq(min(raw_spdf$rain_s), max(raw_spdf$rain_s), len = 100))
  names(rain_lc) <- paste0("rain", 1:100)
  tlu_lc <- inla.make.lincombs(Xcov1 = rep(1, 100),
                               Xcov4 = seq(min(raw_spdf$tlu_s), max(raw_spdf$tlu_s), len = 100))
  names(tlu_lc) <- paste0("tlu", 1:100)
  CA_lc <- inla.make.lincombs(Xcov1 = rep(1, 4),
                              Xcov5 = c(0,1,0,0),
                              Xcov6 = c(0,0,1,0),
                              Xcov7 = c(0,0,0,1))
  names(CA_lc) <- paste0("CA", 1:4)
  all_lc <- c(pop_lc, rain_lc, tlu_lc, CA_lc)
  
  ## fit model:
  res <- inla(formula, data=inla.stack.data(stck),
              control.predictor=list(A = inla.stack.A(stck), compute=TRUE),
              family = family,
              lincomb = all_lc,
              control.inla = list(int.strategy='eb'),
              verbose=FALSE)
return(res)
}

res <- make_model(0.4, 0.12, 0.5) 

summary(res)

for (i in 1:length(res$marginals.fixed)) {
   tmp = inla.tmarginal(function(x) x, res$marginals.fixed[[i]]) 
   plot(tmp, type = "l", xlab = paste("Fixed effect marginal", i, ":", colnames(Xcov)[i]), ylab = "Density")
   abline(v = 0, lty = 2)
}

unscale <- function(x, scale.params = sc.p) {
   return((x * scale.params$'scaled:scale') + scale.params$'scaled:center')
}

# Function to add letter to top left corner
put.fig.letter <- function(label, x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  coords <- c(0.09,0.9)
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, cex = 1.5, ...)
}

#### FIGURE 7
par(mfrow = c(2,2))

plot(raw_spdf$pop + 1, raw_spdf$deg_rate, axes = FALSE, bty = "l",
     xlab = as.expression(bquote("Population Density (" ~km^2~ ")")),
     ylab = "Slope of bare ground scores", pch = 20, col = rgb(0,0,0,0.05),
     log = "x")
box(bty = "o")
sc.p <- attributes(scale(raw_spdf$lpop))
xs <- exp(unscale(seq(min(raw_spdf$pop_s), max(raw_spdf$pop_s), len = 100)))
polygon(c(xs, rev(xs)), c(res$summary.lincomb.derived[grep("pop", rownames(res$summary.lincomb.derived)), "0.025quant"],
                          rev(res$summary.lincomb.derived[grep("pop", rownames(res$summary.lincomb.derived)), "0.975quant"])),
        border = NA, col = rgb(0.7, 0, 0.1, 0.5))
lines(xs, res$summary.lincomb.derived$`0.5quant`[grep("pop", rownames(res$summary.lincomb.derived))], lwd = 2, col = rgb(0.7, 0, 0.1))
axis(2)
axis(1, at = c(1, 6, 11, 51, 101, 501), labels = c(0, 5, 10, 50, 100, 500))
put.fig.letter(label = 'A)', font = 2)

plot(raw_spdf$rain, raw_spdf$deg_rate, 
     xlab = "Rainfall (mm)",
     ylab = "Slope of bare ground scores", pch = 20, col = rgb(0,0,0,0.05))
box(bty = "o")
sc.p <- attributes(scale(raw_spdf$rain))
xs <- unscale(seq(min(raw_spdf$rain_s), max(raw_spdf$rain_s), len = 100))
polygon(c(xs, rev(xs)), c(res$summary.lincomb.derived[grep("rain", rownames(res$summary.lincomb.derived)), "0.025quant"],
                          rev(res$summary.lincomb.derived[grep("rain", rownames(res$summary.lincomb.derived)), "0.975quant"])),
        border = NA, col = rgb(0.7, 0, 0.1, 0.5))
lines(xs, res$summary.lincomb.derived$`0.5quant`[grep("rain", rownames(res$summary.lincomb.derived))], lwd = 2, col = rgb(0.7, 0, 0.1))
put.fig.letter(label = 'B)', font = 2)

plot(raw_spdf$tlu + 1, raw_spdf$deg_rate, axes = FALSE, bty = "l",
     xlab = as.expression(bquote("Livestock Density (TLU /" ~km^2~ ")")),
     ylab = "Slope of bare ground scores", pch = 20, col = rgb(0,0,0,0.05),
     log = "x")
box(bty = "o")
sc.p <- attributes(scale(raw_spdf$ltlu))
xs <- exp(unscale(seq(min(raw_spdf$tlu_s), max(raw_spdf$tlu_s), len = 100)))
polygon(c(xs, rev(xs)), c(res$summary.lincomb.derived[grep("tlu", rownames(res$summary.lincomb.derived)), "0.025quant"],
                          rev(res$summary.lincomb.derived[grep("tlu", rownames(res$summary.lincomb.derived)), "0.975quant"])),
        border = NA, col = rgb(0.7, 0, 0.1, 0.5))
lines(xs, res$summary.lincomb.derived$`0.5quant`[grep("tlu", rownames(res$summary.lincomb.derived))], lwd = 2, col = rgb(0.7, 0, 0.1))
axis(2)
axis(1, at = c(1, 6, 11, 51, 101, 501), labels = c(0, 5, 10, 50, 100, 500))
put.fig.letter(label = 'C)', font = 2)

raw_spdf$CA <- factor(raw_spdf$CA, levels = c(0, 1, 3, 2), labels = c("None", "CCRO", "WMA", "NP "))
boxplot(raw_spdf$deg_rate ~ raw_spdf$CA_labelled, xaxt = "n",
        xlab = "Land use designation",
        ylab = "Slope of bare ground scores")
axis(1, at = 1:4, labels = c("Village", "CCRO", "NP", "WMA"))
for(i in 1:4){
   lines(rep(i, 2), res$summary.lincomb.derived[grep("CA", rownames(res$summary.lincomb.derived)), c(4, 6)][i,],
         col =  rgb(0.7, 0, 0.1))
   points(i, res$summary.lincomb.derived[grep("CA", rownames(res$summary.lincomb.derived)), 5][i],
          col =  rgb(0.7, 0, 0.1), pch = 20)
}
put.fig.letter(label = 'D)', font = 2)

par(mfrow = c(1,1))

