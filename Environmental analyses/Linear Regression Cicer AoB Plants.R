##============================================================================##
## Load the packages ----
##============================================================================##
pacman::p_load(dismo, rgbif, raster, dplyr, viridis, tibble, magrittr, broom,
                 twidlr, leaps, relaimpo, gvlma, MASS)
##============================================================================##

##============================================================================##
## Load the data -----
##============================================================================##
soil_data <- readRDS('Soil data stacked.rds')
current <- readRDS('Current climate data stacked.rds')
current <- stack(current, soil_data)
##============================================================================##


##============================================================================##
## Load the file containing the coords ----
##============================================================================##
cicer <- readxl::read_excel('./Cicer/Cicer genetic distances.xlsx')

all_cretan <- cicer
##============================================================================##


##============================================================================##
## Extract the environmental variables for each taxon ----
##============================================================================##
coordinates(all_cretan) <- c("Longitude","Latitude")

pts.clim_cretan <- raster::extract(current, all_cretan, method = "bilinear")

petromarula.clim_cretan <- data.frame(cbind(coordinates(all_cretan), 
                                            pts.clim_cretan, 
                                            all_cretan@data))
##============================================================================##


##============================================================================##
## Convert to tibble ----
##============================================================================##
petromarula.clim_cretan <- as_tibble(petromarula.clim_cretan)
##============================================================================##


##============================================================================##
## Correlation  ----
##============================================================================##
usdm::vifcor(petromarula.clim_cretan %>% 
               as.data.frame() %>% 
               dplyr::select(aridityIndexThornthwaite:SNDPPT), 
             th = 0.7)

cicer_mlr <- petromarula.clim_cretan %>% 
  dplyr::select(Taxon, H, uH, aridityIndexThornthwaite, bio_3, 
                CLYPPT, ORCDRC)
##============================================================================##


##============================================================================##
## Normality  ----
##============================================================================##
normality <- lapply(cicer_mlr[,2:7], shapiro.test)
normality
##============================================================================##


##============================================================================##
## Log-transform  ----
##============================================================================##
cicer_mlr_log <- lapply(cicer_mlr[,3:7], as.numeric) %>% 
  as_tibble() %>% 
  log10() %>% 
  mutate(Taxon = cicer_mlr$Taxon,
         H = cicer_mlr$H,
         uH = cicer_mlr$uH) %>% 
  dplyr::select(Taxon, H, uH, everything()) 
##============================================================================##


##============================================================================##
## Relative weights function  ----
##============================================================================##
relweights <- function(fit,...){
  R <- cor(fit$model)
  nvar <- ncol(R)
  rxx <- R[2:nvar, 2:nvar]
  rxy <- R[2:nvar, 1]
  svd <- eigen(rxx)
  evec <- svd$vectors
  ev <- svd$values
  delta <- diag(sqrt(ev))
  lambda <- evec %*% delta %*% t(evec)
  lambdasq <- lambda ^ 2
  beta <- solve(lambda) %*% rxy
  rsquare <- colSums(beta ^ 2)
  rawwgt <- lambdasq %*% beta ^ 2
  import <- (rawwgt / rsquare) * 100
  import <- as.data.frame(import)
  row.names(import) <- names(fit$model[2:nvar])
  names(import) <- "Weights"
  import <- import[order(import),1, drop=FALSE]
  dotchart(import$Weights, labels=row.names(import),
           xlab="% of R-Square", pch=19,
           main="Relative Importance of Predictor Variables",
           sub=paste("Total R-Square=", round(rsquare, digits=3)),
           ...)
  return(import)
}
##============================================================================##


##============================================================================##
## Best subsets regression  ----
##============================================================================##
bes.tot <- regsubsets(H ~ aridityIndexThornthwaite + bio_3 + CLYPPT + ORCDRC, 
                      data = cicer_mlr_log, 
                      nbest = 4)

plot(bes.tot,
     scale = "adjr2")

bes.tot <- regsubsets(uH ~ aridityIndexThornthwaite + bio_3 + CLYPPT + ORCDRC, 
                      data = cicer_mlr_log, 
                      nbest = 4)

plot(bes.tot,
     scale = "adjr2")
##============================================================================##


##============================================================================##
## Final regression  ----
##============================================================================##
tot <- lm(H ~ aridityIndexThornthwaite + bio_3 + CLYPPT, 
          data = cicer_mlr_log)  

summary(tot)
##============================================================================##


##============================================================================##
## Relative weights  ----
##============================================================================##
relweights(tot, 
           col = "blue")
##============================================================================##


##============================================================================##
## Check if the IVs have any problem  ----
##============================================================================##
sqrt(car::vif(tot)) > 2
##============================================================================##


##============================================================================##
## Bonferroni correction  ----
##============================================================================##
summary(tot) 
p <- c(0.00095, 0.02620, 0.01434) 
alp <- 0.05 
p.adjust(p, method="bonferroni", n=length(p)) 
##============================================================================##




