##============================================================================##
## Load the libraries ----
##============================================================================##
pacman::p_load(tidyverse, raster, biomod2, biogeo, spThin, rasterVis, rgbif, usdm, 
               sf, xlsx, dismo, spocc, viridis)
##============================================================================##

##============================================================================##
## Load the study area ----
##============================================================================##
study_area <- sf::st_read('./Cicer/Final study area.shp')
study_area <- sf::as_Spatial(study_area$geometry)
##============================================================================##


##============================================================================##
## Load the data ----
##============================================================================##
predictors_Peloponnese_current_subset <- subset(current, c('bio_13', 'bio_4',
                                                           'bio_7', 'CECSOL',
                                                           'CLYPPT', 'CRFVOL',
                                                           'ORCDRC', 'SLTPPT'))

predictors_Peloponnese_future_cc_26_subset <- subset(cc_26, c('bio_13', 'bio_4',
                                                              'bio_7', 'CECSOL',
                                                              'CLYPPT', 'CRFVOL',
                                                              'ORCDRC', 'SLTPPT'))

predictors_Peloponnese_future_cc_85_subset <- subset(cc_85, c('bio_13', 'bio_4',
                                                              'bio_7', 'CECSOL',
                                                              'CLYPPT', 'CRFVOL',
                                                              'ORCDRC', 'SLTPPT'))

predictors_Peloponnese_future_bc_26_subset <- subset(bc_26, c('bio_13', 'bio_4',
                                                              'bio_7', 'CECSOL',
                                                              'CLYPPT', 'CRFVOL',
                                                              'ORCDRC', 'SLTPPT'))

predictors_Peloponnese_future_bc_85_subset <- subset(bc_85, c('bio_13', 'bio_4',
                                                              'bio_7', 'CECSOL',
                                                              'CLYPPT', 'CRFVOL',
                                                              'ORCDRC', 'SLTPPT'))

predictors_Peloponnese_future_he_26_subset <- subset(he_26, c('bio_13', 'bio_4',
                                                              'bio_7', 'CECSOL',
                                                              'CLYPPT', 'CRFVOL',
                                                              'ORCDRC', 'SLTPPT'))

predictors_Peloponnese_future_he_85_subset <- subset(he_85, c('bio_13', 'bio_4',
                                                              'bio_7', 'CECSOL',
                                                              'CLYPPT', 'CRFVOL',
                                                              'ORCDRC', 'SLTPPT'))
##============================================================================##



##============================================================================##
## Load the occurrences ----
##============================================================================##
library(maptools)

KMZs <- list.files(path="./Cicer/", pattern="*.kmz", full.names=FALSE)

LonLat <- sapply(KMZs, 
                 function(x) 
                   getKMLcoordinates(kmlfile = unzip(zipfile = paste0("./Cicer/Cicer/", x),
                                                            exdir = "./Cicer"), 
                                     ignoreAltitude = TRUE)[[1]])

colnames(LonLat) <- gsub(pattern = ".kmz", 
                         replacement = "", 
                         x = colnames(LonLat))

rownames(LonLat) <- c("Longitude", "Latitude")

LonLat <- t(LonLat) %>% 
  as.data.frame() 
LonLat

cicer <- LonLat %>%  
  mutate(Specimen = rownames(LonLat),
         Locality = c(rep('Aigio', 2), 
                      rep('Chelmos', 2), 
                      rep('Zireia', 2), 
                      'Chelmos')) %>% 
  as_tibble() %>% 
  dplyr::select(Specimen, everything())

cicer


cyn_sph <- cicer
##============================================================================##


##============================================================================##
## Format the data in the way biomod2 and ecospat want them
##============================================================================##
Dianthus_Data <- BIOMOD_FormatingData(resp.var = rep(1, nrow( cyn_sph ) ), 
                                                              expl.var = predictors_Peloponnese_current_subset,
                                                              resp.xy = cyn_sph[,c('Longitude', 'Latitude')],
                                                              resp.name = "cicer", 
                                                              PA.nb.rep = 100, 
                                                              PA.nb.absences = nrow( cyn_sph ), 
                                                              PA.strategy = 'disk', 
                                                              PA.dist.min  =  5700, 
                                                              PA.dist.max  =  57000)
##============================================================================##




##============================================================================##
## Set the modelling parameters
##============================================================================##
Dianthus_Options <- BIOMOD_ModelingOptions()


##============================================================================##
## Calibrate the simple bivariate models
##============================================================================##
my_ESM_Dianthus <- ecospat.ESM.Modeling( data = Dianthus_Data,
                                models = c('RF', 'ANN', 'CTA'),
                                models.options = Dianthus_Options,
                                NbRunEval = 10,
                                DataSplit = 80,
                                weighting.score = c('TSS'),
                                parallel = T)
##============================================================================##


##============================================================================##


##============================================================================##
## Evaluate and average the simple bivariate models to ESMs
##============================================================================##
### Evaluation and average of simple bivariate models to ESMs


my_ESM_Dianthus_EF <- ecospat.ESM.EnsembleModeling(my_ESM_Dianthus, 
                                                   weighting.score = c("TSS"), 
                                                   models = 'all', 
                                                   threshold = 0)



my_ESM_Dianthus_EF$weights.EF
my_ESM_Dianthus_EF$weights
my_ESM_Dianthus_EF$ESM.evaluations



Ensemble_evaluations <- tibble::as.tibble(my_ESM_Dianthus_EF$ESM.evaluations)
Ensemble_evaluations$model <- factor(Ensemble_evaluations$model)
levels(Ensemble_evaluations$model)

mean(Ensemble_evaluations$TSS)
mean(Ensemble_evaluations$SomersD)
mean(Ensemble_evaluations$Kappa)
mean(Ensemble_evaluations$AUC)
mean(Ensemble_evaluations$MPA)
##============================================================================##


##============================================================================##
## Projection of the simple bivariate models into new space
##============================================================================##
my_ESM_Dianthus_proj_current <- ecospat.ESM.Projection(ESM.modeling.output = my_ESM_Dianthus,
                                            new.env = predictors_Peloponnese_current_subset,
                                            parallel = T)
##============================================================================##


##============================================================================##
## Projection of the calibrated ESMs into new space
##============================================================================##
my_ESM_Dianthus_EFproj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = my_ESM_Dianthus_proj_current,
                                                        ESM.EnsembleModeling.output = my_ESM_Dianthus_EF)
##============================================================================##


##============================================================================##
## Then, we will simply apply the same function with the same parameters as for
## the current conditions for 2050
##============================================================================##

##--------------
## CCSM4 RCP 26
##--------------

### Projection of simple bivariate models  into 2050
my_ESM_Dianthus_proj_2070_cc26 <- ecospat.ESM.Projection(ESM.modeling.output = my_ESM_Dianthus,
                                           new.env = predictors_Peloponnese_future_cc_26_subset,
                                           parallel = T)

### Projection of calibrated ESMs into 2050
my_ESM_Dianthus_EFproj_2070_cc26 <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = my_ESM_Dianthus_proj_2070_cc26,
                                                     ESM.EnsembleModeling.output = my_ESM_Dianthus_EF)
##============================================================================##

##--------------
## CCSM4 RCP 85
##--------------

### Projection of simple bivariate models  into 2050
my_ESM_Dianthus_proj_2070_cc85 <- ecospat.ESM.Projection(ESM.modeling.output = my_ESM_Dianthus,
                                                            new.env = predictors_Peloponnese_future_cc_85_subset,
                                                            parallel = T)

### Projection of calibrated ESMs into 2050
my_ESM_Dianthus_EFproj_2070_cc85 <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = my_ESM_Dianthus_proj_2070_cc85,
                                                                      ESM.EnsembleModeling.output = my_ESM_Dianthus_EF)
##============================================================================##

##--------------
## BCC RCP 26
##--------------

### Projection of simple bivariate models  into 2050
my_ESM_Dianthus_proj_2070_bc26 <- ecospat.ESM.Projection(ESM.modeling.output = my_ESM_Dianthus,
                                                            new.env = predictors_Peloponnese_future_bc_26_subset,
                                                            parallel = T)

### Projection of calibrated ESMs into 2050
my_ESM_Dianthus_EFproj_2070_bc26 <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = my_ESM_Dianthus_proj_2070_bc26,
                                                                      ESM.EnsembleModeling.output = my_ESM_Dianthus_EF)
##============================================================================##

##--------------
## BCC RCP 85
##--------------

### Projection of simple bivariate models  into 2050
my_ESM_Dianthus_proj_2070_bc85 <- ecospat.ESM.Projection(ESM.modeling.output = my_ESM_Dianthus,
                                                            new.env = predictors_Peloponnese_future_bc_85_subset,
                                                            parallel = T)

### Projection of calibrated ESMs into 2050
my_ESM_Dianthus_EFproj_2070_bc85 <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = my_ESM_Dianthus_proj_2070_bc85,
                                                                      ESM.EnsembleModeling.output = my_ESM_Dianthus_EF)
##============================================================================##

##----------------
## HadGEM2 RCP 26
##----------------

### Projection of simple bivariate models  into 2050
my_ESM_Dianthus_proj_2070_he26 <- ecospat.ESM.Projection(ESM.modeling.output = my_ESM_Dianthus,
                                                            new.env = predictors_Peloponnese_future_he_26_subset,
                                                            parallel = T)

### Projection of calibrated ESMs into 2050
my_ESM_Dianthus_EFproj_2070_he26 <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = my_ESM_Dianthus_proj_2070_he26,
                                                                      ESM.EnsembleModeling.output = my_ESM_Dianthus_EF)
##============================================================================##

##----------------
## HadGEM2 RCP 85
##----------------

### Projection of simple bivariate models  into 2050
my_ESM_Dianthus_proj_2070_he85 <- ecospat.ESM.Projection(ESM.modeling.output = my_ESM_Dianthus,
                                                            new.env = predictors_Peloponnese_future_he_85_subset,
                                                            parallel = T)

### Projection of calibrated ESMs into 2050
my_ESM_Dianthus_EFproj_2070_he85 <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = my_ESM_Dianthus_proj_2070_he85,
                                                                      ESM.EnsembleModeling.output = my_ESM_Dianthus_EF)
##============================================================================##


##============================================================================##
## Create current clamping mask
##============================================================================##
# 
current.list_Dianthus <- list.files(path = 'ESM.BIOMOD.output.cicer',
                           pattern = "_ClampingMask\\.grd$", all.files = T,
                           full.names = T, recursive = T, include.dirs = T)


current <- Filter(function(x) grepl("current", x), current.list_Dianthus)
cc_26 <- Filter(function(x) grepl("cc_26", x), current.list_Dianthus)
cc_85 <- Filter(function(x) grepl("cc_85", x), current.list_Dianthus)
bc_26 <- Filter(function(x) grepl("bc_26", x), current.list_Dianthus)
bc_85 <- Filter(function(x) grepl("bc_85", x), current.list_Dianthus)
he_26 <- Filter(function(x) grepl("he_26", x), current.list_Dianthus)
he_85 <- Filter(function(x) grepl("he_85", x), current.list_Dianthus)




current_stack_Dianthus <- stack(current)
cc_26_stack_Dianthus <- stack(cc_26)
cc_85_stack_Dianthus <- stack(cc_85)
bc_26_stack_Dianthus <- stack(bc_26)
bc_85_stack_Dianthus <- stack(bc_85)
he_26_stack_Dianthus <- stack(he_26)
he_85_stack_Dianthus <- stack(he_85)

current_Clamping_Mask_Dianthus <- mean(current_stack_Dianthus)
cc_26_Clamping_Mask_Dianthus <- mean(cc_26_stack_Dianthus)
cc_85_Clamping_Mask_Dianthus <- mean(cc_85_stack_Dianthus)
bc_26_Clamping_Mask_Dianthus <- mean(bc_26_stack_Dianthus)
bc_85_Clamping_Mask_Dianthus <- mean(bc_85_stack_Dianthus)
he_26_Clamping_Mask_Dianthus <- mean(he_26_stack_Dianthus)
he_85_Clamping_Mask_Dianthus <- mean(he_85_stack_Dianthus)

levelplot(current_Clamping_Mask_Dianthus,
          main = 'Clamping mask current Cicer graecum',
          col.regions = colorRampPalette(c('grey90', 'yellow4', 'green4'))(100))

levelplot(cc_26_Clamping_Mask_Dianthus,
          main = 'Clamping mask CCSM4 RCP2.6 Cicer graecum',
          col.regions = colorRampPalette(c('grey90', 'yellow4', 'green4'))(100))

levelplot(cc_85_Clamping_Mask_Dianthus,
          main = 'Clamping mask CCSM4 RCP8.5 Cicer graecum',
          col.regions = colorRampPalette(c('grey90', 'yellow4', 'green4'))(100))

levelplot(bc_26_Clamping_Mask_Dianthus,
          main = 'Clamping mask BCC RCP2.6 Cicer graecum',
          col.regions = colorRampPalette(c('grey90', 'yellow4', 'green4'))(100))

levelplot(bc_85_Clamping_Mask_Dianthus,
          main = 'Clamping mask BCC RCP8.5 Cicer graecum',
          col.regions = colorRampPalette(c('grey90', 'yellow4', 'green4'))(100))

levelplot(he_26_Clamping_Mask_Dianthus,
          main = 'Clamping mask HadGEM2 RCP2.6 Cicer graecum',
          col.regions = colorRampPalette(c('grey90', 'yellow4', 'green4'))(100))

levelplot(he_85_Clamping_Mask_Dianthus,
          main = 'Clamping mask HadGEM2 RCP8.5 Cicer graecum',
          col.regions = colorRampPalette(c('grey90', 'yellow4', 'green4'))(100))


## Load the predictions -----
current.list_Dianthus <- list.files(path = '/ESM.BIOMOD.output.cicer/',
                                    pattern = "_ESM.BIOMOD.{,2}\\.grd$", all.files = T,
                                    full.names = T, recursive = T, include.dirs = T)


current <- Filter(function(x) grepl("current", x), current.list_Dianthus)
cc_26 <- Filter(function(x) grepl("cc_26", x), current.list_Dianthus)
cc_85 <- Filter(function(x) grepl("cc_85", x), current.list_Dianthus)
bc_26 <- Filter(function(x) grepl("bc_26", x), current.list_Dianthus)
bc_85 <- Filter(function(x) grepl("bc_85", x), current.list_Dianthus)
he_26 <- Filter(function(x) grepl("he_26", x), current.list_Dianthus)
he_85 <- Filter(function(x) grepl("he_85", x), current.list_Dianthus)

current_stack_Dianthus <- stack(current)
cc_26_stack_Dianthus <- stack(cc_26)
cc_85_stack_Dianthus <- stack(cc_85)
bc_26_stack_Dianthus <- stack(bc_26)
bc_85_stack_Dianthus <- stack(bc_85)
he_26_stack_Dianthus <- stack(he_26)
he_85_stack_Dianthus <- stack(he_85)

current_Dianthus <- mean(current_stack_Dianthus)
cc_26_Dianthus <- mean(cc_26_stack_Dianthus)
cc_85_Dianthus <- mean(cc_85_stack_Dianthus)
bc_26_Dianthus <- mean(bc_26_stack_Dianthus)
bc_85_Dianthus <- mean(bc_85_stack_Dianthus)
he_26_Dianthus <- mean(he_26_stack_Dianthus)
he_85_Dianthus <- mean(he_85_stack_Dianthus)
##============================================================================##


##============================================================================##
## Plot
##============================================================================##
cyno_sph_occ <- cyn_sph %>% dplyr::select(Longitude, Latitude)

rasterVis::levelplot(my_ESM_Dianthus_EFproj_current$EF,
                     main = 'Cicer graecum\n Current ensemble projections',
                     col.regions = colorRampPalette(c('grey90', 'yellow4', 'green4'))(100)) +
  layer(sp.points(SpatialPoints(cyno_sph_occ), pch = 20, col = 1))

rasterVis::levelplot(my_ESM_Dianthus_EFproj_2070_cc26,
          main = 'Cicer graecum ensemble projections\nin 2070 with cc26',
          col.regions = colorRampPalette(c('grey90', 'yellow4', 'green4'))(100)) +
  layer(sp.points(SpatialPoints(cyno_sph_occ), pch = 20, col = 1))
##============================================================================##


##============================================================================##
# Construct binary maps from an optimal cutoff and find the range change
##============================================================================##

EF_current_cutoff_Dianthus <- ecospat.mpa(current_Dianthus, 
                                          cyno_sph_occ, perc = 1)
EF_future_cutoff_Dianthus_cc26 <- ecospat.mpa(cc_26_Dianthus, 
                                              cyno_sph_occ, perc = 1)
EF_future_cutoff_Dianthus_cc85 <- ecospat.mpa(cc_85_Dianthus, 
                                              cyno_sph_occ, perc = 1)
EF_future_cutoff_Dianthus_bc26 <- ecospat.mpa(bc_26_Dianthus, 
                                              cyno_sph_occ, perc = 1)
EF_future_cutoff_Dianthus_bc85 <- ecospat.mpa(bc_85_Dianthus, 
                                              cyno_sph_occ, perc = 1)
EF_future_cutoff_Dianthus_he26 <- ecospat.mpa(he_26_Dianthus, 
                                              cyno_sph_occ, perc = 1)
EF_future_cutoff_Dianthus_he85 <- ecospat.mpa(he_85_Dianthus, 
                                              cyno_sph_occ, perc = 1)


binary_current_Dianthus <- ecospat::ecospat.binary.model(current_Dianthus, 
                                                            EF_current_cutoff_Dianthus)
binary_future_Dianthus_cc26 <- ecospat::ecospat.binary.model(cc_26_Dianthus, 
                                                                EF_future_cutoff_Dianthus_cc26)
binary_future_Dianthus_cc85 <- ecospat::ecospat.binary.model(cc_85_Dianthus, 
                                                                EF_future_cutoff_Dianthus_cc85)
binary_future_Dianthus_bc26 <- ecospat::ecospat.binary.model(bc_26_Dianthus, 
                                                                EF_future_cutoff_Dianthus_bc26)
binary_future_Dianthus_bc85 <- ecospat::ecospat.binary.model(bc_85_Dianthus, 
                                                                EF_future_cutoff_Dianthus_bc85)
binary_future_Dianthus_he26 <- ecospat::ecospat.binary.model(he_26_Dianthus, 
                                                                EF_future_cutoff_Dianthus_he26)
binary_future_Dianthus_he85 <- ecospat::ecospat.binary.model(he_85_Dianthus, 
                                                                EF_future_cutoff_Dianthus_he85)
##============================================================================##


##============================================================================##
## Conservative approach
##============================================================================##
## the suitability of all cells showing variable values not experienced during
## the model training (values greater than zero in the clamping mask) was set to
## zero
##============================================================================##
binary_current_Dianthus[current_Clamping_Mask_Dianthus > 0] = 0
binary_future_Dianthus_cc26[cc_26_Clamping_Mask_Dianthus > 0] = 0
binary_future_Dianthus_cc85[cc_85_Clamping_Mask_Dianthus > 0] = 0
binary_future_Dianthus_bc26[bc_26_Clamping_Mask_Dianthus > 0] = 0
binary_future_Dianthus_bc85[bc_85_Clamping_Mask_Dianthus > 0] = 0
binary_future_Dianthus_he26[he_26_Clamping_Mask_Dianthus > 0] = 0
binary_future_Dianthus_he85[he_85_Clamping_Mask_Dianthus > 0] = 0
##============================================================================##


##============================================================================##
## SRC current -> 2050
##============================================================================##

##--------------------------------
## First do it for all the metrics
##--------------------------------
SRC_current_2050_cc26_Dianthus <- BIOMOD_RangeSize( binary_current_Dianthus,
                                           binary_future_Dianthus_cc26 )

SRC_current_2050_cc85_Dianthus <- BIOMOD_RangeSize( binary_current_Dianthus,
                                           binary_future_Dianthus_cc85 )

SRC_current_2050_bc26_Dianthus <- BIOMOD_RangeSize( binary_current_Dianthus,
                                           binary_future_Dianthus_bc26 )

SRC_current_2050_bc85_Dianthus <- BIOMOD_RangeSize( binary_current_Dianthus,
                                           binary_future_Dianthus_bc85 )

SRC_current_2050_he26_Dianthus <- BIOMOD_RangeSize( binary_current_Dianthus,
                                           binary_future_Dianthus_he26 )

SRC_current_2050_he85_Dianthus <- BIOMOD_RangeSize( binary_current_Dianthus,
                                           binary_future_Dianthus_he85 )



SRC_current_2050_cc26_Dianthus$Compt.By.Models
SRC_current_2050_cc85_Dianthus$Compt.By.Models
SRC_current_2050_bc26_Dianthus$Compt.By.Models
SRC_current_2050_bc85_Dianthus$Compt.By.Models
SRC_current_2050_he26_Dianthus$Compt.By.Models
SRC_current_2050_he85_Dianthus$Compt.By.Models
##============================================================================##



##============================================================================##
## Mask by environmental area
##============================================================================##

##--------------------
## On the whole island
##--------------------

Dianthus_src_map_cc26 <- stack(SRC_current_2050_cc26_Dianthus$Diff.By.Pixel)
Dianthus_src_map_cc85 <- stack(SRC_current_2050_cc85_Dianthus$Diff.By.Pixel)
Dianthus_src_map_bc26 <- stack(SRC_current_2050_bc26_Dianthus$Diff.By.Pixel)
Dianthus_src_map_bc85 <- stack(SRC_current_2050_bc85_Dianthus$Diff.By.Pixel)
Dianthus_src_map_he26 <- stack(SRC_current_2050_he26_Dianthus$Diff.By.Pixel)
Dianthus_src_map_he85 <- stack(SRC_current_2050_he85_Dianthus$Diff.By.Pixel)


my.at <- seq(-2.5,1.5,1) 
myColorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     labels=c("lost", "pres", "abs","gain"), ## labels
                     at=my.at[-1]-0.5 ## where to print labels
                   ))

ensemble_src_map <- stack(Dianthus_src_map_cc26,
                          Dianthus_src_map_cc85,
                          Dianthus_src_map_bc26,
                          Dianthus_src_map_bc85,
                          Dianthus_src_map_he26,
                          Dianthus_src_map_he85)

names(ensemble_src_map) <- c('CCSM4 RCP 26', 'CCSM4 RCP 85', 
                             'BCC RCP 26', 'BCC RCP 85', 
                             'HadGEM2 RCP 26', 'HadGEM2 RCP 85')

rasterVis::levelplot( ensemble_src_map,
                      main = "Cicer graecum range change",
                      colorkey = myColorkey,
                      col.regions = colorRampPalette(c('green4', 'grey90', 'black'))(100))
layout = c(2,3) 
##============================================================================##
