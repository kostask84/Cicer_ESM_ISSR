##============================================================================##
## Load the libraries ----
##============================================================================##
pacman::p_load(tidyverse, extrafont, fossil, betapart)
loadfonts(device="win")
##============================================================================##


##============================================================================##
## Create geographical distance matrix ----
##============================================================================##
aegean_geo_data <- readxl::read_excel('./Cicer/Cicer genetic distances.xlsx', 
                                      sheet = 1) %>% 
  dplyr::select(Longitude:Latitude)

coordinates(aegean_geo_data) <- c("Longitude","Latitude")
crs(aegean_geo_data) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
aegean_geo_data_mx <- as.matrix(earth.dist(aegean_geo_data))


genetic_distance <- readxl::read_excel('./Cicer/Cicer genetic distances.xlsx', 
                                       sheet = 3) %>% 
  as.matrix()
genetic_distance
##============================================================================##


##============================================================================##
## Mantel correlograms ----
##============================================================================##

opar <- par()

tiff('Mantel correlograms.tif', 
     units = "cm", 
     width = 20, 
     height = 20, 
     res = 600)

par(mfrow = c(1,1), family = "Times New Roman")

plot(vegan::mantel.correlog(genetic_distance, aegean_geo_data_mx,
                            nperm = 9999,
                            # cutoff = F,
                            r.type = 'spearman'))


dev.off()

par(opar)
##============================================================================##

