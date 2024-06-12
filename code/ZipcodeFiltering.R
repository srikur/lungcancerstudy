rm(list=ls())
library(sf)
library(readxl)
library(data.table)
library(geosphere)
library(svMisc)
library(parallel)
setwd("/Users/srikur/Downloads/tl_2020_us_zcta520/")

# Load the shapefile
shapefile <- st_read("tl_2020_us_zcta520.shp")

lung_cancer_data <- data.table(read_excel("/Users/srikur/Documents/GitHub/medicine/data/lung_cancer_studies_final.xlsx"))
# get Zipcodes variable, split on "|" and unnest
zipcodes <- lung_cancer_data$Zipcodes
zipcodes <- strsplit(zipcodes, "|", fixed = TRUE)
zipcodes <- unlist(zipcodes)
zipcodes <- zipcodes[!is.na(zipcodes)]
zipcodes <- as.character(zipcodes)

# boston_zips = c("02108", "02109", "02110", "02111", "02112", "02113", "02114", "02115", "02116", "02118", "02119", "02120", "02121", "02122", "02124", "02125", "02126", "02127", "02128", "02129", "02130", "02131", "02132", "02134", "02135", "02136", "02151", "02152", "02163", "02199", "02203", "02210", "02215", "02467")
# unique(zipcodes[which(zipcodes %in% boston_zips)])

# Filter the shapefile to only include the zipcodes in the lung cancer data
zipcode_shapefile <- shapefile[shapefile$ZCTA5CE20 %in% zipcodes,]

# Save the filtered shapefile
st_write(zipcode_shapefile, "zipcode_filtered.shp", append = FALSE)

#------ Calculations ------

# For each ZCTA5CE20 in the total shapefile, get the centroid of its polygon
polygon_centroids <- st_centroid(shapefile)
# convert to data.table with variables: ZCTA5CE20, X, Y
polygon_centroid_dt <- data.table(ZCTA5CE20 = shapefile$ZCTA5CE20, 
                                  X = st_coordinates(polygon_centroids)[,"X"], 
                                  Y = st_coordinates(polygon_centroids)[,"Y"])

# For each ZCTA5CE20 in the filtered shapefile, get the centroid of its polygon
filtered_polygon_centroids <- st_centroid(zipcode_shapefile)
# convert to data.table with variables: ZCTA5CE20, X, Y
filtered_polygon_centroid_dt <- data.table(ZCTA5CE20 = zipcode_shapefile$ZCTA5CE20, 
                                          X = st_coordinates(filtered_polygon_centroids)[,"X"], 
                                          Y = st_coordinates(filtered_polygon_centroids)[,"Y"])


# Goal: Get the zipcodes in shapefile that are > 500 miles away from any zipcodes in the filtered shapefile
# For each zipcode in the shapefile, calculate the distance to all zipcodes in the filtered shapefile
# Store the minimum distance for each zipcode in the shapefile

# Processing the data in parallel

# Function to calculate the distance betweena zipcode in the shapefile and all zipcodes in the filtered shapefile
process_zcta <- function(current_zipcode) {
  # current_zipcode <- polygon_centroid_dt[i]
  current_zipcode_coords <- c(current_zipcode$X, current_zipcode$Y)
  
  distances <- sapply(1:nrow(filtered_polygon_centroid_dt), function(j) {
    filtered_zipcode <- filtered_polygon_centroid_dt[j]
    filtered_zipcode_coords <- c(filtered_zipcode$X, filtered_zipcode$Y)
    distm(current_zipcode_coords, filtered_zipcode_coords, fun = distHaversine)
  })
  
  min(distances)
}
# Create a data.table to store the minimum distance for each zipcode in the shapefile
min_distances <- data.table(ZCTA5CE20 = shapefile$ZCTA5CE20, min_distance = rep(0, nrow(shapefile)))

# For each zipcode in the shapefile, calculate the distance to all zipcodes in the filtered shapefile
# Store the minimum distance for each zipcode in the shapefile
# Use mclapply to process the data in parallel with 16 cores
results <- mclapply(1:nrow(polygon_centroid_dt), function(i) {
  process_zcta(polygon_centroid_dt[i])
}, mc.cores = 12)

# Update the min_distances data.table with the calculated minimum distances
min_distances$min_distance <- unlist(results)

names(min_distances) <- c("zipcode", "min_distance")

# get quantiles: 0.05, 0.25, 0.5, 0.75, 0.95, 0.99
quantiles <- quantile(min_distances$min_distance, c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99))

# Save the min_distances data.table
setwd("/Users/srikur/Documents/GitHub/medicine/data/")
fwrite(min_distances, "min_zipcode_distances.csv")

# filter shapefile based on min_distances > median
filtered_shapefile <- shapefile[shapefile$ZCTA5CE20 %in% min_distances$zipcode[min_distances$min_distance > quantiles[5]],]

# Save the filtered shapefile
st_write(filtered_shapefile, "large_distance_zcta_shapefile.shp", append = FALSE)