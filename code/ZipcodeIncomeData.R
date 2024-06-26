rm(list=ls())
library(sf)
library(readxl)
library(data.table)
library(geosphere)
library(svMisc)
library(parallel)
setwd("/Users/srikur/Downloads/tl_2020_us_zcta520/")

# read all columns as character
distances <- fread("min_zipcode_distances.csv", header = TRUE, colClasses = "character")
# get the top 5% largest distances
top_5_percent <- distances[order(-as.numeric(distances$min_distance)),][1:round(0.05*nrow(distances)),]

race_data_2020 <- fread("ACSDP5Y2020.DP05-Data.csv", header = TRUE)[-1,]
income_data_2020 <- fread("ACSST5Y2020.S1901-Data.csv", header = TRUE)[-1,]

# DP05_0001E = total population
# DP05_0002E = number male
# DP05_0003E = number female
# DP05_0037E = number white
# DP05_0038E = number black
# DP05_0039E = number native
# DP05_0044E = number asian
# DP05_0052E = number pacific
# DP05_0057E = number other
# DP05_0058E = number mixed
# DP05_0071E = number hispanic

race_data_2020 <- race_data_2020[, .(NAME, DP05_0001E, DP05_0002E, DP05_0003E, DP05_0037E, DP05_0038E, 
                                    DP05_0039E, DP05_0044E, DP05_0052E, DP05_0057E, DP05_0058E, DP05_0071E)]
names(race_data_2020) <- c("zipcode", "total_pop", "num_male", "num_female", "num_white", "num_black", 
                           "num_native", "num_asian", "num_pacific", "num_other", "num_mixed", "num_hispanic")
race_data_2020[, zipcode := gsub("ZCTA5 ", "", zipcode)]
race_data_2020[, zipcode := as.character(zipcode)]
# make all columns numeric except zipcode
race_data_2020[, c("total_pop", "num_male", "num_female", "num_white", "num_black", 
                   "num_native", "num_asian", "num_pacific", "num_other", "num_mixed", "num_hispanic") := lapply(.SD, as.numeric), 
                   .SDcols = c("total_pop", "num_male", "num_female", "num_white", "num_black", 
                               "num_native", "num_asian", "num_pacific", "num_other", "num_mixed", "num_hispanic")]


# S1901_C01_012E = median income

income_data_2020 <- income_data_2020[, .(NAME, S1901_C01_012E)]
names(income_data_2020) <- c("zipcode", "median_income")
income_data_2020[, zipcode := gsub("ZCTA5 ", "", zipcode)]
income_data_2020[, zipcode := as.character(zipcode)]

# get median_income rows where median_income has a non-digit character
# error_income = income_data_2020[!grepl("^[0-9]+$", median_income)]
# top_5_percent$zipcode[which(top_5_percent$zipcode %in% error_income$zipcode)]

income_data_2020 <- income_data_2020[!median_income == "-",]
income_data_2020[, median_income := gsub("[^0-9]", "", median_income)]


# Filter race_data_2020 to only include zipcodes that are also in income_data_2020
race_data_2020 <- race_data_2020[zipcode %in% income_data_2020$zipcode,]

# Merge
zipcode_data <- merge(income_data_2020, race_data_2020, by = "zipcode")
zipcode_data[, top_5_percent := ifelse(zipcode %in% top_5_percent$zipcode, 1, 0)]

# Save the data
setwd("/Users/srikur/Documents/GitHub/medicine/data/")
# fwrite(zipcode_data, "zipcode_dataset.csv")

# Graph 1: Zipcode by Income, separated into 5 buckets
# - income < 25k
# - 25k <= income < 50k
# - 50k <= income < 75k
# - 75k <= income < 100k
# - income >= 100k

graph1_data <- zipcode_data[, .(zipcode, median_income, top_5_percent)]
graph1_data[, median_income := as.numeric(median_income)]
# create income buckets with count
graph1_data <- graph1_data[, .(count = .N), by = .(zipcode, median_income, top_5_percent)]
graph1_data <- graph1_data[, income_bucket := cut(median_income, breaks = c(0, 25000, 50000, 75000, 100000, Inf), 
                                                  labels = c("< 25k", "25k - 50k", "50k - 75k", "75k - 100k", "> 100k"))]
graph1_data <- graph1_data[, .(count = sum(count)), by = .(zipcode, income_bucket, top_5_percent)]
# group percentages for each income bucket out of total zipcodes for either top 5% or not top 5%
num_top5 <- sum(graph1_data$top_5_percent)
num_not_top5 <- nrow(graph1_data) - num_top5
graph1_data <- graph1_data[, percentage := count/num_top5, by = .(income_bucket, top_5_percent)]
graph1_data <- graph1_data[, percentage := ifelse(top_5_percent == 0, count/num_not_top5, percentage)]
# sum percentages for each income bucket
graph1_data <- graph1_data[, .(percentage = sum(percentage)), by = .(income_bucket, top_5_percent)]

library(ggplot2)

ggplot(graph1_data, aes(x = income_bucket, y = percentage, fill = factor(top_5_percent))) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(title = "Median Income by Zip Code", x = "Income Bucket", y = "Percent of Zipcodes") + 
  scale_fill_discrete(name = "", labels = c("Trial Site Zip Codes", "Top 5% Zip Codes Furthest From Trial Site")) +
  theme(legend.position = "bottom")

# Graph 2: Demographic Percentages by Zipcode for Top 5% Furthest, US Average, and Trial Site

# calculate percentages for each demographic category
graph2_data <- zipcode_data[, .(zipcode, total_pop, top_5_percent, num_white, num_black, num_native, num_pacific, num_asian, num_hispanic)]
graph2_data[, num_white := num_white/total_pop]
graph2_data[, num_black := num_black/total_pop]
graph2_data[, num_native := num_native/total_pop]
graph2_data[, num_asian := num_asian/total_pop]
graph2_data[, num_pacific := num_pacific/total_pop]
graph2_data[, num_hispanic := num_hispanic/total_pop]

# combine native and pacific into other
graph2_data[, num_other := num_native + num_pacific]

# group by top_5_percent and average percentages for each demographic category
graph2_data <- graph2_data[, .(num_white = mean(num_white), num_black = mean(num_black), num_other = mean(num_other), 
                               num_asian = mean(num_asian), num_hispanic = mean(num_hispanic)), by = .(top_5_percent)]
names(graph2_data)[1] <- "group"
# set group to "Trial Site" for 0 and "Top 5%" for 1
graph2_data <- graph2_data[, group := ifelse(group == 0, "Trial Site", "Top 5%")]

# add US average demographic percentages from zip code data
us_avg_demographics <- zipcode_data[, .(num_white = sum(num_white)/sum(total_pop), num_black = sum(num_black)/sum(total_pop), 
                                        num_native = sum(num_native)/sum(total_pop), num_pacific = sum(num_pacific)/sum(total_pop), 
                                        num_asian = sum(num_asian)/sum(total_pop), num_hispanic = sum(num_hispanic)/sum(total_pop))]
us_avg_demographics <- us_avg_demographics[, num_other := num_native + num_pacific]
us_avg_demographics <- cbind(us_avg_demographics, group = "US Average")
us_avg_demographics[, num_pacific := NULL]
us_avg_demographics[, num_native := NULL]

# combine all demographic data
graph2_data <- rbind(graph2_data, us_avg_demographics)

# melt data for ggplot
graph2_data <- melt(graph2_data, id.vars = c("group"), measure.vars = c("num_white", "num_black", "num_other", "num_asian", "num_hispanic"))

# change "variable" values to "White", "Black", "Other", "Asian", "Hispanic"
graph2_data <- graph2_data[, variable := ifelse(variable == "num_white", "White", 
                                                ifelse(variable == "num_black", "Black", 
                                                       ifelse(variable == "num_other", "Other", 
                                                              ifelse(variable == "num_asian", "Asian", "Hispanic"))))]

# create ggplot
ggplot(graph2_data, aes(x = variable, y = value, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(title = "Demographic Percentages by Zip Code", x = "Demographic Category", y = "Percent of Population") + 
  scale_fill_discrete(name = "", labels = c("Trial Site Zip Codes", "Top 5% Furthest Zip Codes", "US Average")) +
  theme(legend.position = "bottom")
