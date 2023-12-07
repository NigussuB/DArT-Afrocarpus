
### Working on the DArT Data

## CONVERTING COORDINATES FROM DMS TO DD
#instal packages and load libraries
install.packages("parzer")
install.packages("readxl")
install.packages("dplyr")
install.packages("writexl")
install.packages("oce")
require(parzer)
require(tidyverse)
require(sf)
require(leaflet)
library(dplyr)
library(parzer)
library(readxl)
library(writexl)
library(oce)

#load data
coords = readxl::read_excel("C:/Users/iDesire Computer/Documents/R/DArT/coord_data.xlsx", sheet = 1)

head(coords)
#convert the coordinates that were in UTM into decimal degrees
#library(oce)
coords_dd <- utm2lonlat(coords$Easting, coords$Northing, zone = 37, hemisphere = "N", km = FALSE)

head(coords_dd)
#change the resulting file into a data frame
coord_1 <- data.frame(coords_dd)
head(coord_1)

# write the file / expport it to excel

write_xlsx(coord_1, "coord_data.xlsx" )

## Reading DArT Files into a genlight Object using read.dart()

library(dartR)
