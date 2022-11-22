# (c) Author: Sophia Burke
# Institute for the Study of Earth Oceans and Space
# University of New Hampshire
# Original File Created : June 5, 2020
# File Modified: 11/15/2022
# contact: sophia.burke@unh.edu; sophieaburke@gmail.com

# This is the source code for supplemental calculations presented in 
# Burke et al. 2022 manuscript entitled: Connecting Methane Ebullition Flux to
# Thaw Pond size using Unpiloted Aerial Systems. In particular: Median Daily Average Flux and Total 
# Daily Average Flux across an eight day moving window that is before, centered 
# around, and after each UAS flight date. We settled on using Median Daily 
# Average Flux across a centered eight day moving window based on the 
# statistical significance of the flux windows against water polygon area. 
# Those statistical tests are included. This code also shows how Cumulative 
# flux was calculated.

# This code also includes how we determined eight days was a good window size.

# Note: The #comments in this code pertain to the code BELOW the comment.

####################################################################################
 
### set working directory *** Change this to suit your own directory ***##########
WD <- 'C:\\Users\\saj82\\OneDrive - USNH\\RSE_ThawPondPaper\\R_Stordalen_UASImagery'
setwd(WD)

# load in necessary packages (what functions are they needed for?)
library(plyr) # (ddply function)
library(dplyr) #(distinct function)
library(zoo) #(rollmedian and rollsum functions)
#####

#load in the UAS data file & format it ###########################################

# StordalenMire_ThawPond_UASPolygons_2014to2018.csv
# 12 columns 
# FlightDate - YYYY-MM-DD 
# Field_ID - pond identifier that matching the fieldbooks (Field_letter)
# Burkeetal2019_ID - pond identifier that matches what was used in 
      # Burke et al.(2019) (letter)
# UAStype -  either 'FWing' or 'Quad'
    # 'FWing' -- Fixed Wing UAS
    # 'Quad'  -- Quadcopter UAS
# PolygonType - either 'pondedge' or 'water'
# PondType - numerical value: 1-4, represents the pond type described in 
      # Burke et al. (2019)
# area_m2 -  the area of the drawn polygon in m2 (drawn in QGIS)
# edge_m - the length of the perimeter (edge) of the drawn polygon in m 
      # (drawn in QGIS)
# edge.area -  the ratio between edge/area, unitless
# Median8dC_BubFlux - the calculated median daily ebullitive flux over an 
      # eight day centered moving window (see below for explaination of how 
      # this is calcualted) in mg CH4 m-2 d-1
# Cumulative_BubFlux - the cumulative ebullitive emission per unit area at the 
      # end of the sampling season for that row in the matrix
# TotPrec_8dpreflight - total precipiation recorded (mm) in the eight days 
    # leading up to and including the Flight date. This was calculated from 
    # 10min to 30min frequency meteorological data measured by ICOS-Sweden 
    # (Integrated Carbon Observation System)at Stordalen Mire, Sweden
# see ReadMe_StordalenMire_ThawPond_UASPolygons_2014to2018.txt file for more 
    # details.

RSEdataframe_all <- read.csv(
          ".\\StordalenMire_ThawPond_UASPolygons_2014to2018.csv", header = TRUE)
# format flight date to a POSIXct object
RSEdataframe_all$FlightDate <- as.Date(RSEdataframe_all$FlightDate,
                                       format = "%Y-%m-%d", tz = "CET")
# create a column that is just 'month'
RSEdataframe_all$month      <- as.numeric(strftime(RSEdataframe_all$FlightDate,
                                                   "%m", tz = "CET"))
# create a column that is the numerical year YYYY
RSEdataframe_all$year       <- as.numeric(strftime(RSEdataframe_all$FlightDate,
                                                   "%Y", tz = "CET"))
# create a column that is the day of year (doy) based on the flight date
RSEdataframe_all$doy        <- as.numeric(strftime(RSEdataframe_all$FlightDate,
                                                   "%j", tz = "CET"))
# create a column that pastes together the Burkeetal2019_ID (pond) and year for 
# that row for matching purposes 
RSEdataframe_all$PondYr     <- paste(RSEdataframe_all$Burkeetal2019_ID,
                                     RSEdataframe_all$year,sep = "")
# same thing but for Date & pond, for matching purposes
RSEdataframe_all$Date_Pond  <- paste(RSEdataframe_all$FlightDate,
                                     RSEdataframe_all$Burkeetal2019_ID, 
                                     sep = "_")

# How did we choose eight days for the size of the moving window?##################
# Since the Quadcopter was flown multiple times a season, we looked 
  # specifically at the time difference between Quadcopter flights 

# lets therefore subset RSEdataframe_all so it only contains quadcopter data
QuadFlightDiff <- subset(RSEdataframe_all,RSEdataframe_all$UAStype == "Quad")
# lets further subset this so as to include only unique flight dates
QuadFlightDiff<- distinct(QuadFlightDiff,QuadFlightDiff$FlightDate)
# give the column a proper name
colnames(QuadFlightDiff) <- c("FlightDate")
# create a column containing just 'year' to separate out the flight dates
QuadFlightDiff$year <- strftime(QuadFlightDiff$FlightDate,"%Y", tz = "CET")
# order the matrix so that the dates are in the correct order
QuadFlightDiff <- QuadFlightDiff[order(QuadFlightDiff$FlightDate),]
# calculate the different in flight dates within each sampling season
QuadFlightDiff$DiffinTime <- ave(as.numeric(QuadFlightDiff$FlightDate), 
                                 QuadFlightDiff$year, 
                                 FUN = function(x) c(NA,diff(x)))
# calculate the average difference in flight dates across 2016,2017 and 2018
AvgTimeDiff_betweenflights <- mean(QuadFlightDiff$DiffinTime, na.rm = TRUE)
# the answer, 7.545 days was rounded up to Eight for moving windows used in 
# this project.
#####

#load in Daily Average Bubble Flux file & format it ##########################################

# StordalenMire_ThawPond_DailyAvgBubbleFlux_2012-2018.csv
# 4 columns 
# Date - YYYY-MM-DD 
# Field_ID - pond identifier that matching the fieldbooks
# Burkeetal2019_ID - pond identifier that matches what was used in 
      # Burke et al. (2019)
# DailyAvgBubbleFlux - daily average ebullitive flux (DAPf - daily 
      # average pond flux)
# per unit area in each pond in mg CH4 m-2 d-1 
# see ReadMe file for more details.

#read in the file  
BubFlux_2012to2018 <- read.csv(
      ".\\StordalenMire_ThawPond_DailyAvgBubbleFlux_2012-2018.csv", header = TRUE)
# format the date column to a POSIXct object
BubFlux_2012to2018$Date   <- as.Date(BubFlux_2012to2018$Date,
                                        format = "%Y-%m-%d", tz = "CET")
# create a collumn that is numerical year YYYY
BubFlux_2012to2018$Year   <- as.numeric(strftime(BubFlux_2012to2018$Date,
                                                 "%Y", tz = "CET"))
# create a column that pastes together the Burkeetal2019_ID (aka pond) and 
  # year for that row for matching purposes
BubFlux_2012to2018$PondYr <- paste(BubFlux_2012to2018$Burkeetal2019_ID, 
                                   BubFlux_2012to2018$Year, sep = "")
# create a column that pastes together the Date column and the Burkeetal2019_ID
  # (aka pond) for that row for matching purposes
BubFlux_2012to2018$Date_Pond  <- paste(BubFlux_2012to2018$Date,
                                       BubFlux_2012to2018$Burkeetal2019_ID, 
                                       sep = "_")
# lets round the daily flux values to the third decimal place so that its 
  # within the accuracy of the measurment (see ReadMe)
# Daily flux values are accurate to three decimal places
BubFlux_2012to2018[,4] <- round(BubFlux_2012to2018[,4],digits = 3) 

#####

## Calculate Cumulative Flux & Eight Day Median & Total Daily Average Flux ######
#Create a dataframe that contains cumulative flux per year (mg CH4 m-2), per 
  # pond, using the PondYr column as the grouping factor
CumulativeFlux <- ddply(BubFlux_2012to2018,.(PondYr), summarize,
                        CumulativeFlux = sum(DailyAvgBubbleFlux, na.rm = TRUE))
# round the Cumulative flux numbers to within the measurement accuracy of the
 # daily measurments. Daily flux values are accurate to three decimal places
CumulativeFlux[,2] <- round(CumulativeFlux[,2],digits = 3)

# BubFlux_2012to2018 contains all the Days that I have measured DAPf for.
# in order to calculate the rollmedians and rollsums correctly I need to ensure
# all days June 1-Sept 30 are present in the matrix (whether I have a daily 
# flux value or not for each day) I'm applying a rollmedian too.
# There are occasions in the middle of the season where I don't have a DAPf 
# (due to a toss sample etc). but I still need that date to be present so that
# the rollmedian can be preformed on the correct dates.

# Create a data.frame that is the correct size:
# what is the correct size? There are eight ponds with DAPf and the easiest 
# thing to do is to create a matrix that has enough rows to fill with dates 
# ranging from June 1 - Sept 30 for each of the seven sampling seasons (though
# we don't have ebullition data to fill this entire matrix).

# create an empty data frame that is 6,832 rows and 5 columns 
ExpandBubFlux <- data.frame(matrix(data=NA,nrow = 6832,ncol = 5))
#name the collumns
names(ExpandBubFlux) <- c("Date","Burkeetal2019_ID","DailyAvgBubbleFlux",
                          "Date_Pond","Year")
#create a column that contains the Burkeetal2019_ID, with each pond ID repeated
# 854 times. Why 854? June 1 - Sept 30 is 122 days, 122 * 7 sampling 
    # seasons = 854
ponds <-  c(rep("A",854),rep("B",854),rep("C",854),rep("D",854),rep("E",854),
            rep("F",854),rep("G",854),rep("H",854))
# repeate June 1 - Sept 30 for each field season 2012 to 2018
dates <- c(seq(as.Date("2012-06-01",format = "%Y-%m-%d",tz = "CET"),by = "day",
               length = 122),
           seq(as.Date("2013-06-01",format = "%Y-%m-%d",tz = "CET"),by = "day", 
               length = 122),
           seq(as.Date("2014-06-01",format = "%Y-%m-%d",tz = "CET"),by = "day",
               length = 122),
           seq(as.Date("2015-06-01",format = "%Y-%m-%d",tz = "CET"),by = "day",
               length = 122),
           seq(as.Date("2016-06-01",format = "%Y-%m-%d",tz = "CET"),by = "day",
               length = 122),
           seq(as.Date("2017-06-01",format = "%Y-%m-%d",tz = "CET"),by = "day",
               length = 122),
           seq(as.Date("2018-06-01",format = "%Y-%m-%d",tz = "CET"),by = "day", 
               length = 122))
# repeat the dates matrix eight times so each pond (of which there are 8 in the 
    # BubFlux_2012to2018 file) will have a date (June 1 - Sept 30)
Dates <- c(rep(dates,8))
# format the dates I just created
ExpandBubFlux$Date <- as.Date(Dates,format = "%Y-%m-%d", tz = "CET")
# fill the Burkeetal2019_ID column with the values in 'ponds'
ExpandBubFlux$Burkeetal2019_ID <- ponds
# create an identifier column that is a combo of Date and pond that is unique to
# each row
  # this column will match the previously created Date_Pond column in 
  # BubFlux_2012to2018
ExpandBubFlux$Date_Pond <- paste(ExpandBubFlux$Date,
                                 ExpandBubFlux$Burkeetal2019_ID, sep = "_")
# match the DailyAvgBubbleFlux listed in BubFlux_2012to2018 to the ExpandBubFlux 
# matrix
  # to the corresponding Date_Pond row.
ExpandBubFlux$DailyAvgBubbleFlux <- BubFlux_2012to2018$DailyAvgBubbleFlux[
  match(ExpandBubFlux$Date_Pond,BubFlux_2012to2018$Date_Pond)]
# create a column that designates the year of the data in each row
ExpandBubFlux$Year <- strftime(ExpandBubFlux$Date, "%Y", tz = "CET")
# create an identifier column that is a combo of the pond ID and the year column
ExpandBubFlux$PondYr <- paste(ExpandBubFlux$Burkeetal2019_ID,
                              ExpandBubFlux$Year,sep = " ")

# calculate a MEDIAN daily average flux (mg CH4 m-2 d-1) across an eight day 
# moving window that is:
# A) the eight days leading up to and including the flight date
ExpandBubFlux$Median8dR_BubFlux <- ave(ExpandBubFlux$DailyAvgBubbleFlux,
                                       ExpandBubFlux$PondYr,
                                       FUN = function(x) rollmedian(x,k = 8,
                                                    fill = NA,align = "right"))
# B) the three days leading up to the flight date, the flight date, and the four
  # days following the flight date
ExpandBubFlux$Median8dC_BubFlux <- ave(ExpandBubFlux$DailyAvgBubbleFlux,
                                       ExpandBubFlux$PondYr,
                                       FUN = function(x) rollmedian(x,k = 8,
                                                    fill = NA,align = "center"))
# C) the eigth days following and including the flight date.
ExpandBubFlux$Median8dL_BubFlux <- ave(ExpandBubFlux$DailyAvgBubbleFlux,
                                   ExpandBubFlux$PondYr,
                                   FUN = function(x) rollmedian(x,k = 8,
                                                    fill = NA,align = "left"))

# calculate a TOTAL bubble flux (mg CH4 m-2) across an eight day moving window 
# that is:
# A) the eight days leading up to and including the flight date
ExpandBubFlux$Sum8dR_BubFlux <- ave(ExpandBubFlux$DailyAvgBubbleFlux,
                                    ExpandBubFlux$PondYr,
                                    FUN = function(x) rollsum(x,k = 8,fill = NA,
                                                              align = "right"))
# B) the three days leading up to the flight date, the flight date, and the four 
 # days following the flight date
ExpandBubFlux$Sum8dC_BubFlux <- ave(ExpandBubFlux$DailyAvgBubbleFlux,
                                    ExpandBubFlux$PondYr,
                                    FUN = function(x) rollsum(x,k = 8,fill = NA,
                                                              align = "center"))
# C) the eigth days following and including the flight date.
ExpandBubFlux$Sum8dL_BubFlux <- ave(ExpandBubFlux$DailyAvgBubbleFlux,
                                    ExpandBubFlux$PondYr,
                                    FUN = function(x) rollsum(x,k = 8,fill = NA,
                                                              align = "left"))
##### 

#match above Cumulative Flux, Median 8d flux, Total 8d Flux to the matching... ########## 
  # ...PondYr in RSEdataframe_all 
#
# the numbers in Cumulativeflux2 collumn should match those in the existing 
  # CumulativeFlux column in the RSEdataframe_all matrix
RSEdataframe_all$CumulativeFlux2 <- CumulativeFlux$CumulativeFlux[
  match(RSEdataframe_all$PondYr,CumulativeFlux$PondYr)]

# add these flux calculations to the UAS matrix by matching PondYr
RSEdataframe_all$Median8dR_BubFlux <- ExpandBubFlux$Median8dR_BubFlux[
  match(RSEdataframe_all$Date_Pond,ExpandBubFlux$Date_Pond)]

## the numbers in Median8dC_BubFlux2 collumn should match those in the existing 
# Median8dC_BubFlux column in the RSEdataframe_all matrix
RSEdataframe_all$Median8dC_BubFlux2 <- ExpandBubFlux$Median8dC_BubFlux[
  match(RSEdataframe_all$Date_Pond,ExpandBubFlux$Date_Pond)]

RSEdataframe_all$Median8dL_BubFlux <- ExpandBubFlux$Median8dL_BubFlux[
  match(RSEdataframe_all$Date_Pond,ExpandBubFlux$Date_Pond)]
# match the Total fluxes calculated above to the correct rows in the 
  # RSEdataframe_all by matching Date_Pond
RSEdataframe_all$Sum8dR_BubFlux <- ExpandBubFlux$Sum8dR_BubFlux[
  match(RSEdataframe_all$Date_Pond,ExpandBubFlux$Date_Pond)]

RSEdataframe_all$Sum8dC_BubFlux <- ExpandBubFlux$Sum8dC_BubFlux[
  match(RSEdataframe_all$Date_Pond,ExpandBubFlux$Date_Pond)]

RSEdataframe_all$Sum8dL_BubFlux <- ExpandBubFlux$Sum8dL_BubFlux[
  match(RSEdataframe_all$Date_Pond,ExpandBubFlux$Date_Pond)]

#####

## Statistical Analysis ##########################################################
# Now to test the statistical significance of the different 8 day moving windows
# of daily average ebullitive flux, we looked specifically at Quadcopter data 
# and the water polygons.

# create a subset of just water polygon data from Quadcopter imagery
Quad_Water_PondPolygons <- subset(RSEdataframe_all,
                                  RSEdataframe_all$UAStype == "Quad" & 
                                    RSEdataframe_all$PolygonType == "water")

#make a simple plot to see the distribution of the data...
plot(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, 
     ylim = c(0,500), pch = 21)

# perform a Kendall Correlation test using cor.test 
cor.test(Quad_Water_PondPolygons$area_m2,
                              Quad_Water_PondPolygons$Median8dC_BubFlux2,
                              method = "kendall" )
# AF8dC: p value .02 tau 0.14
#

cor.test(Quad_Water_PondPolygons$area_m2,
                              Quad_Water_PondPolygons$Median8dR_BubFlux,
                              method = "kendall" )
# AF8dR: p value .07 NOT SIGNIFICANT
#

cor.test(Quad_Water_PondPolygons$area_m2,
                              Quad_Water_PondPolygons$Median8dL_BubFlux, 
                              method = "kendall" )
# AF8dL: .08 NOT SIGNIFICANT
#
cor.test(Quad_Water_PondPolygons$area_m2,
                              Quad_Water_PondPolygons$Sum8dC_BubFlux, 
                              method = "kendall" )
# AF8dC: p value .07 tau 0.11
#                          
cor.test(Quad_Water_PondPolygons$area_m2,
                              Quad_Water_PondPolygons$Sum8dR_BubFlux, 
                              method = "kendall" )
# AF8dR: p value .21 NOT SIGNIFICANT       
#             
cor.test(Quad_Water_PondPolygons$area_m2,
                              Quad_Water_PondPolygons$Sum8dL_BubFlux, 
                              method = "kendall" )
# AF8dL: .13 NOT SIGNIFICANT           
#             

## Median8dC_BubFlux was the most significant with a p value of .02, so we chose
# to focus our further tests on Median8dC_BubFlux
