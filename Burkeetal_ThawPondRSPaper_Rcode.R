###################################################################################
# (c) Author: Sophia Burke
# Institute for the Study of Earth Oceans and Space
# University of New Hampshire
# Original File Created : May 28, 2020
# File Modified : April 16, 2024
# contact: sophia.burke@unh.edu; ruth.varner@unh.edu

# This is the source code for figures and statistics presented in 
# Burke et al. 2024 manuscript entitled: Connecting Methane Flux to
# Thaw Pond size using Unpiloted Aerial Systems.

# Other citations of note: Burke, S.A., Wik, M., Lang, A., Contosta, A.R., 
#     Palace, M., Crill, P.M., Varner, R.K., 2019. Long‐Term Measurements of 
#     Methane Ebullition From Thaw Ponds. J. Geophys. Res. Biogeosci. 
#     2018JG004786. https://doi.org/10.1029/2018JG004786

# the data loaded into this file is available on the EMERGE Github in a public
# repository: 

# The input data file is called: /Stordalen_ThawPond_UASPolygons_2014to2018.csv  
# The columns are:
#     FlightDate: YYYY-MM-DD, time zone: CET
#     Field_ID: original pond identifier used in field notes
#     Burkeetl2019_ID: Pond ID used in figures in this paper, following the ID 
#                       used in Burke et al. (2019)
#     UAStype: Either 'FWing' for Fixed Wing, or 'Quad' for Quadcopter
#     PolygonType: Either 'pondedge' for Pond Depression Polygon, or 'water' for Water
#                  Polygon
#     PondType: Either 1,2,3 or 4. See Burke et al. (2019) for complete description 
#               of each type
#     area_m2: area of the polygon (m2) -- calculated in QGIS
#     edge_m: perimeter (edge) length of the polygon (m) -- calculated in QGIS
#     edge.area: edge divided by area of the polygon -- calculated in QGIS
#     Median8dC_BubFlux: Median daily bubbly flux from each thaw pond calculated 
#           across a eight day window centered on each flight date.
#     Cumulative_BubFlux: Total CH4 flux per m2 emitted in the corresponding 
#           field season of each UAS flight
#     TotPrec_8dpreflight : Total precipitation (mm) measured in the eight days 
#           leading up to and including each flight date

# See ReadMe_StordalenMire_ThawPond_UASPolygons_2014to2018.txt for more details

# Note: any #comment# that appears in this code refers to the line of code below it.

# Note: Reference to Pond Edge Area (PEA) in this code is the same as Pond Depression
# Area as in the manuscript. 

####################################################################################

## load in necessary packages (what functions they provide)##########################
library(tidyr) # (separate function)
library(plyr) # (ddply function)
library(dunn.test) #(dunn.test function)
library(dplyr) #(distinct function)
library(agricolae) #(kruskal function)
library(PMCMRplus) # (dscfAllPairsTest function)
#####

#####

# load in the UAS data into R and format it #######################################
RSEdataframe_all <- read.csv(
  ".\\StordalenMire_ThawPond_UASPolygons_2014to2018.csv", header = TRUE)
# format the column containing Burkeetal2019_ID
RSEdataframe_all$Burkeetal2019_ID <- as.factor(RSEdataframe_all$Burkeetal2019_ID)
# format the column containing FlightDate
RSEdataframe_all$FlightDate <- as.POSIXct(RSEdataframe_all$FlightDate,
                                          format = "%Y-%m-%d", tz = "CET")
# create a column containing only the numerical month of each flight
RSEdataframe_all$month      <- as.numeric(strftime(RSEdataframe_all$FlightDate,
                                                   "%m", tz = "CET"))
# create a column containing only the numerical year of each flight 
RSEdataframe_all$year       <- as.numeric(strftime(RSEdataframe_all$FlightDate,
                                                   "%Y", tz = "CET"))
# create a column containing only the numerical day of year (doy) of each flight
RSEdataframe_all$doy        <- as.numeric(strftime(RSEdataframe_all$FlightDate,
                                                   "%j", tz = "CET"))
# create a unique grouping identifier that combines the Burkeetal2019_ID and the year
RSEdataframe_all$PondYr     <- paste(RSEdataframe_all$Burkeetal2019_ID,
                                     RSEdataframe_all$year,sep = " ")

#create a unique grouping identifier that combines the Month and the Year
RSEdataframe_all$MonthYr    <- paste(RSEdataframe_all$month,
                                     RSEdataframe_all$year, sep = " ")
# create unique grouping identifier that combines the Flight date and Pond and UAS
RSEdataframe_all$FlightDatePondUAS <-paste(RSEdataframe_all$FlightDate,
                                           RSEdataframe_all$Burkeetal2019_ID,
                                           RSEdataframe_all$UAStype, sep = " ")

#######################################################################################


# create appropriate subsets of the RSEdataframe_all for plotting and stats
Quad_PondEdge_PondPolygons  <- subset(RSEdataframe_all, 
                                      RSEdataframe_all$UAStype == "Quad" & 
                                        RSEdataframe_all$PolygonType == "pondedge")
Quad_Water_PondPolygons     <- subset(RSEdataframe_all, 
                                      RSEdataframe_all$UAStype == "Quad" & 
                                        RSEdataframe_all$PolygonType == "water")
FWing_PondEdge_PondPolygons <- subset(RSEdataframe_all, 
                                      RSEdataframe_all$UAStype == "FWing" & 
                                        RSEdataframe_all$PolygonType == "pondedge")
FWing_Water_PondPolygons <- subset(RSEdataframe_all, 
                                   RSEdataframe_all$UAStype == "FWing" & 
                                     RSEdataframe_all$PolygonType == "water")
AllPondEdge_QuadFWing_Polygons <- subset(RSEdataframe_all, 
                                         RSEdataframe_all$PolygonType == "pondedge")
AllWater_QuadFWing_Polygons <- subset(RSEdataframe_all, RSEdataframe_all$PolygonType == "water")


AllWater_QuadFWing_Polygons_PondA <- subset(AllWater_QuadFWing_Polygons,
                                            AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "A")
AllWater_QuadFWing_Polygons_PondB <- subset(AllWater_QuadFWing_Polygons,
                                            AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "B")
AllWater_QuadFWing_Polygons_PondC <- subset(AllWater_QuadFWing_Polygons,
                                            AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "C")
AllWater_QuadFWing_Polygons_PondD <- subset(AllWater_QuadFWing_Polygons,
                                            AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "D")
AllWater_QuadFWing_Polygons_PondE <- subset(AllWater_QuadFWing_Polygons,
                                            AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "E")
AllWater_QuadFWing_Polygons_PondF <- subset(AllWater_QuadFWing_Polygons,
                                            AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "F")
AllWater_QuadFWing_Polygons_PondH <- subset(AllWater_QuadFWing_Polygons,
                                            AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "H")


#subsets for Figure S4,S5,S7
Quad_Water_PondPolygons_PondA <- subset(Quad_Water_PondPolygons,
                                        Quad_Water_PondPolygons$Burkeetal2019_ID == "A")
Quad_Water_PondPolygons_PondB <- subset(Quad_Water_PondPolygons,
                                        Quad_Water_PondPolygons$Burkeetal2019_ID == "B")
Quad_Water_PondPolygons_PondC <- subset(Quad_Water_PondPolygons,
                                        Quad_Water_PondPolygons$Burkeetal2019_ID == "C")
Quad_Water_PondPolygons_PondD <- subset(Quad_Water_PondPolygons,
                                        Quad_Water_PondPolygons$Burkeetal2019_ID == "D")
Quad_Water_PondPolygons_PondE <- subset(Quad_Water_PondPolygons,
                                        Quad_Water_PondPolygons$Burkeetal2019_ID == "E")
Quad_Water_PondPolygons_PondF <- subset(Quad_Water_PondPolygons,
                                        Quad_Water_PondPolygons$Burkeetal2019_ID == "F")
Quad_Water_PondPolygons_PondH <- subset(Quad_Water_PondPolygons,
                                        Quad_Water_PondPolygons$Burkeetal2019_ID == "H")

# subsets for Figure S6, S9
AllPondEdge_QuadFWing_PondA <- subset(AllPondEdge_QuadFWing_Polygons,
                                      AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "A")
AllPondEdge_QuadFWing_PondB <- subset(AllPondEdge_QuadFWing_Polygons,
                                      AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "B")
AllPondEdge_QuadFWing_PondC <- subset(AllPondEdge_QuadFWing_Polygons,
                                      AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "C")
AllPondEdge_QuadFWing_PondD <- subset(AllPondEdge_QuadFWing_Polygons,
                                      AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "D")
AllPondEdge_QuadFWing_PondE <- subset(AllPondEdge_QuadFWing_Polygons,
                                      AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "E")
AllPondEdge_QuadFWing_PondF <- subset(AllPondEdge_QuadFWing_Polygons,
                                      AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "F")
AllPondEdge_QuadFWing_PondH <- subset(AllPondEdge_QuadFWing_Polygons,
                                      AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "H")

#subset the Quad_Water_PondPolygons matrix to only containg data from the month of July
Quad_Water_PondPolygons_JulyOnly <- subset(Quad_Water_PondPolygons,
                                           Quad_Water_PondPolygons$month == 7)

## summarizing the data using ddply################################################################
# summarize the water polygon area data collected with the Quadcopter UAS by month
# calculate the median, maximum, and minimum water polygon area (m2) for
# all ponds across the different sampling seasons.
Summary_QuadWA_month <- ddply(Quad_Water_PondPolygons,.(month),summarize,
                              Median_QWA = median(area_m2,na.rm = TRUE),
                              Max_QWA = max(area_m2,na.rm = TRUE),
                              Min_QWA = min(area_m2,na.rm = TRUE))

# summarize the pond edge polygon area data collected with both the Quadcopter
# and Fixed Wing UAS by year
# calculate the median, maximum, and minimum Pond Edge polygon area (m2)
Summary_AllPEA_year <- ddply(AllPondEdge_QuadFWing_Polygons,.(year),summarize,
                             Median_PEA = median(area_m2,na.rm = TRUE),
                             Max_PEA = max(area_m2,na.rm = TRUE),
                             Min_PEA = min(area_m2,na.rm = TRUE))

# summarize the water polygon area data collected with both the quadcopter and fixed wing UAS
# for each month within each year. 
# calculate the mean and standard deviation polygon area (m2)

Summary_AllWA_monthyear <- ddply(AllWater_QuadFWing_Polygons,.(MonthYr),summarize,
                                 Mean_WA = mean(area_m2,na.rm = TRUE),
                                 STDev_WA = sd(area_m2,na.rm = TRUE))


######################################################################################

#summarize the water polygon area collected with the Quadcopter UAS by pond
# calculate median area (m2)
Summary_QWAbyPond <- ddply(Quad_Water_PondPolygons,.(Burkeetal2019_ID), summarize,
                           MedQWA_pond = median(area_m2,na.rm = TRUE))

# create a subset of the AllPondEdge_QuadFWing_Polygons matrix to include only
# flights that occured in the month of July
AllPondEdge_QuadFWing_Polygons_JulyOnly <- subset(AllPondEdge_QuadFWing_Polygons,
                                                  AllPondEdge_QuadFWing_Polygons$month == 7)        

# In order to test the validity of including the FWing UAS imagery in time series 
# analysis we calculated the average water polygon area and pond edge polygon 
# area from Quadcopter UAS imagery for each pond across each sampling season. 
# These values were then matched to the corresponding PondYr in the matrix 
# containing only FWing UAS imagery. We calculated averages among Quadcopter 
# imagery collected ONLY in July (since the FWing imagery was only collected in 
# July) as well as across the entire field season. 

# calculate the mean of the water polygon area from Quadcopter UAS imagery
# by the PondYr grouping indentifier
Summary_QuadArea_PYr <- ddply(Quad_Water_PondPolygons,.(PondYr), summarize,
                              Avg_WA_Quad = mean(area_m2, na.rm = TRUE))
# calculate the mean of the pond edge polygon area from Quadcopter UAS imagery
# by the PondYr grouping indentifier
Summary_QuadPEArea_Pyr <- ddply(Quad_PondEdge_PondPolygons,.(PondYr), summarize,
                                Avg_PEA_Quad = mean(area_m2, na.rm = TRUE))

# calculate a mean polygon area using the PondYr grouping identifier
Summary_QuadArea_JulyOnly_PYr <- ddply(Quad_Water_PondPolygons_JulyOnly,.(PondYr),
                                       summarize,Avg_JulWA_Quad = mean(area_m2, na.rm = TRUE))

#subset the Quad_PondEdge_PondPolygons matrix to only containg data from the month of July
Quad_PondEdge_PondPolygons_JulyOnly <- subset(Quad_PondEdge_PondPolygons,
                                              Quad_PondEdge_PondPolygons$month == 7)
# calculate a mean polygon area using the PondYr grouping identifier 
Summary_QuadPEArea_JulyOnly_PYr <- ddply(Quad_PondEdge_PondPolygons_JulyOnly,.(PondYr),summarize,
                                         Avg_JulPEA_Quad = mean(area_m2, na.rm = TRUE))


# Match the Avg area from July Only and all Quad to FWing_Water_PondPolygons by PondYr
FWing_Water_PondPolygons$Avg_WA_Quad <- Summary_QuadArea_PYr$Avg_WA_Quad[match(
  FWing_Water_PondPolygons$PondYr,Summary_QuadArea_PYr$PondYr)]

FWing_Water_PondPolygons$Avg_JulWA_Quad <- Summary_QuadArea_JulyOnly_PYr$Avg_JulWA_Quad[
  match(FWing_Water_PondPolygons$PondYr,Summary_QuadArea_JulyOnly_PYr$PondYr)]

# Match the Avg area from the whole field season and all Quad to FWing_Water_PondPolygons by PondYr
FWing_PondEdge_PondPolygons$Avg_PEA_Quad <- Summary_QuadPEArea_Pyr$Avg_PEA_Quad[match(
  FWing_PondEdge_PondPolygons$PondYr,Summary_QuadPEArea_Pyr$PondYr)]

FWing_PondEdge_PondPolygons$Avg_JulPEA_Quad <- Summary_QuadPEArea_JulyOnly_PYr$Avg_JulPEA_Quad[
  match(FWing_PondEdge_PondPolygons$PondYr,Summary_QuadPEArea_JulyOnly_PYr$PondYr)]

##
# calculate the median pond edge area from both Quadcopter and Fixed Wing UAS
# imagery by the grouping identifier PondYr
MedPEA_byPondYr <- ddply(AllPondEdge_QuadFWing_Polygons,.(PondYr), summarize,
                         MedPEA_pondyr = median(area_m2,na.rm = TRUE),
                         MedRatioPEA_pondyr = median(edge.area, na.rm = TRUE))
# match the cumulative bubble flux to the same PondYr in MedPEA_byPondYr
MedPEA_byPondYr$CumFlux <- AllPondEdge_QuadFWing_Polygons$Cumulative_BubFlux[match(MedPEA_byPondYr$PondYr,AllPondEdge_QuadFWing_Polygons$PondYr)]
# separate the group identifier back out into Burkeetal2019_ID and Year
MedPEA_byPondYr <- separate(MedPEA_byPondYr,col = PondYr, into = c("Burkeetal2019_ID","Year"), remove = TRUE)

#summarize the data in MedPEA_byPondYr (which had a separate median area for each field season)
# into a single median area (m2) across all three field seasons. Calculate also the 25% and 75% percentile,
# and the same thing for a median cumulative flux across all three field seasons. 
Summary_MedPEA_MedCumFl_allyrs_allponds <- ddply(MedPEA_byPondYr,.(Burkeetal2019_ID),summarize,
                                                 MedPEA_byPond_allyrs = median(MedPEA_pondyr, na.rm = TRUE),
                                                 Quant25_PEA_allyrs = quantile(MedPEA_pondyr,.25,na.rm = TRUE),
                                                 Quant75_PEA_allyrs = quantile(MedPEA_pondyr,.75,na.rm = TRUE),
                                                 MedRatioPEA_byPond_allyrs = median(MedRatioPEA_pondyr, na.rm = TRUE),
                                                 Quant25_RatioPEA_allyrs = quantile(MedRatioPEA_pondyr,.25,na.rm = TRUE),
                                                 Quant75_RatioPEA_allyrs = quantile(MedRatioPEA_pondyr,.75,na.rm = TRUE),
                                                 MedCumFl_allyrs = median(CumFlux,na.rm = TRUE),
                                                 Quant25_CumFl_allyrs = quantile(CumFlux,.25,na.rm = TRUE),
                                                 Quant75_CumFl_allyrs = quantile(CumFlux,.75,na.rm = TRUE))


#######

## Colors for Plotting ############################################################
# for plotting each pond as a separate boxplot
# this matches the color scheme used in Burke et al. (2019)

# These colors were chosen with the help of the Colorbrewer website:
# https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3 

Pondcolors = c(rgb(254,224,144,maxColorValue = 255),#D (A)
               rgb(224,243,248,maxColorValue = 255),#E (B)
               rgb(69,117,180,maxColorValue = 255), # S (C)
               rgb(253,174,97,maxColorValue = 255),#C (D)
               rgb(215,48,39, maxColorValue = 255),#A (E)
               rgb(244,109,67,maxColorValue = 255),#B (F)
               #rgb(171,217,233, maxColorValue = 255),#F (G)
               rgb(116,173,209,maxColorValue = 255))#R (H)
# plotting each month separately
Monthcolors = c(rgb(237,248,251,maxColorValue = 255), #June
                rgb(179,205,227, maxColorValue = 255), #July
                rgb(140,150,198, maxColorValue = 255)) #August

# color pallets for plotting. Some ponds had no UAS imagery in particular years,
# which is why the pallets are separated here.
Pondcolors_201817 <-c(rgb(254,224,144,maxColorValue = 255),#D (A)
                      rgb(224,243,248,maxColorValue = 255),#E (B)
                      rgb(69,117,180,maxColorValue = 255), # S (C)
                      #rgb(253,174,97,maxColorValue = 255),#C (D)
                      rgb(215,48,39, maxColorValue = 255),#A (E)
                      rgb(244,109,67,maxColorValue = 255),#B (F)
                      #rgb(171,217,233, maxColorValue = 255),#F (G)
                      rgb(116,173,209,maxColorValue = 255))#R (H)
Pondcolors_2016 <-c(rgb(254,224,144,maxColorValue = 255),#D (A)
                    rgb(224,243,248,maxColorValue = 255),#E (B)
                    rgb(69,117,180,maxColorValue = 255), # S (C)
                    rgb(253,174,97,maxColorValue = 255),#C (D)
                    rgb(215,48,39, maxColorValue = 255),#A (E)
                    rgb(244,109,67,maxColorValue = 255),#B (F)
                    #rgb(171,217,233, maxColorValue = 255),#F (G)
                    rgb(116,173,209,maxColorValue = 255))#R (H)

#For plotting each year separately
SamplingSeasons_colors = c(
  rgb(237,248,251, maxColorValue = 255), #2014
  rgb(178,226,226, maxColorValue = 255), #2015
  rgb(102,194,164, maxColorValue = 255), #2016
  rgb(44,162,95, maxColorValue = 255), #2017
  rgb(0,109,44, maxColorValue = 255)) # 2018

######



# create empty data frame to be used to calculate ratio of water area and pond edge area for each fligth date
# there are 163 unique FlightDatePond values
PondsFillUp <- data.frame(matrix(data=NA,nrow = 170,ncol = 1))
names(PondsFillUp) <- c("FlightDatePondUAS")
PondsFillUp$FlightDatePondUAS<- paste(unique(RSEdataframe_all$FlightDatePondUAS))
PondsFillUp <- separate(PondsFillUp,col = FlightDatePondUAS, into = c("FlightDate","Burkeetal2019_ID","UASType"), sep = " ", remove = FALSE)
PondsFillUp$FlightDate <- as.POSIXct(PondsFillUp$FlightDate,
                                     format = "%Y-%m-%d", tz = "CET")
# create a column containing only the numerical month of each flight
PondsFillUp$month      <- as.numeric(strftime(PondsFillUp$FlightDate,
                                              "%m", tz = "CET"))
# create a column containing only the numerical year of each flight 
PondsFillUp$year       <- as.numeric(strftime(PondsFillUp$FlightDate,
                                              "%Y", tz = "CET"))

PondsFillUp$PondEdge_area <- AllPondEdge_QuadFWing_Polygons$area_m2[match(PondsFillUp$FlightDatePondUAS, AllPondEdge_QuadFWing_Polygons$FlightDatePondUAS)]
PondsFillUp$Water_area <- AllWater_QuadFWing_Polygons$area_m2[match(PondsFillUp$FlightDatePondUAS, AllWater_QuadFWing_Polygons$FlightDatePondUAS)]
PondsFillUp$WPERatio_percent <- (PondsFillUp$Water_area/PondsFillUp$PondEdge_area)*100
PondsFillUp$TotPrec_8dpreflight <- RSEdataframe_all$TotPrec_8dpreflight[match(PondsFillUp$FlightDate,RSEdataframe_all$FlightDate)]

PondsFillUp$PondMoYr <- paste(PondsFillUp$Burkeetal2019_ID,PondsFillUp$month, PondsFillUp$year, sep= " ")

PondsFillUp_A <- subset(PondsFillUp,PondsFillUp$Burkeetal2019_ID =="A")
PondsFillUp_B <- subset(PondsFillUp,PondsFillUp$Burkeetal2019_ID =="B")
PondsFillUp_C <- subset(PondsFillUp,PondsFillUp$Burkeetal2019_ID =="C")
PondsFillUp_D <- subset(PondsFillUp,PondsFillUp$Burkeetal2019_ID =="D")
PondsFillUp_E <- subset(PondsFillUp,PondsFillUp$Burkeetal2019_ID =="E")
PondsFillUp_F <- subset(PondsFillUp,PondsFillUp$Burkeetal2019_ID =="F")
PondsFillUp_H <- subset(PondsFillUp,PondsFillUp$Burkeetal2019_ID =="H")


Summary_PondsFillUp <- ddply(PondsFillUp,.(PondMoYr), summarize,
                             Median_WPERatio = median(WPERatio_percent, na.rm = TRUE),
                             Quant25_WPERatio = quantile(WPERatio_percent,0.25, na.rm = TRUE),
                             Quant75_WPERatio = quantile(WPERatio_percent,0.75, na.rm = TRUE))

Summary_PondsFillUp <- separate(Summary_PondsFillUp, col = PondMoYr, into = c("Burkeetal2019_ID","month","year"), sep = " ", remove = FALSE)

Summary_PondsFillUp$MoYr <- paste(Summary_PondsFillUp$month, Summary_PondsFillUp$year, sep = " ")

Summary_PondsFillUp$MoYr_order <- ifelse(Summary_PondsFillUp$MoYr == "7 2014",2,
                                         ifelse(Summary_PondsFillUp$MoYr == "7 2015",6,
                                                ifelse(Summary_PondsFillUp$MoYr == "6 2016",9,
                                                       ifelse(Summary_PondsFillUp$MoYr == "7 2016",10,
                                                              ifelse(Summary_PondsFillUp$MoYr == "8 2016",11,
                                                                     ifelse(Summary_PondsFillUp$MoYr == "6 2017",13,
                                                                            ifelse(Summary_PondsFillUp$MoYr == "7 2017",14,
                                                                                   ifelse(Summary_PondsFillUp$MoYr == "8 2017",15,
                                                                                          ifelse(Summary_PondsFillUp$MoYr == "6 2018",17,
                                                                                                 ifelse(Summary_PondsFillUp$MoYr == "7 2018",18,
                                                                                                        ifelse(Summary_PondsFillUp$MoYr == "8 2018",19,NA)))))))))))




##Code for Plot creation and Statistics#######################################
# the code used for plot creation and related statistical analysis can be found 
# below, any stats presented on particular plots can be found below each plot.

# Figure 3 #####

# Plot Description: This three panel plot shows Water Area (m^2) measured from 
# both Quadcopter and Fixed Wing imagery by pond (panel A), Pond Depression Area measured from 
# both Quadcopter and Fixed Wing imagery by pond (panel B), and the ratio between the two areas (panel C).
# Significant results from the Kruskal-Wallis Ranks Sum test (H, d.f, p) 
# is displayed at the top of each figure figure. Lowercase letters above the boxplots 
# represent significant pairwise differences between ponds and pond types. Numbers
# below the x-axis represent the number of images included below panels A and B.


## save the plot
#pdf("Figure03.pdf",width= 11, height = 8.5)

plot.new()
par(mfrow=c(1,3))

#A.)
# set the dimensions of the top plot
par(mar = c(2.5,4,1,0.5))
par(fig=c(0,1,.68,1), new=TRUE,oma = c(1,1,1,0), cex.axis = 1)

boxplot(area_m2~Burkeetal2019_ID, data = AllWater_QuadFWing_Polygons, col=Pondcolors,
        ylim = c(0,400),ylab = " ", xlab = " ", outpch = 16, axes = FALSE)
axis(1, at = c(1,2,3,4,5,6,7),lab = c("A","B","C","D","E","F","H"), cex.axis = 1.25)
axis(2,at = c(0,100,200,300,400,500,600),las =2,cex.axis=1.25)
mtext(side = 3, expression (bold("A.)")), at  = -0.25, line = 0, cex = 1)
text(7.25, 400, expression("H = 127.56"), cex = 1.25)
text(7.25, 370, expression("d.f. = 6"), cex = 1.25)
text(7.25, 340, expression(italic("p")~"< 0.0001"), cex = 1.25)
text(1, 365, expression("a"), cex =1.25) # pond A
text(2, 110, expression("b"), cex =1.25) # pond B
text(3, 85, expression("c"), cex = 1.25) # pond C
text(4, 195, expression("cd"), cex = 1.25) # pond D
text(5, 205, expression("a"), cex = 1.25) # pond E
text(6, 163, expression("d"), cex = 1.25) # pond F
text(7, 285, expression("e"), cex = 1.25) # pond H
mtext(side = 1, expression(bolditalic("n = ")), line = 2.25, outer = FALSE, 
      at = 0.25, cex = .75 )
mtext(side = 1, expression(bolditalic("25")), line = 2.25, outer = FALSE, at = 1, 
      cex = .75)
mtext(side = 1, expression(bolditalic("24")), line = 2.25, outer = FALSE, at = 2, 
      cex = .75)
mtext(side = 1, expression(bolditalic("28")), line = 2.25, outer = FALSE, at = 3, 
      cex = .75)
mtext(side = 1, expression(bolditalic("20")), line = 2.25, outer = FALSE, at = 4, 
      cex = .75)
mtext(side = 1, expression(bolditalic("24")), line = 2.25, outer = FALSE, at = 5, 
      cex = .75)
mtext(side = 1, expression(bolditalic("24")), line = 2.25, outer = FALSE, at = 6, 
      cex = .75)
mtext(side = 1, expression(bolditalic("25")), line = 2.25, outer = FALSE, at = 7, 
      cex = .75)
mtext(side = 2,expression("Water Area ("*~m^{2}*")"),line = 3,
      outer = FALSE, cex = 1)
box(lty =1)
#
# Use the kruskal function to perform a Kruskal-Wallis test, providing a H
# statistic and a p-value, as well as perform a multiple pairwise comparison
# and provide the final lower case letters signifying significant differences.
kruskal(AllWater_QuadFWing_Polygons$area_m2,AllWater_QuadFWing_Polygons$Burkeetal2019_ID,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# I used the kruskal.test() function to double check what the p value of the test
# is given the kruskal() function cannot approximate very significant p values
kruskal.test(area_m2~Burkeetal2019_ID, data = AllWater_QuadFWing_Polygons)
# Given the above test reports the same H and df values as krusal() I will list
# the p value of this test as <0.0001

# set the dimensions of the middle plot
par(mar = c(2.5,4,1,0.5))
#par(mar = c(5,5,4,1))
#
par(fig=c(0,1,.35,.65),new = TRUE,oma = c(1,1,1,0), cex.axis = 1)
boxplot(area_m2~Burkeetal2019_ID, data = AllPondEdge_QuadFWing_Polygons, 
        col=Pondcolors,ylim = c(0,575),ylab = " ", xlab = " ", outpch = 16, 
        axes = FALSE)
# customize the bottom x axis
axis(1,at = c(1,2,3,4,5,6,7), lab = c("A","B","C","D","E","F","H"), cex.axis = 1.25)
# customize the left y axis
axis(2,at = c(0,100,200,300,400,500,600),las =2, cex.axis= 1.25)
# label the subplot 
mtext(side = 3, expression (bold("B.)")), at  = -0.25, line = 0, cex = 1)
# add the H statistic to the plot
text(7.25, 560, expression("H = 151.68"), cex = 1.25)
# add the df to the plot
text(7.25, 520, expression("d.f. = 6"), cex = 1.25)
# add the p value to the plot
text(7.25, 480, expression(italic("p")~"< 0.0001"), cex = 1.25)
# add the lowercase letters that represent pairwise comparisions to above their
# respective boxplots (see kruskal() call below)
text(1, 225, expression("a"), cex =1.25) # Pond A
text(2, 90, expression("b"), cex =1.25) # Pond B
text(3, 85, expression("c"), cex = 1.25) # Pond C
text(4, 560, expression("d"), cex = 1.25) # Pond D
text(5, 190, expression("e"), cex = 1.25) # Pond E
text(6, 160, expression("e"), cex = 1.25) # Pond F
text(7, 245, expression("f"), cex = 1.25)# Pond H
# add the number of images represented by each boxplot below the x axis label
mtext(side = 1, expression(bolditalic("n =   ")), line = 2.25, outer = FALSE, 
      at = 0.25, cex = .75 )
mtext(side = 1, expression(bolditalic("25")), line = 2.25, outer = FALSE, at = 1, 
      cex = .75)
mtext(side = 1, expression(bolditalic("24")), line = 2.25, outer = FALSE, at = 2, 
      cex = .75)
mtext(side = 1, expression(bolditalic("28")), line = 2.25, outer = FALSE, at = 3, 
      cex = .75)
mtext(side = 1, expression(bolditalic("17")), line = 2.25, outer = FALSE, at = 4, 
      cex = .75)
mtext(side = 1, expression(bolditalic("24")), line = 2.25, outer = FALSE, at = 5, 
      cex = .75)
mtext(side = 1, expression(bolditalic("24")), line = 2.25, outer = FALSE, at = 6, 
      cex = .75)
mtext(side = 1, expression(bolditalic("25")), line = 2.25, outer = FALSE, at = 7, 
      cex = .75)
# label the y axis
mtext(side = 2,expression("Pond Depression Area ("*~m^{2}*")"),line = 3,outer = FALSE,
      cex = 1)
# add a black boarder around the plot
box(lty =1)

kruskal(AllPondEdge_QuadFWing_Polygons$area_m2,AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
kruskal.test(AllPondEdge_QuadFWing_Polygons$area_m2,AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID)
# Given the above test reports the same H and df values as krusal() I will list
# the p value of this test as <0.0001


# set the dimensions of the bottom plot
par(mar = c(3.5,4,1,0.5))
par(fig=c(0,1,0,.32), new=TRUE,oma = c(1,1,1,0), cex.axis = 1)
boxplot(PondsFillUp$WPERatio_percent~PondsFillUp$Burkeetal2019_ID, 
        col = Pondcolors, ylab = " ", xlab = " ", ylim = c(0,178), axes= FALSE)
mtext(side = 3, expression (bold("C.)")), at  = -0.25, line = 0, cex = 1)
mtext(side = 2,expression("Ratio of Water Area" ),line = 3.75,outer = FALSE,
      cex = 1)
mtext(side = 2,expression("to Pond Depression Area (%)" ),line = 2.5,outer = FALSE,
      cex = 1)
axis(2, at = c(0,50,100,150), las = 2, cex.axis = 1.25)
axis(1,at = c(1,2,3,4,5,6,7), lab = c("A","B","C","D","E","F","H"), cex.axis = 1.25)
mtext(side = 1, expression("Pond"), line = 3.5, outer = FALSE, cex = 1)
kruskal(PondsFillUp$WPERatio_percent,PondsFillUp$Burkeetal2019_ID,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
text(1, 178, expression("ab"), cex =1.25) # Pond A
text(2, 120, expression("c"), cex =1.25) # Pond B
text(3, 117, expression("ab"), cex = 1.25) # Pond C
text(4, 60, expression("d"), cex = 1.25) # Pond D
text(5, 130, expression("a"), cex = 1.25) # Pond E
text(6, 135, expression("bc"), cex = 1.25) # Pond F
text(7, 130, expression("a"), cex = 1.25)# Pond H

text(7.25, 175, expression("H = 76.97"), cex = 1.25)
# add the df to the plot
text(7.25, 165, expression("d.f. = 6"), cex = 1.25)
# add the p value to the plot
text(7.25, 150, expression(italic("p")~"< 0.0001"), cex = 1.25)
box(lty =1)
abline(h=100, lty = 2)

# 
# # #####

#Figure 4 ##########################################
#Plot Description: This figure shows the water area (m2) measured from both quadcopter and fixed wing imagery
# during June, July and August in 2016, 2017 and 2018 as the top plot.The months are represented by different colors. The number of 
# individual polygons represented in each boxplot are displayed below the x axis 
# label.Then below the top plot, three plots show Water Area for each month separated out by year. 
## save the plot
#pdf("Figure04.pdf",width= 11, height = 8.5)
par(mfrow =c(1,1))
plot.new()
par(fig=c(0.25,0.75,0.5,0.98), new=TRUE,mar = c(4.5,5,0,0), oma = c(1,1,1,1), 
    cex.axis = 1)

boxplot(area_m2~month, data = AllWater_QuadFWing_Polygons,
        col = Monthcolors, axes = FALSE, ylab = " ", xlab = " ",ylim = c(0,400),
        outpch = 16)
mtext(side = 2,expression("Water Area ("*~m^{2}*")"),line = 4,
      outer = FALSE, cex = 1.5)
#mtext(side = 2, expression(italic("from Quadcopter and Fixed Wing Imagery")), line = 3, outer = FALSE, cex = 1.25)
#mtext(side = 1,expression("Month"), line = 5,outer = FALSE,cex = 1.5)
axis(1, at = c(1,2,3), lab = c("June","July","August"), cex.axis= 1.25)
#mtext(side = 3, expression("A.)"), line = 1, at = 0, outer = FALSE, cex=1.5)
axis(2,at = c(0,100,200,300,400,500,600), las = 2, cex.axis=1.25)
mtext(side = 1, expression(bolditalic("n = ")), line = 3, outer = FALSE,
      at = 0.25, cex = 1.25 )
mtext(side = 1, expression(bolditalic("52")), line = 3, outer = FALSE, at = 1,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("100")), line = 3, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("18")), line = 3, outer = FALSE, at = 3,
      cex = 1.25)
mtext(side = 3, expression("2014 - 2018"), line = 0.5, at = 2, outer = FALSE, cex=1.5)

box(lty = 1)



kruskal(AllWater_QuadFWing_Polygons$area_m2,AllWater_QuadFWing_Polygons$month, 
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)

AllWater_QuadFWing_Polygons_2016 <- subset(AllWater_QuadFWing_Polygons,AllWater_QuadFWing_Polygons$year == '2016')
AllWater_QuadFWing_Polygons_2017 <- subset(AllWater_QuadFWing_Polygons,AllWater_QuadFWing_Polygons$year == '2017')
AllWater_QuadFWing_Polygons_2018 <- subset(AllWater_QuadFWing_Polygons,AllWater_QuadFWing_Polygons$year == '2018')


par(mfrow =c(1,1))
#plot.new()
par(fig=c(0,.35,0,0.5), new=TRUE,mar = c(6,5,2,0), oma = c(1,1,1,1), 
    cex.axis = 1)

#A.)
boxplot(area_m2~month, data = AllWater_QuadFWing_Polygons_2016,
        col = Monthcolors, axes = FALSE, ylab = " ", xlab = " ",ylim = c(0,400),
        outpch = 16)
mtext(side = 2,expression("Water Area ("*~m^{2}*")"),line = 4,
      outer = FALSE, cex = 1.5)
#mtext(side = 2, expression(italic("from Quadcopter and Fixed Wing Imagery")), line = 3, outer = FALSE, cex = 1.25)
axis(1, at = c(1,2,3), lab = c("June","July","August"), cex.axis= 1.25)
mtext(side = 3, expression("2016"), line = 0.5, at = 2, outer = FALSE, cex=1.5)
axis(2,at = c(0,100,200,300,400,500,600), las = 2, cex.axis=1.25)
sum(AllWater_QuadFWing_Polygons_2016$month == "6")
sum(AllWater_QuadFWing_Polygons_2016$month == "7")
sum(AllWater_QuadFWing_Polygons_2016$month == "8")
mtext(side = 1, expression(bolditalic("n = ")), line = 3, outer = FALSE,
      at = 0.25, cex = 1.25 )
mtext(side = 1, expression(bolditalic("32")), line = 3, outer = FALSE, at = 1,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("48")), line = 3, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("7")), line = 3, outer = FALSE, at = 3,
      cex = 1.25)
box(lty = 1)
kruskal(AllWater_QuadFWing_Polygons_2016$area_m2,
        AllWater_QuadFWing_Polygons_2016$month,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# not significant


par(fig=c(0.33,.67,0,0.5), new=TRUE,mar = c(6,5,2,0), oma = c(1,1,1,1), 
    cex.axis = 1)

#B.)
boxplot(area_m2~month, data = AllWater_QuadFWing_Polygons_2017,
        col = Monthcolors, axes = FALSE, ylab = " ", xlab = " ",ylim = c(0,400),
        outpch = 16)
mtext(side = 1,expression("Month"), line = 5, at = 1.5 , outer = FALSE,cex = 1.5)
axis(1, at = c(1,2,3), lab = c("June","July","August"), cex.axis= 1.25)
mtext(side = 3, expression("2017"), line = 0.5, at = 1.5, outer = FALSE, cex=1.5)
axis(2,at = c(0,100,200,300,400,500,600), las = 2, cex.axis=1.25)
sum(AllWater_QuadFWing_Polygons_2017$month == "6")
sum(AllWater_QuadFWing_Polygons_2017$month == "7")
sum(AllWater_QuadFWing_Polygons_2017$month == "8") ## No imagery collected in August 2017
mtext(side = 1, expression(bolditalic("6")), line = 3, outer = FALSE, at = 1,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("16")), line = 3, outer = FALSE, at = 2,
      cex = 1.25)
box(lty = 1)

kruskal(AllWater_QuadFWing_Polygons_2017$area_m2,
        AllWater_QuadFWing_Polygons_2017$month,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# not significant

par(fig=c(0.65,1,0,0.5), new=TRUE,mar = c(6,5,2,0), oma = c(1,1,1,1), 
    cex.axis = 1)

#C.)
boxplot(area_m2~month, data = AllWater_QuadFWing_Polygons_2018,
        col = Monthcolors, axes = FALSE, ylab = " ", xlab = " ",ylim = c(0,400),
        outpch = 16)
axis(1, at = c(1,2,3), lab = c("June","July","August"), cex.axis= 1.25)
mtext(side = 3, expression("2018"), line = 1, at = 2, outer = FALSE, cex=1.5)
axis(2,at = c(0,100,200,300,400,500,600), las = 2, cex.axis=1.25)
sum(AllWater_QuadFWing_Polygons_2018$month == "6")
sum(AllWater_QuadFWing_Polygons_2018$month == "7")
sum(AllWater_QuadFWing_Polygons_2018$month == "8")

mtext(side = 1, expression(bolditalic("14")), line = 3, outer = FALSE, at = 1,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("24")), line = 3, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("11")), line = 3, outer = FALSE, at = 3,
      cex = 1.25)
box(lty = 1)
kruskal(AllWater_QuadFWing_Polygons_2018$area_m2,
        AllWater_QuadFWing_Polygons_2018$month,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# not significant

# Figure 5########################################################################
# Plot Description: This plot shows Pond Depression Polygon Area (m^2)
# by sampling season with a all ponds plotted together measured from both Quadcopter and Fixed Wing imagery. 
# The sampling seasons are represented by different colors. The number of 
# individual polygons represented in each boxplot are displayed below year x axis 
# label.
#
## save the plot
#pdf("Figure05.pdf",width= 11, height = 8.5)

par(mfrow = c(1,1))
plot.new()
par(mar = c(4,4.75,1,0))
par(fig=c(0,0.48,0,1), new=TRUE,oma = c(1,1,1,0), cex.axis = 1)

#A.)
boxplot(area_m2~year, data = AllPondEdge_QuadFWing_Polygons,
        col = SamplingSeasons_colors, ylim = c(0,600), axes = FALSE,
        ylab = " ",xlab = " ",outpch = 16)
mtext(side = 1, expression ("Year"), line = 4, cex = 1.5, at = 6)
mtext(side = 2, expression ("Pond Depression Area ("*~m^{2}*")"), line = 4,
      cex = 1.5)
#mtext(side = 2, expression(italic("from Quadcopter and Fixed Wing Imagery")), line = 3 ,
#                          cex = 1.25)
mtext(side = 3, expression("A.)"), cex = 1.5, line=0, at = -1)

axis(1, at = c(1,2,3,4,5), lab = c("2014","2015","2016","2017","2018"),
     cex.axis = 1.25)
axis(2,at = c(0,100,200,300,400,500,600), las = 2, cex.axis=1.25)

mtext(side = 1, expression(bolditalic("n = ")), line = 2.5, outer = FALSE,
      at = 0.25, cex = 1.25)
mtext(side = 1, expression(bolditalic("6")), line = 2.5, outer = FALSE, at = 1,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("6")), line = 2.5, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("84")), line = 2.5, outer = FALSE, at = 3,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("22")), line = 2.5, outer = FALSE, at = 4,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("49")), line = 2.5, outer = FALSE, at = 5,
      cex = 1.25)
box(lty = 1, lwd =1)
# #
kruskal(AllPondEdge_QuadFWing_Polygons$area_m2,
        AllPondEdge_QuadFWing_Polygons$year,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# not significant


par(mar = c(4,4,1,0.5))
par(fig=c(0.52,1,0,1), new=TRUE,oma = c(1,1,1,0), cex.axis = 1)
#A.)
boxplot(area_m2~year, data = AllWater_QuadFWing_Polygons,
        col = SamplingSeasons_colors, ylim = c(0,600), axes = FALSE,
        ylab = " ",xlab = " ",outpch = 16)
#mtext(side = 1, expression ("Year"), line = 4, cex = 1.5)
mtext(side = 2, expression ("Water Area ("*~m^{2}*")"), line = 4,
      cex = 1.5)
#mtext(side = 2, expression(italic("from Quadcopter and Fixed Wing Imagery")), line = 3 ,
#     cex = 1.25)
mtext(side = 3, expression("B.)"), cex = 1.5, line=0, at = -1)

axis(1, at = c(1,2,3,4,5), lab = c("2014","2015","2016","2017","2018"),
     cex.axis = 1.25)
axis(2,at = c(0,100,200,300,400,500,600), las = 2, cex.axis=1.25)

mtext(side = 1, expression(bolditalic("n = ")), line = 2.5, outer = FALSE,
      at = 0.25, cex = 1.25)
mtext(side = 1, expression(bolditalic("6")), line = 2.5, outer = FALSE, at = 1,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("6")), line = 2.5, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("87")), line = 2.5, outer = FALSE, at = 3,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("22")), line = 2.5, outer = FALSE, at = 4,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("49")), line = 2.5, outer = FALSE, at = 5,
      cex = 1.25)
box(lty = 1, lwd =1)
# #
kruskal(AllWater_QuadFWing_Polygons$area_m2,
        AllWater_QuadFWing_Polygons$year,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# not significant


########################################################################################
# NOTE: Code to produce Supplementary Figure S1 can be found in the 
# Burkeetal_RSEManuscript_SupplementalFluxCalculations.R Script.

# Supplementary Figure S3#######################################################
# Plot Description: Two panel plot with the left panel (A.) plotting the 
# Average Pond Depression Area (m^2) measured from Quadcopter imagery (y axis) 
# versus Pond Depression Area (m^2) measured from Fixed Wing Imagery (x axis).
# There is a inset plot in panel A that is the same as the larger Panel A plot 
# except that the y axis is an average of Pond Edge Area from Quadcopter imagery
# collected only in July. Panel B shows average Water area of Quadcopter Imagery 
# (y axis) vs. Water Area (m^2) measured from the Fixed Wing Imagery with a inset 
# plot in panel B that is the same as the larger Panel B plot except that the 
# y axis is an average of Pond Edge Area from Quadcopter imagery collected only 
# in July. Each pond is represented by different colors and each sampling season 
#(2016-2018) are represented by different shapes. The y axis is an average while
# the x axis is not because the Quadcopter was flown multiple times each year, 
# while the Fixed Wing drone was flown only once per year. At the top of 
# both plots is the tau value and p value from a Kendall Test. A 1:1 line is also
# plotted on both plots (and inset plots) with the 1:1 line in the larger plots
# slightly offset in order to accommodate the inset plots.
#
## save the plot
#pdf("FigureS03.pdf",width= 11, height = 8.5)

# A.)
plot.new()
par(mar=c(8,8,8,8))

par(fig = c(0.02,.52,0,1), new = TRUE,mar = c(7,4,1,1),oma = c(1,1,1,1), cex.axis=1)
plot(Avg_PEA_Quad~area_m2, data = FWing_PondEdge_PondPolygons, 
     subset = FWing_PondEdge_PondPolygons$year == 2018,xlim = c(0,325), 
     ylim =c(0,275), 
     pch = 21,col = "black", bg= Pondcolors_201817, axes = FALSE,cex =2, 
     ylab =" ",xlab = " ")
abline(0,1,lty = 1)
text(255,273, expression("1:1"), cex = 1.25)
mtext(side = 3, expression("A.)"), cex = 1.5, line=0, at = -40)
mtext(side = 2,expression("Average Pond Depression Area " (m^2)),
      cex = 1.5, line = 3.5)
mtext(side = 2,expression(italic("from Quadcopter Imagery")),
      cex = 1.25, line = 2.5)
mtext(side = 1,expression("Pond Depression Area " (m^2)),line = 3, 
      cex = 1.5)
mtext(side = 1, expression(italic("from Fixed Wing Imagery")), line = 4, cex = 1.25)

axis(1,at = c(50,100,150,200,250), cex.axis= 1.25)
axis(2,at = c(50,100,150,200,250), las= 2,cex.axis= 1.25)
points(Avg_PEA_Quad~area_m2, data = FWing_PondEdge_PondPolygons, 
       subset = FWing_PondEdge_PondPolygons$year == 2017,cex =2,xlim = c(0,300),
       ylim = c(0,300),pch = 22, col = "black", bg = Pondcolors_201817)
points(Avg_PEA_Quad~area_m2, data = FWing_PondEdge_PondPolygons, 
       subset = FWing_PondEdge_PondPolygons$year == 2016,cex=2, xlim = c(0,300), 
       ylim = c(0,300), pch = 23, col = "black", bg = Pondcolors_2016)
box(lty =1 )
text(90,275, expression(tau~"= 0.90"),cex = 1.25)
text(180,275, expression(italic("p")~"< 0.0001"), cex = 1.25)

# A Kendall Correlation is performed to look at significat relationships among
# continuous variables.
cor.test(FWing_PondEdge_PondPolygons$area_m2,
         FWing_PondEdge_PondPolygons$Avg_PEA_Quad,method = "kendall")
# strong correlation
#
# inset plot for Panel A

par(fig = c(0.25,0.49,0.06,0.6), new = TRUE,mar = c(7,4,1,0),oma = c(1,1,1,1), cex.axis=1)
plot(Avg_JulPEA_Quad~area_m2, data = FWing_PondEdge_PondPolygons, 
     subset = FWing_PondEdge_PondPolygons$year == 2018,xlim = c(0,275),
     ylim =c(0,275), 
     pch = 21,col = "black", bg= Pondcolors_201817, axes = FALSE, ylab =" ",
     xlab = " ", cex = 2)
abline(0,1,lty = 1)
text(245,273, expression("1:1"), cex = 1.25)

mtext(side = 3,expression("July Only"), cex = 1.25, line = 0)
axis(1,at = c(50,150,250),cex.axis= 1.25)
axis(2,at = c(50,100,150,200,250), las= 2,cex.axis= 1.25)
points(Avg_JulPEA_Quad~area_m2, data = FWing_PondEdge_PondPolygons,
       subset = FWing_PondEdge_PondPolygons$year == 2017,cex = 2, 
       xlim = c(0,300),ylim = c(0,300),pch = 22, col = "black", 
       bg = Pondcolors_201817)
points(Avg_JulPEA_Quad~area_m2, data = FWing_PondEdge_PondPolygons, 
       subset = FWing_PondEdge_PondPolygons$year == 2016,cex = 2,
       xlim = c(0,300), ylim = c(0,300), pch = 23, col = "black", 
       bg = Pondcolors_2016)
box(lty =1 )
text(90,260, expression(tau~"= 0.91"),cex = 1.25)
text(90,230, expression(italic("p")~"< 0.0001"), cex = 1.25)
# strong correlation
cor.test(FWing_PondEdge_PondPolygons$area_m2,
         FWing_PondEdge_PondPolygons$Avg_JulPEA_Quad,method = "kendall")
#
# create a custom legend for the plot
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(-0.7,-0.95, c("A","B","C","D","E","F","H"), xpd = TRUE, horiz = TRUE, 
       inset = c(0,0), bty = "n", pch = c(NA,NA,NA,NA,NA,NA,NA), col = "black",
       fill = Pondcolors, cex = 1.25)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(.3,-0.95, c("2016","2017","2018"), xpd = TRUE, horiz = TRUE, 
       inset = c(0,0), bty = "n", pch = c(23,22,21), col = "black",cex = 1.25)
#

#B.)
par(fig = c(0.52,1,0,1), new = TRUE,mar = c(7,5,1,0),oma = c(1,1,1,1), cex.axis=1)
plot(Avg_WA_Quad~area_m2, data = FWing_Water_PondPolygons, 
     subset = FWing_Water_PondPolygons$year == 2018,xlim = c(0,325), 
     ylim =c(0,275), 
     pch = 21,col = "black", bg= Pondcolors_201817, axes = FALSE,cex =2, 
     ylab =" ",xlab = " ")
abline(0,1,lty = 1)
text(255,273, expression("1:1"), cex = 1.25)
mtext(side = 3, expression("B.)"), cex = 1.5, line=0, at = -40)

mtext(side = 2,expression("Average Water Area " (m^2)),
      cex = 1.5, line = 3.5)
mtext(side = 2, expression(italic("from Quadcopter Imagery")), cex= 1.25, 
      line = 2.5)

mtext(side = 1,expression("Water Area " (m^2)),line = 3,
      cex = 1.5)
mtext(side = 1, expression(italic("from Fixed Wing Imagery")), line = 4, cex =1.25)
axis(1,at = c(50,100,150,200,250),cex.axis= 1.25)
axis(2,at = c(50,100,150,200,250), las= 2,cex.axis= 1.25)
points(Avg_WA_Quad~area_m2, data = FWing_Water_PondPolygons, 
       subset = FWing_Water_PondPolygons$year == 2017,cex =2,xlim = c(0,300),
       ylim = c(0,300),pch = 22, col = "black", bg = Pondcolors_201817)
points(Avg_WA_Quad~area_m2, data = FWing_Water_PondPolygons, 
       subset = FWing_Water_PondPolygons$year == 2016,cex=2,xlim = c(0,300),
       ylim = c(0,300), pch = 23, col = "black", bg = Pondcolors_2016)
box(lty =1 )
text(90,275, expression(tau~"= 0.78"),cex = 1.25)
text(180,275, expression(italic("p")~"< 0.0001"), cex = 1.25)

# strong correlation
cor.test(FWing_Water_PondPolygons$area_m2,
         FWing_Water_PondPolygons$Avg_WA_Quad,method = "kendall")
#
# inset plot in Panel B. 

par(fig = c(0.73,1,0.11,0.58), new = TRUE,mar = c(5,5,.25,1),oma = c(1,1,1,1), 
    cex.axis=1)
plot(Avg_JulWA_Quad~area_m2, data = FWing_Water_PondPolygons, 
     subset = FWing_Water_PondPolygons$year == 2018,xlim = c(0,275), 
     ylim =c(0,275),pch = 21,col = "black", bg= Pondcolors_201817, 
     axes = FALSE, ylab =" ",xlab = " ", cex = 2)
abline(0,1,lty = 1)
text(245,273, expression("1:1"), cex = 1.25)

mtext(side = 3, expression("July Only"), cex = 1.25)

axis(1,at = c(50,150,250),cex.axis= 1.25)
axis(2,at = c(50,100,150,200,250), las= 2,cex.axis= 1.25)
points(Avg_JulWA_Quad~area_m2, data = FWing_Water_PondPolygons,
       subset = FWing_Water_PondPolygons$year == 2017,cex = 2, 
       xlim = c(0,300),ylim = c(0,300),pch = 22, col = "black", 
       bg = Pondcolors_201817)
points(Avg_JulWA_Quad~area_m2, data = FWing_Water_PondPolygons,
       subset = FWing_Water_PondPolygons$year == 2016,cex = 2,
       xlim = c(0,300), ylim = c(0,300), pch = 23, col = "black", 
       bg = Pondcolors_2016)
box(lty =1 )
text(90,260, expression(tau~"= 0.84"),cex = 1.25)
text(90,230, expression(italic("p")~"< 0.0001"), cex = 1.25)
# strong correlation
cor.test(FWing_Water_PondPolygons$area_m2,
         FWing_Water_PondPolygons$Avg_WA_Quad,method = "kendall")

#dev.off()


######
# Supplementary Figure S4#####################################################################        
# Plot Description: This figure is an expansion on Figure S4. Each pond is 
# plotted separately with boxplots of the area of Water Polygons (m^2) measured 
# from Quadcopter imagery AND Fixed Wing Images broken down by month. Ponds with significant differences 
# in Water Polygon Area by month have a * by their names. The H statistic, d.f. 
# and p values are also plotted (resulting from a Kruskal Wallis Test).
## save the plot
#pdf("FigureS04.pdf",width= 11, height = 8.5)
#
par(mfrow =c(1,1))
plot.new()
# Pond A
par(fig=c(.125,.375,.55,1), new=TRUE,mar = c(3,4,2,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~month, data = AllWater_QuadFWing_Polygons, 
        subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "A",
        col = Monthcolors, ylab = " ", xlab = " ", ylim = c(0,420), 
        axes = FALSE, outpch = 16)
mtext(side = 3, expression ("Pond A*"), cex = 1.5)
axis(1,at =c(1,2,3), lab = c("June","July","August"),cex.axis= 1.25)
axis(2,at = c(0,100,200,300,400), las = 2,cex.axis= 1.25)
mtext(side = 1, expression(bolditalic("n = ")), line = 2.25, outer = FALSE, 
      at = 0.25, cex = 1.25 )
mtext(side = 1, expression(bolditalic("7")), line = 2.25, outer = FALSE, at = 1,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("15")), line = 2.25, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("3")), line = 2.25, outer = FALSE, at = 3,
      cex = 1.25)
text(1,410, expression("H = 6.49"), cex = 1.25)
text(2,410,expression("2 d.f."), cex= 1.25)
text(3,410, expression(italic("p")~"= 0.04"), cex = 1.26)
text(1,370, expression("a"), cex = 1.25)
text(2,310, expression("ab"), cex = 1.25)
text(3,330, expression("ab"), cex = 1.25)
box(lty = 1)
#
kruskal(AllWater_QuadFWing_Polygons_PondA$area_m2,
        AllWater_QuadFWing_Polygons_PondA$month,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# significant
#Pond B
par(fig=c(.375,.625,.55,1), new=TRUE,mar = c(3,4,2,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~month, data = AllWater_QuadFWing_Polygons, 
        subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "B", 
        col = Monthcolors, ylab = " ", xlab = " ", ylim = c(0,100), 
        axes = FALSE, outpch = 16)
mtext(side = 3, expression ("Pond B"), cex = 1.5)
mtext(side = 1, expression(bolditalic("7")), line = 2.25, outer = FALSE, at = 1, 
      cex = 1.25)
mtext(side = 1, expression(bolditalic("15")), line = 2.25, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("2")), line = 2.25, outer = FALSE, at = 3,
      cex = 1.25)
axis(1,at =c(1,2,3), lab = c("June","July","August"),cex.axis= 1.25)
axis(2,at = c(0,25,50,75,100),las = 2,cex.axis= 1.25)
box(lty = 1, lwd =1)
# 
kruskal(AllWater_QuadFWing_Polygons_PondB$area_m2,
        AllWater_QuadFWing_Polygons_PondB$month,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
#not significant
# Pond C 
par(fig=c(.625,.875,.55,1), new=TRUE,mar = c(3,4,2,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~month, data = AllWater_QuadFWing_Polygons, 
        subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "C", 
        col = Monthcolors, ylab = " ", xlab = " ", ylim = c(0,100),
        axes = FALSE, outpch = 16)
mtext(side = 3, expression ("Pond C*"), cex = 1.5)
mtext(side = 1, expression(bolditalic("9")), line = 2.25, outer = FALSE, at = 1, 
      cex = 1.25)
mtext(side = 1, expression(bolditalic("16")), line = 2.25, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("3")), line = 2.25, outer = FALSE, at = 3, 
      cex = 1.25)
axis(1,at =c(1,2,3), lab = c("June","July","August"),cex.axis= 1.25)
axis(2,at = c(0,25,50,75,100),las = 2,cex.axis= 1.25)
text(1,98, expression("H = 8.17"), cex = 1.25)
text(2,98,expression("2 d.f."), cex= 1.25)
text(3,98, expression(italic("p")~"= 0.02"), cex = 1.26)
text(1,70, expression("a"), cex = 1.25)
text(2,65, expression("b"), cex = 1.25)
text(3,60, expression("b"), cex = 1.25)
box(lty = 1, lwd = 1)
#
kruskal(AllWater_QuadFWing_Polygons_PondC$area_m2,
        AllWater_QuadFWing_Polygons_PondC$month,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# significant
# # 
# Pond D
par(fig=c(0.02,.27,0,.5), new=TRUE,mar = c(4,4,2,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~month, data = AllWater_QuadFWing_Polygons,
        subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "D", 
        col = Monthcolors, ylab = " ", xlab = " ", ylim = c(0,220), 
        axes = FALSE, outpch = 16)
mtext(side = 3, expression ("Pond D"), cex = 1.5)
mtext(side = 1, expression(bolditalic("n = ")), line = 2.25, outer = FALSE, 
      at = 0.25, cex = 1.25 )
mtext(side = 2, expression ("Water Area ("*~m^{2}*")"), line = 4, 
      at = 250 , cex = 1.5)
#mtext(side = 2, expression (italic("from Quadcopter & Fixed Wing Imagery")), line = 3, 
#      at = 250 , cex = 1.25)

mtext(side = 1, expression(bolditalic("7")), line = 2.25, outer = FALSE, at = 1, 
      cex = 1.25)
mtext(side = 1, expression(bolditalic("11")), line = 2.25, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("2")), line = 2.25, outer = FALSE, at = 3, 
      cex = 1.25)
axis(1,at =c(1,2,3), lab = c("June","July","August"),cex.axis= 1.25)
axis(2,at = c(0,50,100,150,200),las =2,cex.axis= 1.25)
box(lty = 1, lwd=1)
# 
kruskal(AllWater_QuadFWing_Polygons_PondD$area_m2,
        AllWater_QuadFWing_Polygons_PondD$month,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# not significant
# 
#Pond E
par(fig=c(0.26,.51,0,.5), new=TRUE,mar = c(4,4,2,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~month, data = AllWater_QuadFWing_Polygons, 
        subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "E", 
        col = Monthcolors, ylab = " ", xlab = " ", ylim = c(0,220), 
        axes = FALSE)
mtext(side = 3, expression ("Pond E"), cex = 1.5)
mtext(side = 1, expression(bolditalic("n = ")), line = 2.25, outer = FALSE, 
      at = 0.25, cex = 1.25 )
mtext(side = 1, expression(bolditalic("7")), line = 2.25, outer = FALSE, at = 1,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("14")), line = 2.25, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("3")), line = 2.25, outer = FALSE, at = 3,
      cex = 1.25)
axis(1,at =c(1,2,3), lab = c("June","July","August"),cex.axis= 1.25)
axis(2,at = c(0,50,100,150,200),las = 2,cex.axis= 1.25)
box(lty=1,lwd =1)
# 
kruskal(AllWater_QuadFWing_Polygons_PondE$area_m2,
        AllWater_QuadFWing_Polygons_PondE$month,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# not significant
# Pond F
par(fig=c(.5,.75,0,.5), new=TRUE,mar = c(4,4,2,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~month, data = AllWater_QuadFWing_Polygons, 
        subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "F", 
        col = Monthcolors, ylab = " ", xlab = " ", ylim = c(0,220), 
        axes = FALSE)
mtext(side = 3, expression ("Pond F"), cex = 1.5)

mtext(side = 1, expression(bolditalic("8")), line = 2.25, outer = FALSE, at = 1, 
      cex = 1.25)
mtext(side = 1, expression(bolditalic("13")), line = 2.25, outer = FALSE, at = 2, 
      cex = 1.25)
mtext(side = 1, expression(bolditalic("3")), line = 2.25, outer = FALSE, at = 3, 
      cex = 1.25)
axis(1,at =c(1,2,3), lab = c("June","July","August"),cex.axis= 1.25)
axis(2,at = c(0,50,100,150,200), las = 2,cex.axis= 1.25)
mtext(side = 1, expression ("Month"), line = 4, at =-.25 , cex = 1.5)
box(lty =1, lwd =1)
# 
kruskal(AllWater_QuadFWing_Polygons_PondF$area_m2,
        AllWater_QuadFWing_Polygons_PondF$month,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# not significant
#Pond H
par(fig=c(.75,1,0,.5), new=TRUE, mar = c(4,4,2,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~month, data = AllWater_QuadFWing_Polygons, 
        subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "H",
        col = Monthcolors, ylab = " ", xlab = " ", ylim = c(0,420),
        outpch = 16, axes = FALSE)
axis(1,at =c(1,2,3), lab = c("June","July","August"),cex.axis= 1.25)
mtext(side = 1, expression(bolditalic("7")), line = 2.25, outer = FALSE, at = 1,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("16")), line = 2.25, outer = FALSE, at = 2, 
      cex = 1.25)
mtext(side = 1, expression(bolditalic("2")), line = 2.25, outer = FALSE, at = 3, 
      cex = 1.25)
axis(2,at = c(0,100,200,300,400), las = 2,cex.axis= 1.25)
mtext(side = 3, expression ("Pond H"), cex = 1.5)
box(lty = 1, lwd =1)
# 
kruskal(AllWater_QuadFWing_Polygons_PondH$area_m2,
        AllWater_QuadFWing_Polygons_PondH$month,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)

#dev.off()
# not significant
#####


# 
##########
# Supplementary Figure S5#####################################################################
# Plot Description: This plot shows the relationship between total precipitation
# measured in the eight days leading up to and including the flight date and
# the ratio of water polygon area to pond depression area measured from both Quadcopter and 
# fixed wing Imagery
## save the plot
#pdf("FigureS05.pdf",width= 11, height = 8.5)
par(mfrow = c(1,1))
plot.new()
par(mar = c(5,5,5,5))
par(fig=c(0.02,.27,.5,1), new=TRUE,mar = c(3,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
# All Ponds
plot(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp, subset = PondsFillUp$Burkeetal2019_ID == "A",pch = 21, 
     col = "black", bg = Pondcolors[1], ylim = c(0,175),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
points(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp, subset = PondsFillUp$Burkeetal2019_ID == "B",pch = 21, 
       col = "black", bg = Pondcolors[2], ylim = c(0,175),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
points(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp, subset = PondsFillUp$Burkeetal2019_ID == "C",pch = 21, 
       col = "black", bg = Pondcolors[3], ylim = c(0,175),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
points(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp, subset = PondsFillUp$Burkeetal2019_ID == "D",pch = 21, 
       col = "black", bg = Pondcolors[4], ylim = c(0,175),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
points(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp, subset = PondsFillUp$Burkeetal2019_ID == "E",pch = 21, 
       col = "black", bg = Pondcolors[5], ylim = c(0,175),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
points(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp, subset = PondsFillUp$Burkeetal2019_ID == "F",pch = 21, 
       col = "black", bg = Pondcolors[6], ylim = c(0,175),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
points(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp, subset = PondsFillUp$Burkeetal2019_ID == "H",pch = 21, 
       col = "black", bg = Pondcolors[7], ylim = c(0,175),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
axis(1, at = c(0,20,40,60,80,100), cex.axis = 1.25)
axis(2, at = c(0,25,50,75,100,125,150,175), las = 2, cex.axis = 1.25)
mtext(side = 3, expression("All Ponds*"),cex = 1.5)
mtext(side = 2, expression ("Ratio of Water Area to Pond Depression Area (%)"), line = 3, at = -20, 
      cex = 1.25)
text(80,170 , expression(tau~" = 0.11"), cex = 1)
text(80,160 , expression(italic("p")~"= 0.04"), cex = 1)
abline(h=100, lty= 1)
box(lty =1, lwd = 1)
# the cor.test provides a tau value and p value
cor.test(PondsFillUp$WPERatio_percent,
         PondsFillUp$TotPrec_8dpreflight,method = "kendall")
# not significant

par(fig=c(0.26,.51,.5,1), new=TRUE,mar = c(3,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
# Plot A
plot(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp, 
     subset = PondsFillUp$Burkeetal2019_ID == "A", pch = 21,
     col = "black", bg=rgb(254,224,144, maxColorValue = 255), ylab = " ", 
     xlab = " ",
     axes = FALSE, xlim = c(0,100), ylim = c(0,175))
axis(1,at =c(0,20,40,60,80,100),cex.axis= 1.25)
axis(2, at = c(0,25,50,75,100,125,150,175), las = 2, cex.axis = 1.25)
mtext(side = 3, expression("Pond A*"), cex = 1.5)
text(80,170 , expression(tau~" = 0.27"), cex = 1)
text(80,160 , expression(italic("p")~"= 0.05"), cex = 1)
box(lty = 1,lwd =1)
abline(h=100, lty= 1)
cor.test(PondsFillUp_A$WPERatio_percent,
         PondsFillUp_A$TotPrec_8dpreflight,method = "kendall")
#   significant.
#
# Pond B
par(fig=c(0.5,.75,.5,1), new=TRUE,mar = c(3,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp, 
     subset = PondsFillUp$Burkeetal2019_ID == "B", pch = 21,
     col = "black", bg=rgb(224,243,248, maxColorValue = 255), ylab = " ", 
     xlab = " ", axes = FALSE,
     ylim = c(0,175), xlim = c(0,100))
mtext(side = 3, expression("Pond B*"), cex = 1.5)
axis(1,at =c(0,20,40,60,80,100),cex.axis= 1.25)
axis(2, at = c(0,25,50,75,100,125,150,175), las = 2, cex.axis = 1.25)
text(80,170 , expression(tau~" = 0.44"), cex = 1)
text(80,160 , expression(italic("p")~"= 0.002"), cex = 1)
box(lty = 1,lwd =1)
abline(h=100, lty= 1)
# strong correlation
cor.test(PondsFillUp_B$WPERatio_percent,
         PondsFillUp_B$TotPrec_8dpreflight,method = "kendall")
#   significant
#
# Pond C
par(fig=c(.75,1,.5,1), new=TRUE,mar = c(3,3,0.5,1), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp,
     subset = PondsFillUp$Burkeetal2019_ID == "C", pch = 21,
     col = "black", bg=rgb(69,117,180, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, xlim = c(0,100), ylim = c(0,175))
axis(1,at =c(0,20,40,60,80,100),cex.axis= 1.25)
axis(2, at = c(0,25,50,75,100,125,150,175), las = 2, cex.axis = 1.25)
mtext(side = 3, expression("Pond C"), cex = 1.5)
abline(h=100, lty= 1)
box(lty = 1,lwd =1)
#
#cor.test(PondsFillUp_C$WPERatio_percent,
#         PondsFillUp_C$TotPrec_8dpreflight,method = "kendall")
# not significant
#
#Pond D
par(fig=c(0.02,.27,0,.5), new=TRUE,mar = c(3,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp,
     subset = PondsFillUp$Burkeetal2019_ID == "D", pch = 21,
     col = "black", bg=rgb(253,174,97, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, xlim = c(0,100), ylim = c(0,175))
axis(1,at =c(0,20,40,60,80,100),cex.axis= 1.25)
axis(2, at = c(0,25,50,75,100,125,150,175), las = 2, cex.axis = 1.25)
mtext(side = 3, expression("Pond D"), cex = 1.5)
abline(h=100, lty= 1)
box(lty = 1,lwd =1)
#
#cor.test(PondsFillUp_D$WPERatio_percent,
#         PondsFillUp_D$TotPrec_8dpreflight,method = "kendall")
# not significant
#
#PondE 
par(fig=c(0.26,.51,0,.5), new=TRUE,mar = c(3,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp,
     subset = PondsFillUp$Burkeetal2019_ID == "E", pch = 21,
     col = "black", bg=rgb(215,48,39, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, ylim = c(0,175), xlim= c(0,100))
axis(1,at =c(0,20,40,60,80,100),cex.axis= 1.25)
axis(2, at = c(0,25,50,75,100,125,150,175), las = 2, cex.axis = 1.25)
abline(h=100, lty= 1)
mtext(side = 3, expression("Pond E"), cex =1.5)
box(lty = 1,lwd =1)
#
#cor.test(PondsFillUp_E$WPERatio_percent,
#        PondsFillUp_E$TotPrec_8dpreflight,method = "kendall")
# not significant 
#
# Pond F
par(fig=c(.5,.75,0,.5), new=TRUE,mar = c(3,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp,
     subset = PondsFillUp$Burkeetal2019_ID == "F", pch = 21,
     col = "black", bg=rgb(244,109,67, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, ylim = c(0,175), xlim =c(0,100))
axis(1,at =c(0,20,40,60,80,100),cex.axis= 1.25)
axis(2, at = c(0,25,50,75,100,125,150,175), las = 2, cex.axis = 1.25)
mtext(side = 3, expression("Pond F"), cex = 1.5)
box(lty = 1,lwd =1)
abline(h=100, lty= 1)
mtext(side = 1, expression ("Total Precipitation (mm) in the 8 days before Flight"), line = 3,at = -10 , cex = 1.25)
#
#cor.test(PondsFillUp_F$WPERatio_percent,
#        PondsFillUp_F$TotPrec_8dpreflight,method = "kendall")
# not significant
#
#PondH
par(fig=c(0.75,1,0,.5), new=TRUE,mar = c(3,3,0.5,1), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(WPERatio_percent~TotPrec_8dpreflight, data = PondsFillUp,
     subset = PondsFillUp$Burkeetal2019_ID == "H",pch = 21,
     col = "black", bg=rgb(116,173,209, maxColorValue = 255), ylab = " ", 
     xlab = " ",
     axes = FALSE, ylim = c(0,175), xlim = c(0,100))
axis(1,at =c(0,20,40,60,80,100),cex.axis= 1.25)
axis(2, at = c(0,25,50,75,100,125,150,175), las = 2, cex.axis = 1.25)
mtext(side = 3, expression("Pond H"), cex = 1.5)
abline(h=100, lty= 1)
box(lty = 1,lwd =1)
#
#cor.test(PondsFillUp_H$WPERatio_percent,
#        PondsFillUp_H$TotPrec_8dpreflight,method = "kendall")
# not significant

#Figure S6################################
# Plot Description: This plot shows the median ratio of water area to Pond Depression Area (%)
# by month from 2014 to 2018. error bars are displayed on data in 2016-2018 that has multiple images available
# (this varies by pond). Each pond is displayed as a different color and/or symbol. 
# Total precipitation (mm) measured at Stordalen Mire is displayed below the x axis.
## save the plot
#pdf("FigureS06.pdf",width= 11, height = 8.5)

par(mfrow =c(1,1))
plot.new()
par(fig=c(0,1,0,1), new=TRUE,mar = c(6,5,1,0), oma = c(1,1,1,1), 
    cex.axis = 1)
#Pond A
plot(Summary_PondsFillUp$Median_WPERatio~Summary_PondsFillUp$MoYr_order, subset = Summary_PondsFillUp$Burkeetal2019_ID == "A"
     , pch = 22,cex =2,col = "black",bg= Pondcolors[1],ylim = c(0,160), xlim = c(1,19), axes = FALSE, xlab = " ", ylab = " ")
#Pond B
points(Summary_PondsFillUp$Median_WPERatio~Summary_PondsFillUp$MoYr_order, subset = Summary_PondsFillUp$Burkeetal2019_ID == "B"
       , pch = 16,cex =2,col= Pondcolors[2], ylim = c(0,160), xlim = c(1,19))
points(Summary_PondsFillUp$Median_WPERatio~Summary_PondsFillUp$MoYr_order, subset = Summary_PondsFillUp$Burkeetal2019_ID == "B"
       , pch = 10,cex =2,col= "black", xlim = c(1,7))
#Pond C
points(Summary_PondsFillUp$Median_WPERatio~Summary_PondsFillUp$MoYr_order, subset = Summary_PondsFillUp$Burkeetal2019_ID == "C"
       , pch = 24,cex =2,col= "black", bg= Pondcolors[3], ylim = c(0,160), xlim = c(1,19))
#Pond D
points(Summary_PondsFillUp$Median_WPERatio~Summary_PondsFillUp$MoYr_order, subset = Summary_PondsFillUp$Burkeetal2019_ID == "D"
       , pch = 23,cex =2,bg= Pondcolors[4], ylim = c(0,160), xlim = c(1,19))
# Pond E
points(Summary_PondsFillUp$Median_WPERatio~Summary_PondsFillUp$MoYr_order, subset = Summary_PondsFillUp$Burkeetal2019_ID == "E"
       , pch = 21,cex =2,col= "black",bg =Pondcolors[5],ylim = c(0,160), xlim = c(1,19))
#Pond F
points(Summary_PondsFillUp$Median_WPERatio~Summary_PondsFillUp$MoYr_order, subset = Summary_PondsFillUp$Burkeetal2019_ID == "F"
       , pch = 25,cex =2,col="black" ,bg= Pondcolors[6], ylim = c(0,160), xlim = c(1,19))
#Pond H
points(Summary_PondsFillUp$Median_WPERatio~Summary_PondsFillUp$MoYr_order, subset = Summary_PondsFillUp$Burkeetal2019_ID == "H"
       , pch = 22,cex =2,col= "black", bg= Pondcolors[7], ylim = c(0,160), xlim = c(1,19))
# #vertical error bars
arrows(x0=Summary_PondsFillUp$MoYr_order,
       y0=Summary_PondsFillUp$Quant25_WPERatio,
       x1=Summary_PondsFillUp$MoYr_order,
       y1=Summary_PondsFillUp$Quant75_WPERatio,
       code=3, lty=1, angle=90, length=0.1)
abline(h=100, lty= 2)
abline(v=4,lty =1)
abline(v=8,lty =1)
abline(v=12,lty=1)
abline(v=16,lty=1)

axis(2,at = c(25,50,75,100,125,150,175), las = 2, cex=1.25)
axis(1, at=c(2,6,9,10,11,13,14,15,17,18,19),lab = c("July","July","June","July",
                                                    "August","June","July","August","June","July","August"), las = 3)
mtext(side = 2, expression("Median Ratio of Water Area to Pond Depression Area (%)"), line = 3, outer = FALSE, cex = 1.25 )
mtext(side = 1, expression("Sampling Month"), line = 6, outer = FALSE, cex = 1.25)
mtext(side = 1, expression(italic("Total Precipitation")), line = 3.5, outer = FALSE, at = 0)
mtext(side = 1, expression(italic("(mm)")), line = 4.25, outer = FALSE, at = 0)
mtext(side = 1, expression(italic("48.1")), line = 4, outer = FALSE, at = 2)
mtext(side = 1, expression(italic("76.4")), line = 4, outer = FALSE, at = 6)
mtext(side = 1, expression(italic("62.0")), line = 4, outer = FALSE, at = 9)
mtext(side = 1, expression(italic("109.5")), line = 4, outer = FALSE, at = 10)
mtext(side = 1, expression(italic("35.2")), line = 4, outer = FALSE, at = 11)
mtext(side = 1, expression(italic("9.5")), line = 4, outer = FALSE, at = 13)
mtext(side = 1, expression(italic("106.9")), line = 4, outer = FALSE, at = 14)
mtext(side = 1, expression(italic("85.4")), line = 4, outer = FALSE, at = 15)
mtext(side = 1, expression(italic("80.2")), line = 4, outer = FALSE, at = 17)
mtext(side = 1, expression(italic("51.0")), line = 4, outer = FALSE, at = 18)
mtext(side = 1, expression(italic("56.4")), line = 4, outer = FALSE, at = 19)
text(2,160, expression("2014"), cex = 1.25)
text(6,160, expression("2015"), cex = 1.25)
text(10,160, expression("2016"), cex = 1.25)
text(14,160, expression("2017"), cex = 1.25)
text(18,160, expression("2018"), cex = 1.25)
box(lty=1)
mtext(side = 3, expression("Ponds:"), line = 0.75, outer = FALSE, at = -0.5, cex = 1.25)
legend(x = 0, y = 180, legend=c("A", "B", "C", "D", "E", "F", "H")
       ,pch = c(22,21,24,23,21,25,22),col = "black",pt.bg=Pondcolors, 
       cex=1.25, bty = "n", pt.cex = 2, horiz = TRUE, xpd = NA)
legend(x = 0, y = 180,legend = c("A", "B", "C", "D", "E", "F", "H")
       ,pch = c(22,10,24,23,21,25,22),col = "black", 
       cex=1.25, bty = "n", pt.cex = 2, horiz = TRUE, xpd = NA)




# Supplementary Figure S7#########################################################################        
# Plot Description: This plot is an expansion of Figure 6A with each pond plotted
# separately. Each sampling season is colored coded. If a kruskal() was 
# significant, an * is added to the subplot label and lowercase letters signify
# significant pairwise comparisons.
#
## save the plot
#pdf("FigureS07.pdf",width= 11, height = 8.5)

par(mfrow = c(1,1))
plot.new()
#Pond A
par(fig=c(.125,.375,.5,1), new=TRUE,mar = c(6,4,1,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~year, data = AllPondEdge_QuadFWing_Polygons, 
        subset = AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "A", 
        col = SamplingSeasons_colors, ylab = " ", xlab = " ",ylim = c(100,270), 
        axes = FALSE, outpch = 16)
mtext(side = 3,expression("Pond A"), cex = 1.5)
mtext(side = 1, expression(bolditalic("n = ")), line = 3.5, outer = FALSE, 
      at = 0.25, cex = 1.25 )
mtext(side = 1, expression(bolditalic("1")), line = 3.5, outer = FALSE, 
      at = 1, cex = 1.25)
mtext(side = 1, expression(bolditalic("1")), line = 3.5, outer = FALSE, 
      at = 2, cex = 1.25)
mtext(side = 1, expression(bolditalic("13")), line = 3.5, outer = FALSE, 
      at = 3, cex = 1.25)
mtext(side = 1, expression(bolditalic("2")), line = 3.5, outer = FALSE, 
      at = 4, cex = 1.25)
mtext(side = 1, expression(bolditalic("8")), line = 3.5, outer = FALSE, 
      at = 5, cex = 1.25)
axis(1, at = c(1,2,3,4,5), lab = c("2014","2015","2016","2017","2018"), 
     las = 2,cex.axis= 1.25)
axis(2, at = c(100,150,200,250), las =2,cex.axis= 1.25)
box(lty = 1, lwd = 1)

#kruskal(AllPondEdge_QuadFWing_PondA$area_m2,AllPondEdge_QuadFWing_PondA$year,
#        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# not significant
#
#Pond B
par(fig=c(.375,.625,.5,1), new=TRUE,mar = c(6,4,1,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~year, data = AllPondEdge_QuadFWing_Polygons, 
        subset = AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "B",
        col = SamplingSeasons_colors, ylab = " ", xlab = " ", ylim = c(0,120),
        axes = FALSE, outpch = 16)
text(1.5,115, expression("H = 17.6"), cex = 1.25)
text(3,115,expression("d.f. = 4"), cex = 1.25)
text(4.5,115, expression(italic("p")~"= 0.001"), cex = 1.25)
mtext(side = 3, expression("Pond B*"), cex = 1.5) 
text(1,32,expression("a"), cex=1.25)# 2014 
text(2,37,expression("a"), cex=1.25)# 2015 
text(3,41,expression("a"), cex=1.25)#2016 
text(4,43,expression("ab"), cex=1.25)#2017 
text(5,71,expression("b"), cex=1.25)# 2018 
axis(1, at = c(1,2,3,4,5), lab = c("2014","2015","2016","2017","2018"), las = 2,
     cex.axis= 1.25)
axis(2, at = c(0,20,40,60,80,100), las =2,cex.axis= 1.25)

mtext(side = 1, expression(bolditalic("1")), line = 3.5, outer = FALSE, at = 1,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("1")), line = 3.5, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("13")), line = 3.5, outer = FALSE, at = 3,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("2")), line = 3.5, outer = FALSE, at = 4,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("7")), line = 3.5, outer = FALSE, at = 5,
      cex = 1.25)
box(lty =1, lwd =1)
#
kruskal(AllPondEdge_QuadFWing_PondB$area_m2,AllPondEdge_QuadFWing_PondB$year,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
#significant
#
#Pond C
par(fig=c(.625,0.875,.5,1), new=TRUE,mar = c(6,4,1,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~year, data = AllPondEdge_QuadFWing_Polygons, 
        subset = AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID =="C", 
        outerpch = 16, col = SamplingSeasons_colors, ylab = " ", xlab = " ", 
        axes =  FALSE, ylim = c(0,120), outpch = 16)
mtext (side = 3, expression("Pond C"), cex = 1.5)
axis(1, at = c(1,2,3,4,5), lab = c("2014","2015","2016","2017","2018"), las= 2,
     cex.axis= 1.25)
axis(2, at = c(0,20,40,60,80,100), las =2,cex.axis= 1.25)

mtext(side = 1, expression(bolditalic("1")), line = 3.5, outer = FALSE, at = 1, 
      cex = 1.25)
mtext(side = 1, expression(bolditalic("1")), line = 3.5, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("13")), line = 3.5, outer = FALSE, at = 3,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("5")), line = 3.5, outer = FALSE, at = 4, 
      cex = 1.25)
mtext(side = 1, expression(bolditalic("8")), line = 3.5, outer = FALSE, at = 5,
      cex = 1.25)
box(lty = 1, lwd =1)
#
#kruskal(AllPondEdge_QuadFWing_PondC$area_m2,AllPondEdge_QuadFWing_PondC$year,
#        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# not significant
#
# Pond D
par(fig=c(0.02,.27,0,.5), new=TRUE,mar = c(6,4,1,0), oma = c(1,1,1,1),
    cex.axis = 1)
boxplot(area_m2~year, data = AllPondEdge_QuadFWing_Polygons, 
        subset = AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "D",
        col = SamplingSeasons_colors[3:5], ylab = " ", xlab = " ", 
        ylim =c(400,630), axes = FALSE)
mtext(side = 3, expression("Pond D*"), cex = 1.5)
mtext(side = 2, expression ("Pond Depression Area ("*~m^{2}*")"), line = 4, 
      at = 700 , cex = 1.5)
#mtext(side = 2, expression (italic("from Quadcopter and Fixed Wing Imagery")), line = 3, 
#     at = 700 , cex = 1.35)
text(1,620, expression("H = 12.4"), cex = 1.25)
text(2,620, expression("d.f. = 2"), cex = 1.25)
text(3,620, expression(italic("p")~"= 0.002"), cex = 1.25)
text(1,490, expression("a"), cex = 1.25) # 2016 
text(2,515, expression("b"), cex = 1.25) # 2017 
text(3,545, expression("b"), cex = 1.25) # 2018 
axis(1, at = c(1,2,3), lab = c("2016","2017","2018"), las= 2,cex.axis= 1.25)
axis(2, at = c(400,450,500,550,600), las =2,cex.axis= 1.25)
mtext(side = 1, expression( bolditalic("n = ")), line = 3.5, outer = FALSE, 
      at = 0.25, cex = 1.25)
mtext(side = 1, expression(bolditalic("10")), line = 3.5, outer = FALSE, 
      at = 1, cex = 1.25)
mtext(side = 1, expression(bolditalic("2")), line = 3.5, outer = FALSE, 
      at = 2, cex = 1.25)
mtext(side = 1, expression(bolditalic("5")), line = 3.5, outer = FALSE, 
      at = 3, cex = 1.25)
box(lty = 1, lwd =1)
#
kruskal(AllPondEdge_QuadFWing_PondD$area_m2,AllPondEdge_QuadFWing_PondD$year,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
#significant
#
# Pond E
par(fig=c(.27,.52,0,.5), new=TRUE,mar = c(6,4,1,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~year, data = AllPondEdge_QuadFWing_Polygons, 
        subset = AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "E",
        col = SamplingSeasons_colors, ylab = " ", xlab = " ", ylim = c(50,220), 
        axes = FALSE, outpch = 16)
mtext(side = 3, expression ("Pond E*"), cex = 1.5)
text(1.5,215, expression("H = 18.7"), cex = 1.25)
text(3,215, expression("d.f. = 4"), cex = 1.25)
text(4.5,215, expression(italic("p")~"= 0.0009"), cex = 1.25)
axis(1, at = c(1,2,3,4,5), lab = c("2014","2015","2016","2017","2018"), las= 2,
     cex.axis= 1.25)
axis(2, at = c(50,100,150,200), las =2,cex.axis= 1.25)
text(1,80,expression("a"), cex=1.25)#2014 
text(2,120,expression("ab"), cex=1.25)#2015 
text(3,115,expression("a"), cex=1.25)#2016 
text(4,160,expression("b"), cex=1.25)#2017 
text(5,175,expression("b"), cex=1.25)#2018

mtext(side = 1, expression(bolditalic("1")), line = 3.5, outer = FALSE, at = 1, 
      cex = 1.25)
mtext(side = 1, expression(bolditalic("1")), line = 3.5, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("10")), line = 3.5, outer = FALSE, at = 3,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("5")), line = 3.5, outer = FALSE, at = 4, 
      cex = 1.25)
mtext(side = 1, expression(bolditalic("7")), line = 3.5, outer = FALSE, at = 5, 
      cex = 1.25)
box(lty =1, lwd =1)
#
kruskal(AllPondEdge_QuadFWing_PondE$area_m2,AllPondEdge_QuadFWing_PondE$year,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# significant
#
#Pond F
par(fig=c(.52,.76,0,.5), new=TRUE,mar = c(6,4,1,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~year, data = AllPondEdge_QuadFWing_Polygons, 
        subset = AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "F",
        col = SamplingSeasons_colors, ylab = " ", xlab = " ", ylim = c(50,220), 
        axes = FALSE)
mtext(side = 3, expression("Pond F"), cex = 1.5)
axis(1, at = c(1,2,3,4,5), lab = c("2014","2015","2016","2017","2018"), 
     las = 2,cex.axis= 1.25)
axis(2, at = c(50,100,150,200), las =2,cex.axis= 1.25)
mtext(side = 1, expression ("Year"), line = 5.5, at = -1, cex = 1.5)

mtext(side = 1, expression(bolditalic("1")), line = 3.5, outer = FALSE, at = 1,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("1")), line = 3.5, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("13")), line =3.5, outer = FALSE, at = 3,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("2")), line = 3.5, outer = FALSE, at = 4,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("7")), line = 3.5, outer = FALSE, at = 5,
      cex = 1.25)
box(lty =1, lwd =1)
#
kruskal(AllPondEdge_QuadFWing_PondF$area_m2,AllPondEdge_QuadFWing_PondF$year,
        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# not significant
#
# Pond H
par(fig=c(.75,1,0,.5), new=TRUE,mar = c(6,4,1,0), oma = c(1,1,1,1), 
    cex.axis = 1)
boxplot(area_m2~year, data = AllPondEdge_QuadFWing_Polygons, 
        subset = AllPondEdge_QuadFWing_Polygons$Burkeetal2019_ID == "H",  
        col = SamplingSeasons_colors, ylab = " ", xlab = " ", ylim = c(100,270),
        axes = FALSE)
mtext(side = 3, expression ("Pond H"), cex = 1.5)
axis(1, at = c(1,2,3,4,5), lab = c("2014","2015","2016","2017","2018"),las = 2,
     cex.axis= 1.25)
axis(2, at = c(100,150,200,250), las =2,cex.axis= 1.25)

mtext(side = 1, expression(bolditalic("1")), line = 3.5, outer = FALSE, at = 1, 
      cex = 1.25)
mtext(side = 1, expression(bolditalic("1")), line = 3.5, outer = FALSE, at = 2,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("12")), line = 3.5, outer = FALSE, at = 3,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("4")), line = 3.5, outer = FALSE, at = 4,
      cex = 1.25)
mtext(side = 1, expression(bolditalic("7")), line = 3.5, outer = FALSE, at = 5,
      cex = 1.25)
box(lty =1, lwd = 1)
#
#kruskal(AllPondEdge_QuadFWing_PondH$area_m2,AllPondEdge_QuadFWing_PondH$year,
#        alpha = 0.05, p.adj = c("bonferroni"), group = TRUE, console = TRUE)
# not significant

#dev.off()
#

#######
# Supplementary Figure S8 ####
# Plot Description: This set of subplots presents the median daily ebullitive flux
# (mg CH4 m2 d-1) compared to Water Area (m^2) from Quadcopter Imagery collected
# in 2016, 2017, 2018. All ponds are plotted together in the top left plot, followed
# by each pond individually (differentiated by color).
## save the plot
#pdf("FigureS08.pdf",width= 11, height = 8.5)

par(mfrow = c(1,1))
plot.new()
par(mar = c(8,8,8,8))
par(fig=c(0.02,.27,.5,1), new=TRUE,mar = c(4,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, subset = Quad_PondEdge_PondPolygons$Burkeetal2019_ID == "A", pch = 21, 
     col = "black", bg = Pondcolors [1], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,400), cex = 1)
points(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, subset = Quad_PondEdge_PondPolygons$Burkeetal2019_ID == "B", pch = 21, 
       col = "black", bg = Pondcolors [2], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,400), cex = 1)
points(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, subset = Quad_PondEdge_PondPolygons$Burkeetal2019_ID == "C", pch = 21, 
       col = "black", bg = Pondcolors [3], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,400), cex = 1)
points(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, subset = Quad_PondEdge_PondPolygons$Burkeetal2019_ID == "D", pch = 21, 
       col = "black", bg = Pondcolors [4], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,400), cex = 1)
points(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, subset = Quad_PondEdge_PondPolygons$Burkeetal2019_ID == "E", pch = 21, 
       col = "black", bg = Pondcolors [5], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,400), cex = 1)
points(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, subset = Quad_PondEdge_PondPolygons$Burkeetal2019_ID == "F", pch = 21, 
       col = "black", bg = Pondcolors [6], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,400), cex = 1)
points(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, subset = Quad_PondEdge_PondPolygons$Burkeetal2019_ID == "H", pch = 21, 
       col = "black", bg = Pondcolors [7], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,400), cex = 1)

axis(1, at = c(0,100,200,300,400), cex.axis = 1.25)
axis(2, at = c(0,100,200,300,400), las = 2, cex.axis = 1.25)
mtext(side = 3, expression("All Ponds*"),cex = 1.5)
text(100,375 , expression(tau~" = 0.15"), cex = 1.25)
text(300, 375, expression(italic("p")~"= 0.03"), cex = 1.25)
mtext(side = 2, expression ("Median Daily Ebullitive Flux ("*mg*~CH[4]*~m^{2}*~d^{-1}*")"),
      line = 3, at = -20, cex = 1.5)

box(lty =1, lwd = 1)
# the cor.test provides a tau value and p value
cor.test( Quad_Water_PondPolygons$Median8dC_BubFlux,
          Quad_Water_PondPolygons$area_m2,method = "kendall")

par(fig=c(0.26,.51,.5,1), new=TRUE,mar = c(4,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
# Plot A
plot(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, 
     subset = Quad_Water_PondPolygons$Burkeetal2019_ID == "A", pch = 21,
     col = "black", bg=rgb(254,224,144, maxColorValue = 255), ylab = " ", 
     xlab = " ",
     axes = FALSE, xlim = c(0,400), ylim = c(0,1))
axis(1,at =c(0,100,200,300,400),cex.axis= 1.25)
axis(2,at = c(0,0.25,0.5,.75,1), las = 2,cex.axis= 1.25)
mtext(side = 3, expression("Pond A"), cex = 1.5)
box(lty = 1,lwd =1)

#cor.test(Quad_Water_PondPolygons_PondA$Median8dC_BubFlux,
#         Quad_Water_PondPolygons_PondA$area_m2,method = "kendall")
#  NOT significant
#
# Pond B
par(fig=c(0.5,.75,.5,1), new=TRUE,mar = c(4,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)

plot(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, 
     subset = Quad_Water_PondPolygons$Burkeetal2019_ID == "B", pch = 21,
     col = "black", bg=rgb(224,243,248, maxColorValue = 255), ylab = " ",
     xlab = " ", axes = FALSE,
     ylim = c(0,1), xlim = c(0,100))
mtext(side = 3, expression("Pond B"), cex = 1.5)
axis(1,at =c(0,50,100),cex.axis= 1.25)
axis(2,at = c(0,0.25,0.5,0.75,1),las =2,cex.axis= 1.25)
box(lty = 1,lwd =1)
# 
#cor.test(Quad_Water_PondPolygons_PondB$Median8dC_BubFlux,
#         Quad_Water_PondPolygons_PondB$area_m2,method = "kendall")
#  NOT significant 
#
# Pond C
par(fig=c(.75,1,.5,1), new=TRUE,mar = c(4,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, 
     subset = Quad_Water_PondPolygons$Burkeetal2019_ID == "C", pch = 21,
     col = "black", bg=rgb(69,117,180, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, xlim = c(0,100), ylim = c(0,20))
axis(1,at =c(0,25,50,75,100),cex.axis= 1.25)
axis(2,at = c(0,5,10,15,20), las = 2,cex.axis= 1.25)
mtext(side = 3, expression("Pond C"), cex = 1.5)
box(lty = 1,lwd =1)
#
cor.test(Quad_Water_PondPolygons_PondC$Median8dC_BubFlux,
         Quad_Water_PondPolygons_PondC$area_m2,method = "kendall")
# Not significant
#
#Pond D
par(fig=c(0.02,.27,0,.5), new=TRUE,mar = c(5,4,0,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, 
     subset = Quad_Water_PondPolygons$Burkeetal2019_ID == "D", pch = 21,
     col = "black", bg=rgb(253,174,97, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, xlim = c(0,200), ylim = c(0,5))
axis(1,at =c(0,50,100,150,200),cex.axis= 1.25)
axis(2,at = c(0,1,2,3,4,5), las = 2,cex.axis= 1.25)
mtext(side = 3, expression("Pond D"), cex = 1.5)
box(lty = 1,lwd =1)
#
#cor.test(Quad_Water_PondPolygons_PondD$Median8dC_BubFlux,
#         Quad_Water_PondPolygons_PondD$area_m2,method = "kendall")
# NOT signifigant
#
#PondE 
par(fig=c(0.26,.51,0,.5), new=TRUE,mar = c(5,4,0,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, 
     subset = Quad_Water_PondPolygons$Burkeetal2019_ID == "E", pch = 21,
     col = "black", bg=rgb(215,48,39, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, ylim = c(0,150), xlim= c(0,200))
axis(1,at =c(0,50,100,150,200),cex.axis= 1.25)
axis(2,at = c(0,25,50,75,100,125,150), las = 2,cex.axis= 1.25)
mtext(side = 3, expression("Pond E"), cex =1.5)
box(lty = 1,lwd =1)
#
#cor.test(Quad_Water_PondPolygons_PondE$Median8dC_BubFlux,
#         Quad_Water_PondPolygons_PondE$area_m2,method = "kendall")# not significant 
# NOT significant
# Pond F
par(fig=c(.5,.75,0,.5), new=TRUE,mar = c(5,4,0,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, 
     subset = Quad_Water_PondPolygons$Burkeetal2019_ID == "F", pch = 21,
     col = "black", bg=rgb(244,109,67, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, ylim = c(0,400), xlim =c(0,150))
axis(1,at =c(0,50,100,150),cex.axis= 1.25)
axis(2,at = c(0,100,200,300,400), las = 2,cex.axis= 1.25)
mtext(side = 3, expression("Pond F"), cex = 1.5)
box(lty = 1,lwd =1)
mtext(side = 1, expression ("Water Area ("*m^{2}*")"), line = 3,at = -10 , cex = 1.5)
mtext(side = 1, expression (italic("from Quadcopter Imagery")), line = 4, at = -20, cex = 1.25)


#
cor.test(Quad_Water_PondPolygons_PondF$Median8dC_BubFlux,
         Quad_Water_PondPolygons_PondF$area_m2,method = "kendall")
# NOT significant 
#
#PondH
#par(fig=c(.62,.87,0,.5), new=TRUE,mar = c(5,4,2,0), oma = c(1,1,1,1), cex.axis = 1)
par(fig=c(0.75,1,0,.5), new=TRUE,mar = c(5,4,0,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(Median8dC_BubFlux~area_m2, data = Quad_Water_PondPolygons, 
     subset = Quad_Water_PondPolygons$Burkeetal2019_ID == "H",pch = 21,
     col = "black", bg=rgb(116,173,209, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, ylim = c(0,75), xlim = c(0,300))
axis(1,at =c(0,100,200,300),cex.axis= 1.25)
axis(2,at = c(0,25,50,75), las = 2,cex.axis= 1.25)
mtext(side = 3, expression("Pond H*"), cex = 1.5)
box(lty = 1,lwd =1)
text(75,70, expression(tau~" = -0.49"), cex = 1.25)
text(225,70, expression(italic("p")~"= 0.004"), cex = 1.25)

#
cor.test(Quad_Water_PondPolygons_PondH$Median8dC_BubFlux,
         Quad_Water_PondPolygons_PondH$area_m2,method = "kendall")
# significant
#
#dev.off()

######

# Supplementary Figure S9####
#Plot Description: This plot shows the relationship between Total Precip (mm) in
# the 8 days before UAS flight (x axis) and Median Daily ebulltive flux 
# (mg CH4 m2 d-1) across the whole study period. There are 8 subplots, the top
# left plot shows all ponds plotted together, distinguished by color, and the
# following plots are each pond plotted on its own. Ponds that showed a significant
# relationship have a * after their name and the tau and p values are displayed
# at the top of the plot.
## save the plot
#pdf("FigureS09.pdf",width= 11, height = 8.5)

par(mfrow = c(1,1))
plot.new()
par(mar = c(5,5,5,5))
par(fig=c(0.02,.27,.5,1), new=TRUE,mar = c(3,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
# All Ponds
plot(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons, 
     subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "A",pch = 21, 
     col = "black", bg = Pondcolors[1], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
points(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons, 
       subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "B",pch = 21, 
       col = "black", bg = Pondcolors[2], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
points(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons, 
       subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "C",pch = 21, 
       col = "black", bg = Pondcolors[3], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
points(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons, 
       subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "D",pch = 21, 
       col = "black", bg = Pondcolors[4], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
points(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons, 
       subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "E",pch = 21, 
       col = "black", bg = Pondcolors[5], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
points(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons, 
       subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "F",pch = 21, 
       col = "black", bg = Pondcolors[6], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
points(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons, 
       subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "H",pch = 21, 
       col = "black", bg = Pondcolors[7], ylim = c(0,400),ylab = " ",xlab = " ",axes = FALSE, xlim =c (0,100), cex = 1)
axis(1, at = c(0,20,40,60,80,100), cex.axis = 1.25)
axis(2, at = c(0,100,200,300,400), las = 2, cex.axis = 1.25)
mtext(side = 3, expression("All Ponds"),cex = 1.5)
mtext(side = 2, expression ("Median Daily Ebullitive Flux ("*mg*~CH[4]~m^{2}*~d^{-1}*")"), line = 4, at = -20, cex = 1.5)
box(lty =1, lwd = 1)
# the cor.test provides a tau value and p value
#cor.test(AllWater_QuadFWing_Polygons$Median8dC_BubFlux,
#        AllWater_QuadFWing_Polygons$TotPrec_8dpreflight,method = "kendall")
#not significant
par(fig=c(0.26,.51,.5,1), new=TRUE,mar = c(3,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
# Plot A
plot(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons, 
     subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "A", pch = 21,
     col = "black", bg=rgb(254,224,144, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, xlim = c(0,100), ylim = c(0,1))
axis(1,at =c(0,25,50,75,100),cex.axis= 1.25)
axis(2,at = c(0,.25,0.5,.75,1), las = 2,cex.axis= 1.25)
mtext(side = 3, expression("Pond A"), cex = 1.5)
box(lty = 1,lwd =1)

#cor.test(AllWater_QuadFWing_Polygons_PondA$Median8dC_BubFlux,
#         AllWater_QuadFWing_Polygons_PondA$TotPrec_8dpreflight,method = "kendall")
#   not significant.
#
# Pond B
par(fig=c(0.5,.75,.5,1), new=TRUE,mar = c(3,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons, 
     subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "B", pch = 21,
     col = "black", bg=rgb(224,243,248, maxColorValue = 255), ylab = " ", xlab = " ", axes = FALSE,
     ylim = c(0,1), xlim = c(0,100))
mtext(side = 3, expression("Pond B"), cex = 1.5)
axis(1,at =c(0,25,50,75,100),cex.axis= 1.25)
axis(2,at = c(0,.25,.5,.75,1),las =2,cex.axis= 1.25)
box(lty = 1,lwd =1)

#cor.test(AllWater_QuadFWing_Polygons_PondB$Median8dC_BubFlux,
#         AllWater_QuadFWing_Polygons_PondB$TotPrec_8dpreflight,method = "kendall")
#  not significant

#
# Pond C
par(fig=c(.75,1,.5,1), new=TRUE,mar = c(3,3,0.5,1), oma = c(1,1,1,1), 
    cex.axis = 1)

plot(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons,
     subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "C", pch = 21,
     col = "black", bg=rgb(69,117,180, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, xlim = c(0,100), ylim = c(0,10))
axis(1,at =c(0,25,50,75,100),cex.axis= 1.25)
axis(2,at = c(0,2.5, 5, 7.5, 10), las = 2,cex.axis= 1.25)
mtext(side = 3, expression("Pond C"), cex = 1.5)
box(lty = 1,lwd =1)
#
#cor.test(AllWater_QuadFWing_Polygons_PondC$Median8dC_BubFlux,
#         AllWater_QuadFWing_Polygons_PondC$TotPrec_8dpreflight,method = "kendall")
# not significant
#
#Pond D
par(fig=c(0.02,.27,0,.5), new=TRUE,mar = c(3,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons,
     subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "D", pch = 21,
     col = "black", bg=rgb(253,174,97, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, xlim = c(0,100), ylim = c(0,10))
axis(1,at =c(0,25,50,75,100),cex.axis= 1.25)
axis(2,at = c(0,2.5,5,7.5,10), las = 2,cex.axis= 1.25)
mtext(side = 3, expression("Pond D"), cex = 1.5)
box(lty = 1,lwd =1)
#
#cor.test(AllWater_QuadFWing_Polygons_PondD$Median8dC_BubFlux,
#         AllWater_QuadFWing_Polygons_PondD$TotPrec_8dpreflight,method = "kendall")
# not significant
#
#PondE 
par(fig=c(0.26,.51,0,.5), new=TRUE,mar = c(3,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons,
     subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "E", pch = 21,
     col = "black", bg=rgb(215,48,39, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, ylim = c(0,150), xlim= c(0,100))
axis(1,at =c(0,25,50,75,100),cex.axis= 1.25)
axis(2,at = c(0,50,100,150), las = 2,cex.axis= 1.25)
mtext(side = 3, expression("Pond E"), cex =1.5)
box(lty = 1,lwd =1)
#
#cor.test(AllWater_QuadFWing_Polygons_PondE$Median8dC_BubFlux,
#         AllWater_QuadFWing_Polygons_PondE$TotPrec_8dpreflight,method = "kendall")
# not significant 
#
# Pond F
par(fig=c(.5,.75,0,.5), new=TRUE,mar = c(3,4,0.5,0), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons,
     subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "F", pch = 21,
     col = "black", bg=rgb(244,109,67, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, ylim = c(0,400), xlim =c(0,100))
axis(1,at =c(0,20,40,60,80,100),cex.axis= 1.25)
axis(2,at = c(0,100,200,300,400), las = 2,cex.axis= 1.25)
mtext(side = 3, expression("Pond F"), cex = 1.5)
box(lty = 1,lwd =1)
mtext(side = 1, expression ("Total Precipitation (mm) in the 8 days before Flight"),
      line = 3,at = -10 , cex = 1.5)
#
#cor.test(AllWater_QuadFWing_Polygons_PondF$Median8dC_BubFlux,
#         AllWater_QuadFWing_Polygons_PondF$TotPrec_8dpreflight,method = "kendall")
# not significant
#
#PondH
par(fig=c(0.75,1,0,.5), new=TRUE,mar = c(3,3,0.5,1), oma = c(1,1,1,1), 
    cex.axis = 1)
plot(Median8dC_BubFlux~TotPrec_8dpreflight, data = AllWater_QuadFWing_Polygons,
     subset = AllWater_QuadFWing_Polygons$Burkeetal2019_ID == "H",pch = 21,
     col = "black", bg=rgb(116,173,209, maxColorValue = 255), ylab = " ", xlab = " ",
     axes = FALSE, ylim = c(0,75), xlim = c(0,100))
axis(1,at =c(0,25,50,75,100),cex.axis= 1.25)
axis(2,at = c(0,25,50,75), las = 2,cex.axis= 1.25)
mtext(side = 3, expression("Pond H*"), cex = 1.5)
text(20,70 , expression(tau~" = 0.31"), cex = 1.25)
text(75,70, expression(italic("p")~"= 0.052"), cex = 1.25)
box(lty = 1,lwd =1)
#
cor.test(AllWater_QuadFWing_Polygons_PondH$Median8dC_BubFlux,
         AllWater_QuadFWing_Polygons_PondH$TotPrec_8dpreflight,method = "kendall")
# mildly significant
#dev.off()

#####
#####