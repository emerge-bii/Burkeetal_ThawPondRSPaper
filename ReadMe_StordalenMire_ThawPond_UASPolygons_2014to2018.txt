Read Me file for StordalenMire_ThawPond_UASPolygons_2014to2018.csv
Comma (,) eliminated text file

Point of contact:
Dr. Sophia Burke
sophia.burke@unh.edu
sophieaburke@gmail.com


Details on Sampling and Data Processing Methods can be found in the following manuscript:
Burke et al. (2022) entitled Connecting Methane Ebullitive Flux to Thaw Pond size using Unpiloted Aerial Systems.

Note: Reference to Pond Edge Area is the same as Pond Depression Area in the manuscript.

1 header row, 11 headings: 
FlightDate,Field_ID,Burkeetal2019_ID,UAStype,PolygonType,PondType,area_m2,edge_m,edge.area,8dMedianC_BubFlux,Cumulative_BubFlux

Column Headings:
FlightDate		-- YYYY-MM-DD, time zone: CET
Field_ID		-- original pond identifier used in field notes
Burkeetl2019_ID 	-- Pond ID used in figures in this paper, following the ID 
	                   used in Burke et al. (2019)
UAStype			-- Either 'FWing' for Fixed Wing, or 'Quad' for Quadcopter
PolygonType		-- Either 'pondedge' for Pond Edge Polygon, or 'water' for Water Polygon
PondType		-- Either 1,2,3 or 4. See Burke et al. (2019) for complete description of each type
area_m2			-- area of the polygon (m2) -- calculated in QGIS
edge_m			-- perimeter (edge) length of the polygon (m) -- calculated in QGIS
edge.area		-- edge divided by area of the polygon -- calculated in QGIS
Median8dC_BubFlux	-- Median daily bubbly flux (mg CH4 m-2 d-1) from each thaw pond calculated across an eight day window centered on each flight date.
Cumulative_BubFlux	-- Cumulative CH4 (mg CH4 m-2) emitted via ebullition calculated for that sampling season

TotPrec_8dpreflight 	-- Total precipitation measured in the eight days leading up to and including the UAS flight date, in mm.	

NOTE: Reference to specific ponds in this ReadMe file use the Burkeetal(2019)_ID identifier.

Detailed Description:

FlightDate: corresponds to the date during which the UAS (either FWing or Quad) were flown. On average a UAS was flown over each pond every eight days. UAS (both FWing and Quad) were only flown if the weather cooperated: no rain, low wind). More 		detailed spread of the imagery collection:
		2014: FWing was flown once over each pond in July
		2015: FWing was flown once over each pond in July
		2016: Quad was flown 13 days during the field season (~ 6 days apart), FWing was flown once over each pond in July
		2017: Quad was flown 4 days during the field season (~14 days apart)
		2018: Quad was flown 9 days during the field season (~ 7 days apart)

Pond (Field_ID and Burkeetal2019_ID): seven ponds were the focus of this study. see Burke et al. (2019) for a complete description of each pond.

GPS Locations (Lat,Long) of each pond -- the following GPS coordinates are of a point on the shore of each pond, no more than 3 meters from the pond itself.
	A -- 68.35698483,19.04960763
	B -- 68.35688503,19.05062797
	C -- 68.35344712,19.04727789
	D -- 68.35773563,19.05368401
	E -- 68.35270479,19.04733059
	G -- 68.35428965,19.04707764
	H -- 68.35391821,19.04773924

UAStype: either 'FWing' for Fixed-Wing UAS or 'Quad' for Quadcopter UAS. see Burke et al. (2022) for complete description of UAS imagery acquisition and processing.
	FWing: Triton XL, Robota, Lancaster, TX, with a built in camera, Goose™ autopilot program, Robota. Flown at an altitude of 70 m across the full extent of Stordalen Mire. see Palace et al., 2018 for complete methodology

	Quadcopter: Yuneec® Q500 quadcopter UAS, a gimbal RGB camera. Flown at an altitude of 15 - 20 m

PolygonType: refers to which polygon the area/edge/edge.area values correspond to, see Burke et al. (2020) for more details
	pondedge: represents the extent of thaw for each pond
	water: represents where the water collected in each pond on that particular flight date

PondType: This numerical value (1 to 4) represents the pond type described in Burke et al.(2019). Pond types were developed to help explain the statistical differences between daily ebullitive flux of the eight ponds measured in Burke et al. (2019), which fell into four groups. Pond types vary from one another based on vegetation dominance, average pond depth, and hydrologic connectivity (or lack thereof). See table 1 in Burke et al. (2019) and table 1 in Burke et al. (2022) for a complete description of each pond type.

area_m2: The area of the polygon drawn in QGIS, in m2.

edge_m: The length of the perimeter (edge) of the polygon drawn in QGIS, in m.

edge.area: The ratio of the edge_m/area_m2 for the polygon drawn in QGIS, unitless.

8dMedianC_BubFlux: This is the median value of daily average ebullitive flux (mg CH4 m-2 d-1) measured in the eight days centered around each flight date. More specifically, it is the median value of the three days before the flight, the flight date, and the four days following. This was calculated using the ave() function, along with the rollmedian(align = "center") functions in R.

Cumulative_BubFlux: This is the cumulative (total) ebullitive emission per m2 measured from each pond during the field season associated with each UAS flight.
 
TotPrec_8dpreflight : This was calculated from 10 minute to 30 minute total precipitation data measured by ICOS-Sweden (Integrated Carbon Observation System) at Stordalen Mire, Abisko, Sweden using a WeatherHawk instrumentation system placed 4 m a.g.l on the top of an instrumentation building. https://www.icos-sweden.se/station_stordalen.html 


Other helpful sources providing details on sampling and data processing:

Burke, S.A., Wik, M., Lang, A., Contosta, A.R., Palace, M., Crill, P.M., Varner, R.K., 2019. Long‐Term Measurements of Methane Ebullition From Thaw Ponds. J. Geophys. Res. Biogeosci. 2018JG004786. https://doi.org/10.1029/2018JG004786

Palace, M., Herrick, C., DelGreco, J., Finnell, D., Garnello, A., McCalley, C., McArthur, K., Sullivan, F., Varner, R., 2018. Determining Subarctic Peatland Vegetation Using an Unmanned Aerial System (UAS). Remote Sensing 10, 1498. https://doi.org/10.3390/rs10091498


