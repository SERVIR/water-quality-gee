import ee
from ee.ee_exception import EEException
import datetime

####Still under development#########
try:
    ee.Initialize()
except EEException as e:
    from oauth2client.service_account import ServiceAccountCredentials
    credentials = ServiceAccountCredentials.from_p12_keyfile(
    service_account_email='',
    filename='',
    private_key_password='notasecret',
    scopes=ee.oauth.SCOPE + ' https://www.googleapis.com/auth/drive ')
    ee.Initialize(credentials)

geometry = ee.Geometry.Polygon([[[78.41079711914062, 17.465297700926097],
          [78.41629028320312, 17.334252606951086],
          [78.59893798828125, 17.33490806580805],
          [78.55087280273438, 17.494769882318828]]])

# This script processes and charts WQ parameter values from a collection of MODIS L3 OC archive for a particular pixel throughout time.
# WQ parameters such as chlorophyll-a, secchi depth and the trophic state index
# How to use:
# (1) If a "geometry" iable exists in the imports window, delete it.
# (2) Within the map, select the polygon button near the top left side.
# (3) Create a new polygon by drawing around the area of interest
# (4) Adjust the time frame you wish to browse by adding here:
# begin date
startDate = '2016-04-01'
# end date
endDate = '2016-10-31'
# (5) Click Run
# (5) The "available imagery" ImageCollection within the console displays all available imagery. If you wish to investigate a single
#     image from this collection, highlight the FILE_ID and copy and paste it into the proper location within the "Ls8 Single Image" script.
# (6) The charts will appear on the right hand side of the interface, and will display the mean map in the map user interface.
# (6) Click "Run"
# (7) Export each chart as a csv by clicking on the export icon near the top-right of the chart.
# (8) If you are interested in a single image from the collection, copy and paste the FILE_ID here:
singleImage = ee.Image('NASA/OCEANDATA/MODIS-Aqua/L3SMI/A2017183').clip(geometry)

# Author: Benjamin Page #
###################################################################/
# Import Colelction
MODIS = ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI')

# filter collection
FC = MODIS.filterDate(startDate, endDate).filterBounds(geometry)


###################################################################/
### Map Functions ###

def secchi(img):
    Rrs_488 = img.select('Rrs_488')
    Rrs_667 = img.select('Rrs_667')
    ln_blueRed = (Rrs_488.divide(Rrs_667)).log()
    lnMOSD = (ee.Image(1.4856).multiply(ln_blueRed)).add(ee.Image(0.2734))  # R2 = 0.8748 with Anthony's in-situ data
    MOSD = ee.Image(10).pow(lnMOSD)

    SD = (ee.Image(0.1777).multiply(MOSD)).add(ee.Image(1.0813))

    return (SD.set('system:time_start', img.get('system:time_start')))


def trophicState(img):
    TSI = ee.Image(60).subtract((ee.Image(14.41)).multiply((img.log())))
    return (TSI.set('system:time_start', img.get('system:time_start')))

### Models ###/

# chlor_a
chlor_a = FC.select('chlor_a')

# SD
SD = ee.ImageCollection(FC.map(secchi).copyProperties(FC, ['system:time_start']))

# TSI
TSI = SD.map(trophicState)





###################################################################/
# Time Series

# chlorTimeSeries = ui.Chart.image.seriesByRegion(
# chlor_a, geometry, ee.Reducer.mean())
# .setChartType('ScatterChart')
#     .setOptions({
#     title: 'Mean chlor_a',
#     vAxis: {title: 'chlor_a (ug / L)'},
#     lineWidth: 1,
#     pointSize: 4,
# })
#
# sdTimeSeries = ui.Chart.image.seriesByRegion(
# SD, geometry, ee.Reducer.mean())
# .setChartType('ScatterChart')
#     .setOptions({
#     title: 'Mean Secchi Depth',
#     vAxis: {title: 'Zsd [m]'},
#     lineWidth: 1,
#     pointSize: 4,
# })
#
# tsiTimeSeries = ui.Chart.image.seriesByRegion(
# TSI, geometry, ee.Reducer.mean())
# .setChartType('ScatterChart')
#     .setOptions({
#     title: 'Mean TSI',
#     vAxis: {title: 'TSI value)'},
#     lineWidth: 1,
#     pointSize: 4,
# })
#
# print(chlorTimeSeries)
# print(sdTimeSeries)
# print(tsiTimeSeries)

###################################################################/

# Map layers #

chlor_a_mean = chlor_a.mean()
# chlor_a_mean = chlor_a_mean.clip(geometry)
chlor_a_vis = {'min': 0, 'max': 100, 'palette': 'darkblue,blue,limegreen,yellow,orange,orangered,darkred'}
chlor_a_mapid = chlor_a_mean.getMapId(chlor_a_vis)

SD_mean = SD.mean()
# SD_mean = SD_mean.clip(geometry)
SD_vis =  {'min': 0, 'max': 2, 'palette': 'red,orangered,orange,yellow,limegreen,blue,darkblue'}
SD_mapid = SD_mean.getMapId(SD_vis)

TSI_mean = TSI.mean()
# TDI_mean = TSI_mean.clip(geometry)
TSI_vis =  {'min': 0, 'max': 2, 'palette': 'red,orangered,orange,yellow,limegreen,blue,darkblue'}
TSI_mapid = TSI_mean.getMapId(TSI_vis)

singleImageLyr = singleImage.select('chlor_a')
singleImage_vis = {'min': 0, 'max': 100, 'palette': 'darkblue,blue,limegreen,yellow,orange,orangered,darkred'}
singleImage_mapid = singleImageLyr.getMapId(singleImage_vis)

print "CHLOR A",chlor_a_mapid
print "SD",SD_mapid
print "TSI",TSI_mapid
print "Single", singleImage_mapid
# Map.addLayer(chlor_a.mean().clip(geometry),
#              {min: 0, max: 100, palette: ['darkblue', 'blue', 'limegreen', 'yellow', 'orange', 'orangered', 'darkred']},
#              'chlor_a', false)
# Map.addLayer(SD.mean().clip(geometry),
#              {min: 0, max: 2, palette: ['red', 'orangered', 'orange', 'yellow', 'limegreen', 'blue', 'darkblue']}, 'SD',
#              false)
# Map.addLayer(TSI.mean().clip(geometry),
#              {min: 30, max: 80, palette: ['darkblue', 'blue', 'limegreen', 'yellow', 'orange', 'orangered', 'red']},
#              'TSI', true)
# Map.addLayer(singleImage.select('chlor_a'),
#              {min: 0, max: 100, palette: ['darkblue', 'blue', 'limegreen', 'yellow', 'orange', 'orangered', 'darkred']},
#              'singleImage_chlor_a', false)

################################################################/
# Exporting Images #

# # Export mean chlor_a
# Export.image.toDrive({
#     image: chlor_a.mean().clip(geometry),
#     description: 'Mean Chlorophyll-a',
#     scale: 250,
#     region: geometry
# })
#
# # Export mean SD
# Export.image.toDrive({
#     image: SD.mean().clip(geometry),
#     description: 'Mean Secchi Depth',
#     scale: 250,
#     region: geometry
# })
#
# # Export mean TSI
# Export.image.toDrive({
#     image: TSI.mean().clip(geometry),
#     description: 'Mean Trophic State (from index)',
#     scale: 250,
#     region: geometry
# })