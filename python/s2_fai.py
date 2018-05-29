import ee
from ee.ee_exception import EEException
import datetime

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

# geometry = ee.Geometry.Polygon([[[78.41079711914062, 17.465297700926097],
#           [78.41629028320312, 17.334252606951086],
#           [78.59893798828125, 17.33490806580805],
#           [78.55087280273438, 17.494769882318828]]])

geometry = ee.Geometry.Polygon([[[30.76171875, 0.9049611504960419],
          [30.8935546875, -3.487377195492663],
          [35.5517578125, -3.2680324702882952],
          [35.5517578125, 1.9593043032313748]]])

# This script processes a single Sentinel-2 TOA to Rayleigh corrected reflectances.
# The floating algal index (FAI) is calculated from Rrc over water bodies.
# How to use:
# (1) If a "geometry" iable exists in the imports window, delete it.
# (2) Within the map, select the point button near the top left side.
# (3) Create a new point by clicking on a location of interest and click run.
# (4) The "available imagery" ImageCollection within the console displays all available imagery
#     over that location with the necessary filters described by the user below.
# (5) Expand the features within the "available imagery" image collection and select an image
#     by highlighting the FILE_ID and pasting it here:
s2 = ee.Image('COPERNICUS/S2/20180216T075011_20180216T080852_T36MXE')
# (6) Click "Run"
# (7) Export the image to your Google Drive by clicking on the "tasks" tab and clicking "RUN", be sure to specify
#     the proper folder.

# Author: Benjamin Page #
# Citations:
# Page, B.P., Kumar, A. and Mishra, D.R., 2018. A novel cross-satellite based assessment of the spatio-temporal development of a cyanobacterial harmful algal bloom. International Journal of Applied Earth Observation and Geoinformation, 66, pp.69-81.

#########################################################/
#########################################################/
#########################################################/
# User Input #

# begin date
iniDate = '2015-01-01'

# end date
endDate = '2018-02-28'

cloudPerc = 5

#########################################################/
#########################################################/
#########################################################/
# Import Collections #

# sentinel-2
MSI = ee.ImageCollection('COPERNICUS/S2')

# toms / omi
ozone = ee.ImageCollection('TOMS/MERGED')

SRP = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')

#########################################################/
#########################################################/
#########################################################/
# water mask
startMonth = 5
endMonth = 9
startYear = 2013
endYear = 2017

forMask = SRP.filterBounds(geometry).select('B6').filterMetadata('CLOUD_COVER', "less_than", 10).filter(
    ee.Filter.calendarRange(startMonth, endMonth, 'month')).filter(ee.Filter.calendarRange(startYear, endYear, 'year'))
mask = ee.Image(forMask.select('B6').median().lt(300))
mask = mask.updateMask(mask)

# constants
pi = ee.Image(3.141592)
imDate = s2.date()

# filter sentinel 2 collection
FC = MSI.filterDate(iniDate, endDate).filterBounds(geometry).filterMetadata('CLOUDY_PIXEL_PERCENTAGE', "less_than",
                                                                            cloudPerc)


#########################################################/
#########################################################/
#########################################################/
# MSI Atmospheric Correction #

bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']

# rescale
rescale = ee.Image(s2.divide(10000).copyProperties(s2)).select(bands)

# tile footprint
footprint = s2.geometry()

# dem
DEM = ee.Image('USGS/SRTMGL1_003').clip(footprint)

# ozone
DU = ee.Image(ozone.filterDate(iniDate, endDate).filterBounds(geometry).mean())

# Julian Day
imgDate = ee.Date(s2.get('system:time_start'))
FOY = ee.Date.fromYMD(imgDate.get('year'), 1, 1)
JD = imgDate.difference(FOY, 'day').int().add(1)

# earth-sun distance
myCos = ((ee.Image(0.0172).multiply(ee.Image(JD).subtract(ee.Image(2)))).cos()).pow(2)
cosd = myCos.multiply(pi.divide(ee.Image(180))).cos()
d = ee.Image(1).subtract(ee.Image(0.01673)).multiply(cosd).clip(footprint)

# sun azimuth
SunAz = ee.Image.constant(s2.get('MEAN_SOLAR_AZIMUTH_ANGLE')).clip(footprint)

# sun zenith
SunZe = ee.Image.constant(s2.get('MEAN_SOLAR_ZENITH_ANGLE')).clip(footprint)
cosdSunZe = SunZe.multiply(pi.divide(ee.Image(180))).cos()  # in degrees
sindSunZe = SunZe.multiply(pi.divide(ee.Image(180))).sin()  # in degrees

# sat zenith
SatZe = ee.Image.constant(s2.get('MEAN_INCIDENCE_ZENITH_ANGLE_B5')).clip(footprint)
cosdSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).cos()
sindSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).sin()

# sat azimuth
SatAz = ee.Image.constant(s2.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B5')).clip(footprint)

# relative azimuth
RelAz = SatAz.subtract(SunAz)
cosdRelAz = RelAz.multiply(pi.divide(ee.Image(180))).cos()

# Pressure
P = ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM)).pow(5.25588)).multiply(0.01)
Po = ee.Image(1013.25)

# esun
ESUN = ee.Image(ee.Array([ee.Image(s2.get('SOLAR_IRRADIANCE_B1')),
                          ee.Image(s2.get('SOLAR_IRRADIANCE_B2')),
                          ee.Image(s2.get('SOLAR_IRRADIANCE_B3')),
                          ee.Image(s2.get('SOLAR_IRRADIANCE_B4')),
                          ee.Image(s2.get('SOLAR_IRRADIANCE_B5')),
                          ee.Image(s2.get('SOLAR_IRRADIANCE_B6')),
                          ee.Image(s2.get('SOLAR_IRRADIANCE_B7')),
                          ee.Image(s2.get('SOLAR_IRRADIANCE_B8')),
                          ee.Image(s2.get('SOLAR_IRRADIANCE_B8A')),
                          ee.Image(s2.get('SOLAR_IRRADIANCE_B11')),
                          ee.Image(s2.get('SOLAR_IRRADIANCE_B2'))]
                         )).toArray().toArray(1)

ESUN = ESUN.multiply(ee.Image(1))

ESUNImg = ESUN.arrayProject([0]).arrayFlatten([bands])

# create empty array for the images
imgArr = rescale.select(bands).toArray().toArray(1)

# pTOA to Ltoa
Ltoa = imgArr.multiply(ESUN).multiply(cosdSunZe).divide(pi.multiply(d.pow(2)))

# band centers
bandCenter = ee.Image(443).divide(1000).addBands(ee.Image(490).divide(1000))\
.addBands(ee.Image(560).divide(1000))\
.addBands(ee.Image(665).divide(1000))\
.addBands(ee.Image(705).divide(1000))\
.addBands(ee.Image(740).divide(1000))\
.addBands(ee.Number(783).divide(1000))\
.addBands(ee.Number(842).divide(1000))\
.addBands(ee.Number(865).divide(1000))\
.addBands(ee.Number(1610).divide(1000))\
.addBands(ee.Number(2190).divide(1000))\
.toArray().toArray(1)

# ozone coefficients
koz = ee.Image(0.0039).addBands(ee.Image(0.0213))\
.addBands(ee.Image(0.1052))\
.addBands(ee.Image(0.0505))\
.addBands(ee.Image(0.0205))\
.addBands(ee.Image(0.0112))\
.addBands(ee.Image(0.0075))\
.addBands(ee.Image(0.0021))\
.addBands(ee.Image(0.0019))\
.addBands(ee.Image(0))\
.addBands(ee.Image(0))\
.toArray().toArray(1)

###################

# Calculate ozone optical thickness
Toz = koz.multiply(DU).divide(ee.Image(1000))

# Calculate TOA radiance in the absense of ozone
Lt = Ltoa.multiply(((Toz)).multiply((ee.Image(1).divide(cosdSunZe)).add(ee.Image(1).divide(cosdSatZe))).exp())

# Rayleigh optical thickness
Tr = (P.divide(Po)).multiply(ee.Image(0.008569).multiply(bandCenter.pow(-4))).multiply((ee.Image(1).add(
    ee.Image(0.0113).multiply(bandCenter.pow(-2))).add(ee.Image(0.00013).multiply(bandCenter.pow(-4)))))

# Specular reflection (s- and p- polarization states)
theta_V = ee.Image(0.0000000001)
sin_theta_j = sindSunZe.divide(ee.Image(1.333))

theta_j = sin_theta_j.asin().multiply(ee.Image(180).divide(pi))

theta_SZ = SunZe

R_theta_SZ_s = (
    ((theta_SZ.multiply(pi.divide(ee.Image(180)))).subtract(theta_j.multiply(pi.divide(ee.Image(180))))).sin().pow(
        2)).divide(
    (((theta_SZ.multiply(pi.divide(ee.Image(180)))).add(theta_j.multiply(pi.divide(ee.Image(180))))).sin().pow(2)))

R_theta_V_s = ee.Image(0.0000000001)

R_theta_SZ_p = (((theta_SZ.multiply(pi.divide(180))).subtract(theta_j.multiply(pi.divide(180)))).tan().pow(2)).divide(
    (((theta_SZ.multiply(pi.divide(180))).add(theta_j.multiply(pi.divide(180)))).tan().pow(2)))

R_theta_V_p = ee.Image(0.0000000001)

R_theta_SZ = ee.Image(0.5).multiply(R_theta_SZ_s.add(R_theta_SZ_p))

R_theta_V = ee.Image(0.5).multiply(R_theta_V_s.add(R_theta_V_p))

# Sun-sensor geometry
theta_neg = ((cosdSunZe.multiply(ee.Image(-1))).multiply(cosdSatZe)).subtract(
    (sindSunZe).multiply(sindSatZe).multiply(cosdRelAz))

theta_neg_inv = theta_neg.acos().multiply(ee.Image(180).divide(pi))

theta_pos = (cosdSunZe.multiply(cosdSatZe)).subtract(sindSunZe.multiply(sindSatZe).multiply(cosdRelAz))

theta_pos_inv = theta_pos.acos().multiply(ee.Image(180).divide(pi))

cosd_tni = theta_neg_inv.multiply(pi.divide(180)).cos()  # in degrees

cosd_tpi = theta_pos_inv.multiply(pi.divide(180)).cos()  # in degrees

Pr_neg = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tni.pow(2))))

Pr_pos = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tpi.pow(2))))

# Rayleigh scattering phase function
Pr = Pr_neg.add((R_theta_SZ.add(R_theta_V)).multiply(Pr_pos))  # for water Rayleigh correction
# Pr = ee.Image(1) # for terrestrial Rayleigh correction

# rayleigh radiance contribution
denom = ee.Image(4).multiply(pi).multiply(cosdSatZe)
Lr = (ESUN.multiply(Tr)).multiply(Pr.divide(denom))

# rayleigh corrected radiance
Lrc = Lt.subtract(Lr)
LrcImg = Lrc.arrayProject([0]).arrayFlatten([bands])

# rayleigh corrected reflectance
prc = (Lrc.multiply(pi).multiply(d.pow(2)).divide(ESUN.multiply(cosdSunZe)))
prcImg = prc.arrayProject([0]).arrayFlatten([bands])


########################
# models #

# Calculate FAI
NIRprime = (prcImg.select('B4')).add((prcImg.select('B11').subtract(prcImg.select('B4'))).multiply(
    (ee.Image(865).subtract(ee.Image(665))).divide((ee.Image(1610).subtract(ee.Image(665))))))
FAI = ((prcImg.select('B8A').subtract(NIRprime))).multiply(mask)

prcLyr = prcImg.multiply(mask)
prc_vis = {'bands': 'B4,B3,B2', 'min': 0, 'max': 0.4}
prc_mapid = prcLyr.getMapId(prc_vis)

FAI_vis = {'min': -0.05, 'max': 0.2, 'palette': '#000080,#0080FF,#7BFF7B,#FF9700,#800000'}
FAI_mapid = FAI.getMapId(FAI_vis)

print 'PRC',prc_mapid
print 'FAI',FAI_mapid
# Map Layers
# Map.addLayer(footprint, {}, 'footprint', true)
# Map.addLayer(prcImg.multiply(mask), {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.4}, 'prc rgb', false)
# Map.addLayer(FAI, {min: -0.05, max: 0.2, palette: ['000080', '0080FF', '7BFF7B', 'FF9700', '800000']}, 'FAI', true)
#
# # export image #
#
# Export.image.toDrive({
#     image: FAI,
#     description: 's2_FAI',
#     scale: 30,
#     region: footprint
# })

