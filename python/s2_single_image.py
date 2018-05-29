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

# This script processes a single Sentinel-2 image of interest.
# WQ parameters such chlor-a, secchi depth, trophic state index are calculated
# How to use:
# (1) If a "geometry" iable exists in the imports window, delete it.
# (2) Within the map, select the point button near the top left side.
# (3) Create a new point by clicking on a location of interest.
# (4) Adjust the time frame you wish to browse by adding here:
# begin date
iniDate = '2015-05-01'
# end date
endDate = '2018-03-31'
# (5) Adjust a cloud % threshold here:
cloudPerc = 5
# (6) Click Run
# (7) The "available imagery" ImageCollection within the console displays all available imagery. If you wish to investigate a single
#     image from this collection, highlight the FILE_ID in the available imagery features and copy and paste it here:
s2 = ee.Image('COPERNICUS/S2/20171128T075251_20171128T080644_T36MXE')
# (8) The charts will appear on the right hand side of the interface, and will display the mean map in the map user interface.
# (9) Click "Run"
# (10) Export each image by clicking on the run button within the "tasks" tab.

# Author: Benjamin Page #
# Citations:
# Page, B.P. and Mishra, D.R., 2018. A modified atmospheric correction when coupling sentinel-2 and landsat-8 for inland water quality monitoring

#########################################################/
#########################################################/
#########################################################/
# Import Collections #

# sentinel-2
MSI = ee.ImageCollection('COPERNICUS/S2')

# landsat-8 surface reflactance product (for masking purposes)
SRP = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')

# toms / omi
ozone = ee.ImageCollection('TOMS/MERGED')

#########################################################/
#########################################################/
#########################################################/

pi = ee.Image(3.141592)
imDate = s2.date()

footprint = s2.geometry()

# water mask
startMonth = 5
endMonth = 9
startYear = 2013
endYear = 2017

forMask = SRP.filterBounds(geometry).select('B6').filterMetadata('CLOUD_COVER', "less_than", 10).filter(
    ee.Filter.calendarRange(startMonth, endMonth, 'month')).filter(ee.Filter.calendarRange(startYear, endYear, 'year'))
mask = ee.Image(forMask.select('B6').median().lt(300))
mask = mask.updateMask(mask).clip(footprint)

# filter sentinel 2 collection
FC = MSI.filterDate(iniDate, endDate).filterBounds(geometry).filterMetadata('CLOUDY_PIXEL_PERCENTAGE', "less_than",
                                                                            cloudPerc)


#########################################################/
#########################################################/
#########################################################/
# MSI Atmospheric Correction #

bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']

# rescale*
rescale = ee.Image(s2.divide(10000).multiply(mask).copyProperties(s2)).select(bands)

# dem*
DEM = ee.Image('USGS/SRTMGL1_003').clip(footprint)

# ozone*
DU = ee.Image(ozone.filterDate(iniDate, endDate).filterBounds(footprint).mean())

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
P = (ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM)).pow(5.25588)).multiply(
    0.01)).multiply(mask)
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
Pr = Pr_neg.add((R_theta_SZ.add(R_theta_V)).multiply(Pr_pos))

# rayleigh radiance contribution
denom = ee.Image(4).multiply(pi).multiply(cosdSatZe)
Lr = (ESUN.multiply(Tr)).multiply(Pr.divide(denom))

# rayleigh corrected radiance
Lrc = Lt.subtract(Lr)
LrcImg = Lrc.arrayProject([0]).arrayFlatten([bands])

# Aerosol Correction #

# Bands in nm
bands_nm = ee.Image(443).addBands(ee.Image(490))\
.addBands(ee.Image(560))\
.addBands(ee.Image(665))\
.addBands(ee.Image(705))\
.addBands(ee.Image(740))\
.addBands(ee.Image(783))\
.addBands(ee.Image(842))\
.addBands(ee.Image(865))\
.addBands(ee.Image(0))\
.addBands(ee.Image(0))\
.toArray().toArray(1)

# Lam in SWIR bands
Lam_10 = LrcImg.select('B11')
Lam_11 = LrcImg.select('B12')

# Calculate aerosol type
eps = (
    (((Lam_11).divide(ESUNImg.select('B12'))).log()).subtract(((Lam_10).divide(ESUNImg.select('B11'))).log())).divide(
    ee.Image(2190).subtract(ee.Image(1610)))

# Calculate multiple scattering of aerosols for each band
Lam = (Lam_11).multiply(((ESUN).divide(ESUNImg.select('B12')))).multiply(
    (eps.multiply(ee.Image(-1))).multiply((bands_nm.divide(ee.Image(2190)))).exp())

# diffuse transmittance
trans = Tr.multiply(ee.Image(-1)).divide(ee.Image(2)).multiply(ee.Image(1).divide(cosdSatZe)).exp()

# Compute water-leaving radiance
Lw = Lrc.subtract(Lam).divide(trans)

# water-leaving reflectance
pw = (Lw.multiply(pi).multiply(d.pow(2)).divide(ESUN.multiply(cosdSunZe)))
pwImg = pw.arrayProject([0]).arrayFlatten([bands])

# remote sensing reflectance
Rrs = (pw.divide(pi).arrayProject([0]).arrayFlatten([bands]).slice(0, 9))


# / Bio optical Models #/
# chlor_a
NDCI = (Rrs.select('B5').subtract(Rrs.select('B4'))).divide(Rrs.select('B5').add(Rrs.select('B4')))  # ndci
chlor_a = ee.Image(14.039).add(ee.Image(86.115).multiply(NDCI)).add(
    ee.Image(194.325).multiply(NDCI.pow(ee.Image(2))))  # chlor_a

# SD
ln_BlueRed = (Rrs.select('B2').divide(Rrs.select('B4'))).log()
lnMOSD = (ee.Image(1.4856).multiply(ln_BlueRed)).add(ee.Image(0.2734))  # R2 = 0.8748 with in-situ
MOSD = ee.Image(10).pow(lnMOSD)  # log space to (m)
SD = (ee.Image(0.1777).multiply(MOSD)).add(ee.Image(1.0813))

# tsi
TSI_c = ee.Image(30.6).add(ee.Image(9.81)).multiply(chlor_a.log())
TSI_s = ee.Image(60.0).subtract(ee.Image(14.41)).multiply(SD.log())
TSI = (TSI_c.add(TSI_s)).divide(ee.Image(2))

# tsi reclassified
# Create conditions
mask1 = TSI.lt(30)  # classical oligitrophy (1)
mask2 = TSI.gte(30).And(TSI.lt(40))  # (2)
mask3 = TSI.gte(40).And(TSI.lt(50))  # (3)
mask4 = TSI.gte(50).And(TSI.lt(60))  # (4)
mask5 = TSI.gte(60).And(TSI.lt(70))  # (5)
mask6 = TSI.gte(70).And(TSI.lt(80))  # (6)
mask7 = TSI.gte(80)  # (7)

# Reclassify conditions into new values
img1 = TSI.where(mask1.eq(1), 1).mask(mask1)
img2 = TSI.where(mask2.eq(1), 2).mask(mask2)
img3 = TSI.where(mask3.eq(1), 3).mask(mask3)
img4 = TSI.where(mask4.eq(1), 4).mask(mask4)
img5 = TSI.where(mask5.eq(1), 5).mask(mask5)
img6 = TSI.where(mask6.eq(1), 6).mask(mask6)
img7 = TSI.where(mask7.eq(1), 7).mask(mask7)

# Ouput of reclassified image
TSI_R = img1.unmask(img2).unmask(img3).unmask(img4).unmask(img5).unmask(img6).unmask(img7)

s2Lyr = s2.multiply(mask)
s2_vis = {'bands': 'B4,B3,B2', 'min': 0, 'max': 1000}
s2_mapid = s2Lyr.getMapId(s2_vis)

chlor_a_vis = {'min': 0, 'max': 40, 'palette': 'blue,cyan,limegreen,yellow,orange,darkred'}
chlor_a_mapid = chlor_a.getMapId(chlor_a_vis)

TSI_vis = {'min': 30, 'max': 70, 'palette': 'blue,cyan,limegreen,yellow,orange,darkred'}
TSI_mapid = TSI.getMapId(TSI_vis)

TSI_R_vis = {'min': 1, 'max': 7, 'palette': 'purple,blue,limegreen,yellow,orange,orangered,darkred'}
TSI_R_mapid = TSI_R.getMapId(TSI_R_vis)

print 'S2',s2_mapid
print 'Chlor A',chlor_a_mapid
print 'TSI',TSI_mapid
print 'TSI R',TSI_R_mapid
# Map Layers #
# Map.centerObject(footprint, 7)
# Map.addLayer(footprint, {}, 'footprint', true)  # footprint
# Map.addLayer(s2.multiply(mask), {bands: ['B4', 'B3', 'B2'], min: 0, max: 1000}, 'rgb', false)  # rgb
# Map.addLayer(SD, {min: 0, max: 3, palette: ['darkred', 'orange', 'yellow', 'limegreen', 'cyan', 'blue']}, 'SD', false)
# Map.addLayer(chlor_a, {min: 0, max: 40, palette: ['blue', 'cyan', 'limegreen', 'yellow', 'orange', 'darkred']},
#              'chlor_a', false)
# Map.addLayer(TSI, {min: 30, max: 70, palette: ['blue', 'cyan', 'limegreen', 'yellow', 'orange', 'darkred']}, 'TSI',
#              false)
# Map.addLayer(TSI_R,
#              {min: 1, max: 7, palette: ['purple', 'blue', 'limegreen', 'yellow', 'orange', 'orangered', 'darkred']},
#              'TSI_R', true)

# EXPORT IMAGES #

# # Export SD
# Export.image.toDrive({
#     image: SD,
#     description: 'SD',
#     scale: 20,
#     region: footprint
# })
#
# # Export chlor_a
# Export.image.toDrive({
#     image: chlor_a,
#     description: 'chlor_a',
#     scale: 20,
#     region: footprint
# })
#
# # Export TSI
# Export.image.toDrive({
#     image: TSI,
#     description: 'TSI',
#     scale: 20,
#     region: footprint
# })
#
# # Export TSI_R
# Export.image.toDrive({
#     image: TSI_R,
#     description: 'TSI_R',
#     scale: 20,
#     region: footprint
# })