import ee
from ee.ee_exception import EEException

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
# This script processes a single Landsat 8 image of interest.
# WQ parameters such chlor-a, secchi depth, trophic state index are calculated
# How to use:
# (1) If a "geometry" variable exists in the imports window, delete it.
# (2) Within the map, select the point button near the top left side.
# (3) Create a new point by clicking on a location of interest.
# (4) Adjust the time frame you wish to browse by adding here:

geometry = ee.Geometry.Polygon([[[30.76171875, 0.9049611504960419],
          [30.8935546875, -3.487377195492663],
          [35.5517578125, -3.2680324702882952],
          [35.5517578125, 1.9593043032313748]]])
# begin date
iniDate = '2015-05-01'
# end date
endDate = '2018-03-31'
# (5) Adjust a cloud % threshold here:
oliCloudPerc = 5
# (6) Click Run
# (7) The "available imagery" ImageCollection within the console displays all available imagery. If you wish to investigate a single
#     image from this collection, find the FILE_ID within the features and copy and paste it here:
l8 = ee.Image('LANDSAT/LC08/C01/T1/LC08_170060_20171226')
# (8) Click "Run"
# (9) Export each image by clicking on the run button within the "tasks" tab.

# Author: Benjamin Page #
# Citations:
# Page, B.P. and Mishra, D.R., 2018. A modified atmospheric correction when coupling sentinel-2 and landsat-8 for inland water quality monitoring, In Review

#########################################################/
#########################################################/
#########################################################/

# Import Collections #

# landsat 8 raw dn
OLI_DN = ee.ImageCollection('LANDSAT/LC08/C01/T1')

# landsat-8 surface reflactance product (for masking purposes)
SRP = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')

# toms / omi
ozone = ee.ImageCollection('TOMS/MERGED')

#########################################################/
#########################################################/
#########################################################/
# Filtering Collection and Masking #

pi = ee.Image(3.141592)

# water mask
startMonth = 5
endMonth = 9
startYear = 2013
endYear = 2017

forMask = SRP.filterBounds(geometry).select('B6').filterMetadata('CLOUD_COVER', "less_than", 10).filter \
    (ee.Filter.calendarRange(startMonth, endMonth, 'month')).filter(ee.Filter.calendarRange(startYear, endYear, 'year'))
mask = ee.Image(forMask.select('B6').median().lt(600))
mask = mask.updateMask(mask)

# filter landsat 8 collection
FC_OLI = OLI_DN.filterDate(iniDate, endDate).filterBounds(geometry).filterMetadata('CLOUD_COVER', "less_than", oliCloudPerc)

# oli image date
oliDate = l8.date()

footprint = l8.geometry()

#########################################################/
#########################################################/
#########################################################/
################################################################################
# single OLI image ATMOSPHERIC CORRECTION #

# dem
DEM_OLI = ee.Image('USGS/SRTMGL1_003').clip(footprint)

# ozone
DU_OLI = ee.Image(ozone.filterDate(iniDate ,endDate).filterBounds(footprint).mean())

# Julian Day
imgDate_OLI = ee.Date(l8.get('system:time_start'))
FOY_OLI = ee.Date.fromYMD(imgDate_OLI.get('year') ,1 ,1)
JD_OLI = imgDate_OLI.difference(FOY_OLI ,'day').int().add(1)

# Earth-Sun distance
d_OLI = ee.Image.constant(l8.get('EARTH_SUN_DISTANCE'))

# Sun elevation
SunEl_OLI = ee.Image.constant(l8.get('SUN_ELEVATION'))

# Sun azimuth
SunAz_OLI = ee.Image.constant(l8.get('SUN_AZIMUTH'))

# Satellite zenith
SatZe_OLI = ee.Image(0.0)
cosdSatZe_OLI = (SatZe_OLI).multiply(pi.divide(ee.Image(180))).cos()
sindSatZe_OLI = (SatZe_OLI).multiply(pi.divide(ee.Image(180))).sin()

# Satellite azimuth
SatAz_OLI = ee.Image(0.0)

# Sun zenith
SunZe_OLI = ee.Image(90).subtract(SunEl_OLI)
cosdSunZe_OLI = SunZe_OLI.multiply(pi.divide(ee.Image.constant(180))).cos() # in degrees
sindSunZe_OLI = SunZe_OLI.multiply(pi.divide(ee.Image(180))).sin() # in degrees

# Relative azimuth
RelAz_OLI = ee.Image(SunAz_OLI)
cosdRelAz_OLI = RelAz_OLI.multiply(pi.divide(ee.Image(180))).cos()

# Pressure calculation
P_OLI = ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM_OLI)).pow(5.25588)).multiply \
    (0.01)
Po_OLI = ee.Image(1013.25)

# Radiometric Calibration #
# define bands to be converted to radiance
bands_OLI = ['B1' ,'B2' ,'B3' ,'B4' ,'B5' ,'B6' ,'B7']

# radiance_mult_bands
rad_mult_OLI = ee.Image(ee.Array([ee.Image(l8.get('RADIANCE_MULT_BAND_1')),
                                  ee.Image(l8.get('RADIANCE_MULT_BAND_2')),
                                  ee.Image(l8.get('RADIANCE_MULT_BAND_3')),
                                  ee.Image(l8.get('RADIANCE_MULT_BAND_4')),
                                  ee.Image(l8.get('RADIANCE_MULT_BAND_5')),
                                  ee.Image(l8.get('RADIANCE_MULT_BAND_6')),
                                  ee.Image(l8.get('RADIANCE_MULT_BAND_7'))]
                                 )).toArray(1)

# radiance add band
rad_add_OLI = ee.Image(ee.Array([ee.Image(l8.get('RADIANCE_ADD_BAND_1')),
                                 ee.Image(l8.get('RADIANCE_ADD_BAND_2')),
                                 ee.Image(l8.get('RADIANCE_ADD_BAND_3')),
                                 ee.Image(l8.get('RADIANCE_ADD_BAND_4')),
                                 ee.Image(l8.get('RADIANCE_ADD_BAND_5')),
                                 ee.Image(l8.get('RADIANCE_ADD_BAND_6')),
                                 ee.Image(l8.get('RADIANCE_ADD_BAND_7'))]
                                )).toArray(1)

# create an empty image to save new radiance bands to
imgArr_OLI = l8.select(bands_OLI).toArray().toArray(1)
Ltoa_OLI = imgArr_OLI.multiply(rad_mult_OLI).add(rad_add_OLI)

# esun
ESUN_OLI = ee.Image.constant(197.24790954589844)\
.addBands(ee.Image.constant(201.98426818847656))\
.addBands(ee.Image.constant(186.12677001953125))\
.addBands(ee.Image.constant(156.95257568359375))\
.addBands(ee.Image.constant(96.04714965820312))\
.addBands(ee.Image.constant(23.8833221450863))\
.addBands(ee.Image.constant(8.04995873449635)).toArray().toArray(1)
ESUN_OLI = ESUN_OLI.multiply(ee.Image(1))

ESUNImg_OLI = ESUN_OLI.arrayProject([0]).arrayFlatten([bands_OLI])

# Ozone Correction #
# Ozone coefficients
koz_OLI = ee.Image.constant(0.0039).addBands(ee.Image.constant(0.0218))\
.addBands(ee.Image.constant(0.1078))\
.addBands(ee.Image.constant(0.0608))\
.addBands(ee.Image.constant(0.0019))\
.addBands(ee.Image.constant(0))\
.addBands(ee.Image.constant(0))\
.toArray().toArray(1)

# Calculate ozone optical thickness
Toz_OLI = koz_OLI.multiply(DU_OLI).divide(ee.Image.constant(1000))

# Calculate TOA radiance in the absense of ozone
Lt_OLI = Ltoa_OLI.multiply(((Toz_OLI)).multiply
    ((ee.Image.constant(1).divide(cosdSunZe_OLI)).add(ee.Image.constant(1).divide(cosdSatZe_OLI))).exp())

# Rayleigh optical thickness
bandCenter_OLI = ee.Image(443).divide(1000).addBands(ee.Image(483).divide(1000))\
.addBands(ee.Image(561).divide(1000))\
.addBands(ee.Image(655).divide(1000))\
.addBands(ee.Image(865).divide(1000))\
.addBands(ee.Image(1609).divide(1000))\
.addBands(ee.Number(2201).divide(1000))\
.toArray().toArray(1)

# create an empty image to save new Tr values to
Tr_OLI = (P_OLI.divide(Po_OLI)).multiply(ee.Image(0.008569).multiply(bandCenter_OLI.pow(-4))).multiply((ee.Image(1).add\
        (ee.Image(0.0113).multiply(bandCenter_OLI.pow(-2))).add(ee.Image(0.00013).multiply(bandCenter_OLI.pow(-4)))))

# Fresnel Reflection #
# Specular reflection (s- and p- polarization states)
theta_V_OLI = ee.Image(0.0000000001)
sin_theta_j_OLI = sindSunZe_OLI.divide(ee.Image(1.333))

theta_j_OLI = sin_theta_j_OLI.asin().multiply(ee.Image(180).divide(pi))

theta_SZ_OLI = SunZe_OLI

R_theta_SZ_s_OLI = (((theta_SZ_OLI.multiply(pi.divide(ee.Image(180)))).subtract
    (theta_j_OLI.multiply(pi.divide(ee.Image(180))))).sin().pow(2)).divide((((theta_SZ_OLI.multiply
    (pi.divide(ee.Image(180)))).add(theta_j_OLI.multiply(pi.divide(ee.Image(180))))).sin().pow(2)))

R_theta_V_s_OLI = ee.Image(0.0000000001)

R_theta_SZ_p_OLI = (((theta_SZ_OLI.multiply(pi.divide(180))).subtract(theta_j_OLI.multiply(pi.divide(180)))).tan().pow(2)).\
    divide((((theta_SZ_OLI.multiply(pi.divide(180))).add(theta_j_OLI.multiply(pi.divide(180)))).tan().pow(2)))

R_theta_V_p_OLI = ee.Image(0.0000000001)

R_theta_SZ_OLI = ee.Image(0.5).multiply(R_theta_SZ_s_OLI.add(R_theta_SZ_p_OLI))

R_theta_V_OLI = ee.Image(0.5).multiply(R_theta_V_s_OLI.add(R_theta_V_p_OLI))

# Rayleigh scattering phase function #
# Sun-sensor geometry
theta_neg_OLI = ((cosdSunZe_OLI.multiply(ee.Image(-1))).multiply(cosdSatZe_OLI)).\
    subtract((sindSunZe_OLI).multiply(sindSatZe_OLI).multiply(cosdRelAz_OLI))

theta_neg_inv_OLI = theta_neg_OLI.acos().multiply(ee.Image(180).divide(pi))

theta_pos_OLI = (cosdSunZe_OLI.multiply(cosdSatZe_OLI)).\
    subtract(sindSunZe_OLI.multiply(sindSatZe_OLI).multiply(cosdRelAz_OLI))

theta_pos_inv_OLI = theta_pos_OLI.acos().multiply(ee.Image(180).divide(pi))

cosd_tni_OLI = theta_neg_inv_OLI.multiply(pi.divide(180)).cos()  # in degrees

cosd_tpi_OLI = theta_pos_inv_OLI.multiply(pi.divide(180)).cos()  # in degrees

Pr_neg_OLI = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tni_OLI.pow(2))))

Pr_pos_OLI = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tpi_OLI.pow(2))))

# Rayleigh scattering phase function
Pr_OLI = Pr_neg_OLI.add((R_theta_SZ_OLI.add(R_theta_V_OLI)).multiply(Pr_pos_OLI))

# Calulate Lr,
denom_OLI = ee.Image(4).multiply(pi).multiply(cosdSatZe_OLI)
Lr_OLI = (ESUN_OLI.multiply(Tr_OLI)).multiply(Pr_OLI.divide(denom_OLI))

# Rayleigh corrected radiance
Lrc_OLI = (Lt_OLI.divide(ee.Image(10))).subtract(Lr_OLI)
LrcImg_OLI = Lrc_OLI.arrayProject([0]).arrayFlatten([bands_OLI])

# Rayleigh corrected reflectance
prc_OLI = Lrc_OLI.multiply(pi).multiply(d_OLI.pow(2)).divide(ESUN_OLI.multiply(cosdSunZe_OLI))
prcImg_OLI = prc_OLI.arrayProject([0]).arrayFlatten([bands_OLI])

# Aerosol Correction #
# Bands in nm
bands_nm_OLI = ee.Image(443).addBands(ee.Image(483))\
.addBands(ee.Image(561))\
.addBands(ee.Image(655))\
.addBands(ee.Image(865))\
.addBands(ee.Image(0))\
.addBands(ee.Image(0))\
.toArray().toArray(1)

# Lam in SWIR bands
Lam_6_OLI = LrcImg_OLI.select('B6')
Lam_7_OLI = LrcImg_OLI.select('B7')

# Calculate aerosol type
eps_OLI = (((((Lam_7_OLI).divide(ESUNImg_OLI.select('B7'))).log()).\
    subtract(((Lam_6_OLI).divide(ESUNImg_OLI.select('B6'))).log())).divide(ee.Image(2201).subtract(ee.Image(1609))))\
    .multiply(mask)

# Calculate multiple scattering of aerosols for each band
Lam_OLI = (Lam_7_OLI).multiply(((ESUN_OLI).divide(ESUNImg_OLI.select('B7'))))\
    .multiply((eps_OLI.multiply(ee.Image(-1))).multiply((bands_nm_OLI.divide(ee.Image(2201)))).exp())

# diffuse transmittance
trans_OLI = Tr_OLI.multiply(ee.Image(-1)).divide(ee.Image(2)).multiply(ee.Image(1).divide(cosdSatZe_OLI)).exp()

# Compute water-leaving radiance
Lw_OLI = Lrc_OLI.subtract(Lam_OLI).divide(trans_OLI)

# water-leaving reflectance
pw_OLI = (Lw_OLI.multiply(pi).multiply(d_OLI.pow(2)).divide(ESUN_OLI.multiply(cosdSunZe_OLI)))
pwImg_OLI = pw_OLI.arrayProject([0]).arrayFlatten([bands_OLI])

# Rrs
Rrs = (pw_OLI.divide(pi).arrayProject([0]).arrayFlatten([bands_OLI]).slice(0, 5))

###########################################################/
# Models #

# Surface water temperature
TIRS_1 = l8.select('B10')

b10_add_band = ee.Number(TIRS_1.get('RADIANCE_ADD_BAND_10'))
b10_mult_band = ee.Number(TIRS_1.get('RADIANCE_MULT_BAND_10'))

Oi = ee.Number(0.29)

TIRS_cal = (TIRS_1.multiply(b10_mult_band).add(b10_add_band).subtract(Oi))  # .multiply(lakes)

K1_b10 = ee.Number(TIRS_1.get('K1_CONSTANT_BAND_10'))
K2_b10 = ee.Number(TIRS_1.get('K2_CONSTANT_BAND_10'))

LST = ((ee.Image(K2_b10).divide(((ee.Image(K1_b10).divide(ee.Image(TIRS_cal))).add(ee.Image(1))).log()))\
    .subtract(ee.Image(273))).multiply(mask)  # Celsius
LST = (ee.Image(0.7745).multiply(LST)).add(ee.Image(9.6502))  # calibration R2 = 0.8599

# chlor_a
# Chlorophyll-a OC3
a0 = ee.Image(0.2412)
a1 = ee.Image(-2.0546)
a2 = ee.Image(1.1776)
a3 = ee.Image(-0.5538)
a4 = ee.Image(-0.4570)

log_BG = (Rrs.select('B1').divide(Rrs.select('B3'))).log10()

a1a = a1.multiply(log_BG.pow(1))
a2a = a2.multiply(log_BG.pow(2))
a3a = a3.multiply(log_BG.pow(3))
a4a = a4.multiply(log_BG.pow(4))

sum = a1a.add(a2a).add(a3a).add(a4a)

log10_chlor_a = a0.add(sum)

chlor_a = ee.Image(10).pow(log10_chlor_a)
chlor_a_cal = ee.Image(4.0752).multiply(chlor_a).subtract(ee.Image(3.9617))

# SD
ln_BlueRed = (Rrs.select('B2').divide(Rrs.select('B4'))).log()
lnMOSD = (ee.Image(1.4856).multiply(ln_BlueRed)).add(ee.Image(0.2734))  # R2 = 0.8748 with in-situ
MOSD = ee.Image(10).pow(lnMOSD)  # log space to (m)
SD = (ee.Image(0.1777).multiply(MOSD)).add(ee.Image(1.0813))

# tsi
TSI_c = ee.Image(30.6).add(ee.Image(9.81)).multiply(chlor_a_cal.log())
TSI_s = ee.Image(60.0).subtract(ee.Image(14.41)).multiply(SD.log())
TSI = (TSI_c.add(TSI_s)).divide(ee.Image(2))

# Reclassify TSI

# Create conditions
mask1 = TSI.lt(30)  # (1)
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

l8Lyr = l8.multiply(mask)
l8VisParams = {'bands': 'B4,B3,B2', 'min': 0, 'max': 15000}
l8MapId = l8Lyr.getMapId(l8VisParams)

LSTvisParams = {'min': 20, 'max': 30, 'palette': 'darkblue,blue,white,red,darkred'}
LSTMapId = LST.getMapId(LSTvisParams)

chlor_a_cal_visParams = {'min': 0, 'max': 80, 'palette':'blue,cyan,limegreen,yellow,orange,darkred'}
chlor_a_cal_mapid = chlor_a_cal.getMapId(chlor_a_cal_visParams)

SDvisParams = {'min': 0, 'max': 3, 'palette': 'darkred,orange,yellow,limegreen,cyan,blue'}
SDMapId = SD.getMapId(SDvisParams)

TSIvisParams = {'min': 30, 'max': 70, 'palette': 'blue,cyan,limegreen,yellow,orange,darkred'}
TSIMapId = TSI.getMapId(TSIvisParams)

TSI_R_visParams = {'min': 1, 'max': 7, 'palette': 'purple,blue,limegreen,yellow,orange,orangered,darkred'}
TSI_R_MapId = TSI_R.getMapId(TSI_R_visParams)

print 'L8',l8MapId
print 'LST',LSTMapId
print 'Chlor A',chlor_a_cal_mapid
print 'SD',SDMapId
print 'TSI',TSIMapId
print 'TSI R',TSI_R_MapId
# Map Layers

# Map.addLayer(l8.multiply(mask), {bands: ['B4', 'B3', 'B2'], min: 0, max: 15000}, 'rgb', false)  # rgb
# Map.addLayer(LST, {min: 20, max: 30, palette: ['darkblue', 'blue', 'white', 'red', 'darkred']}, 'LST', false)
# Map.addLayer(chlor_a_cal, {min: 0, max: 80, palette: ['blue', 'cyan', 'limegreen', 'yellow', 'orange', 'darkred']},
#              'chlor_a', false)
# Map.addLayer(SD, {min: 0, max: 3, palette: ['darkred', 'orange', 'yellow', 'limegreen', 'cyan', 'blue']}, 'SD', false)
# Map.addLayer(TSI, {min: 30, max: 70, palette: ['blue', 'cyan', 'limegreen', 'yellow', 'orange', 'darkred']}, 'TSI',
#              false)
# Map.addLayer(TSI_R,
#              {min: 1, max: 7, palette: ['purple', 'blue', 'limegreen', 'yellow', 'orange', 'orangered', 'darkred']},
#              'TSI_R', true)
#
#
