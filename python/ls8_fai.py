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

geometry = ee.Geometry.Polygon([[[30.76171875, 0.9049611504960419],
          [30.8935546875, -3.487377195492663],
          [35.5517578125, -3.2680324702882952],
          [35.5517578125, 1.9593043032313748]]])
l8 = ee.Image('LANDSAT/LC08/C01/T1/LC08_170060_20171226')

#Start date
iniDate = '2015-05-01'

#End Date
endDate = '2018-03-31'

oliCloudPerc = 5

#Landsat  8 raw dn
OLI_DN = ee.ImageCollection('LANDSAT/LC08/C01/T1')

#Landsat-8 surface reflactance product (for masking purposes)
SRP = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')

#toms / omi
ozone = ee.ImageCollection('TOMS/MERGED')

#Filtering Collection and Masking
pi = ee.Image(3.141592)

#Water Mask
startMonth = 5
endMonth = 9
startYear = 2013
endYear = 2017

forMask = SRP.filterBounds(geometry).select('B6').filterMetadata('CLOUD_COVER', "less_than", 10).filter(ee.Filter.calendarRange(startMonth, endMonth, 'month')).filter(ee.Filter.calendarRange(startYear, endYear, 'year'))
mask = ee.Image(forMask.select('B6').median().lt(600))
mask = mask.updateMask(mask)
FC_OLI = OLI_DN.filterDate(iniDate, endDate).filterBounds(geometry).filterMetadata('CLOUD_COVER', "less_than", oliCloudPerc)

#OLI image date
oliDate = l8.date()
footprint = l8.geometry()

#DEM
DEM_OLI = ee.Image('USGS/SRTMGL1_003').clip(footprint)

#Ozone
DU_OLI = ee.Image(ozone.filterDate(iniDate,endDate).filterBounds(footprint).mean())

#Julian Day
imgDate_OLI = ee.Date(l8.get('system:time_start'))
FOY_OLI = ee.Date.fromYMD(imgDate_OLI.get('year'),1,1)
JD_OLI = imgDate_OLI.difference(FOY_OLI,'day').int().add(1)

#Earth-Sun distance
d_OLI = ee.Image.constant(l8.get('EARTH_SUN_DISTANCE'))

#Sun elevation
SunEl_OLI = ee.Image.constant(l8.get('SUN_ELEVATION'))

#Sun azimuth
SunAz_OLI = ee.Image.constant(l8.get('SUN_AZIMUTH'))

#Satellite Zenith
SatZe_OLI = ee.Image(0.0)
cosdSatZe_OLI = (SatZe_OLI).multiply(pi.divide(ee.Image(180))).cos()
sindSatZe_OLI = (SatZe_OLI).multiply(pi.divide(ee.Image(180))).sin()

#Satellite Azimuth
SatAz_OLI = ee.Image(0.0)

#Sun zenith
SunZe_OLI = ee.Image(90).subtract(SunEl_OLI)
cosdSunZe_OLI = SunZe_OLI.multiply(pi.divide(ee.Image.constant(180))).cos()
sindSunZe_OLI = SunZe_OLI.multiply(pi.divide(ee.Image(180))).sin()

#Relative azimuth
RelAz_OLI = ee.Image(SunAz_OLI)
cosdRelAz_OLI = RelAz_OLI.multiply(pi.divide(ee.Image(180))).cos()

#Pressure calculation
P_OLI = ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM_OLI)).pow(5.25588)).multiply(0.01)
Po_OLI = ee.Image(1013.25)

#Radiometric Calibration
#define bands to be converted to radiance
bands_OLI = ['B1','B2','B3','B4','B5','B6','B7']

#Radiance Mult Bands
rad_mult_OLI = ee.Image(ee.Array([ee.Image(l8.get('RADIANCE_MULT_BAND_1')),
                        ee.Image(l8.get('RADIANCE_MULT_BAND_2')),
                        ee.Image(l8.get('RADIANCE_MULT_BAND_3')),
                        ee.Image(l8.get('RADIANCE_MULT_BAND_4')),
                        ee.Image(l8.get('RADIANCE_MULT_BAND_5')),
                        ee.Image(l8.get('RADIANCE_MULT_BAND_6')),
                        ee.Image(l8.get('RADIANCE_MULT_BAND_7'))]
                        )).toArray(1)

#Radiance add band
rad_add_OLI = ee.Image(ee.Array([ee.Image(l8.get('RADIANCE_ADD_BAND_1')),
                        ee.Image(l8.get('RADIANCE_ADD_BAND_2')),
                        ee.Image(l8.get('RADIANCE_ADD_BAND_3')),
                        ee.Image(l8.get('RADIANCE_ADD_BAND_4')),
                        ee.Image(l8.get('RADIANCE_ADD_BAND_5')),
                        ee.Image(l8.get('RADIANCE_ADD_BAND_6')),
                        ee.Image(l8.get('RADIANCE_ADD_BAND_7'))]
                        )).toArray(1)

# Create an empty image to save new radiance bands to
imgArr_OLI = l8.select(bands_OLI).toArray().toArray(1)
Ltoa_OLI = imgArr_OLI.multiply(rad_mult_OLI).add(rad_add_OLI)

#esun
ESUN_OLI = ee.Image.constant(197.24790954589844) \
                .addBands(ee.Image.constant(201.98426818847656))\
                .addBands(ee.Image.constant(186.12677001953125))\
                .addBands(ee.Image.constant(156.95257568359375))\
                .addBands(ee.Image.constant(96.04714965820312))\
                .addBands(ee.Image.constant(23.8833221450863))\
                .addBands(ee.Image.constant(8.04995873449635)).toArray().toArray(1)

ESUN_OLI = ESUN_OLI.multiply(ee.Image(1))

ESUNImg_OLI = ESUN_OLI.arrayProject([0]).arrayFlatten([bands_OLI])

#Ozone correction
#Ozone coefficients
koz_OLI = ee.Image.constant(0.0039).addBands(ee.Image.constant(0.0218))\
                          .addBands(ee.Image.constant(0.1078))\
                          .addBands(ee.Image.constant(0.0608))\
                          .addBands(ee.Image.constant(0.0019))\
                          .addBands(ee.Image.constant(0))\
                          .addBands(ee.Image.constant(0))\
                          .toArray().toArray(1)

#Calculate ozone optical thickness
Toz_OLI = koz_OLI.multiply(DU_OLI).divide(ee.Image.constant(1000))
# Calculate TOA radiance in the absense of ozone
Lt_OLI = Ltoa_OLI.multiply(((Toz_OLI)).multiply((ee.Image.constant(1).divide(cosdSunZe_OLI)).add(ee.Image.constant(1).divide(cosdSatZe_OLI))).exp())

# Rayleigh optical thickness
bandCenter_OLI = ee.Image(443).divide(1000).addBands(ee.Image(483).divide(1000))\
                                          .addBands(ee.Image(561).divide(1000))\
                                          .addBands(ee.Image(655).divide(1000))\
                                          .addBands(ee.Image(865).divide(1000))\
                                          .addBands(ee.Image(1609).divide(1000))\
                                          .addBands(ee.Number(2201).divide(1000))\
                                          .toArray().toArray(1)

# create an empty image to save new Tr values to
Tr_OLI = (P_OLI.divide(Po_OLI)).multiply(ee.Image(0.008569).multiply(bandCenter_OLI.pow(-4))).multiply((ee.Image(1).add(ee.Image(0.0113).multiply(bandCenter_OLI.pow(-2))).add(ee.Image(0.00013).multiply(bandCenter_OLI.pow(-4)))))

# Fresnel Reflection #
# Specular reflection (s- and p- polarization states)
theta_V_OLI = ee.Image(0.0000000001)
sin_theta_j_OLI = sindSunZe_OLI.divide(ee.Image(1.333))

theta_j_OLI = sin_theta_j_OLI.asin().multiply(ee.Image(180).divide(pi))

theta_SZ_OLI = SunZe_OLI

R_theta_SZ_s_OLI = (((theta_SZ_OLI.multiply(pi.divide(ee.Image(180)))).subtract(theta_j_OLI.multiply(pi.divide(ee.Image(180))))).sin().pow(2)).divide((((theta_SZ_OLI.multiply(pi.divide(ee.Image(180)))).add(theta_j_OLI.multiply(pi.divide(ee.Image(180))))).sin().pow(2)))

R_theta_V_s_OLI = ee.Image(0.0000000001)
R_theta_SZ_p_OLI = (((theta_SZ_OLI.multiply(pi.divide(180))).subtract(theta_j_OLI.multiply(pi.divide(180)))).tan().pow(2)).divide((((theta_SZ_OLI.multiply(pi.divide(180))).add(theta_j_OLI.multiply(pi.divide(180)))).tan().pow(2)))
R_theta_V_p_OLI = ee.Image(0.0000000001)
R_theta_SZ_OLI = ee.Image(0.5).multiply(R_theta_SZ_s_OLI.add(R_theta_SZ_p_OLI))
R_theta_V_OLI = ee.Image(0.5).multiply(R_theta_V_s_OLI.add(R_theta_V_p_OLI))

# Rayleigh scattering phase function #
# Sun-sensor geometry
theta_neg_OLI = ((cosdSunZe_OLI.multiply(ee.Image(-1))).multiply(cosdSatZe_OLI)).subtract((sindSunZe_OLI).multiply(sindSatZe_OLI).multiply(cosdRelAz_OLI))

theta_neg_inv_OLI = theta_neg_OLI.acos().multiply(ee.Image(180).divide(pi))

theta_pos_OLI = (cosdSunZe_OLI.multiply(cosdSatZe_OLI)).subtract(sindSunZe_OLI.multiply(sindSatZe_OLI).multiply(cosdRelAz_OLI))

theta_pos_inv_OLI = theta_pos_OLI.acos().multiply(ee.Image(180).divide(pi))

cosd_tni_OLI = theta_neg_inv_OLI.multiply(pi.divide(180)).cos() # in degrees

cosd_tpi_OLI = theta_pos_inv_OLI.multiply(pi.divide(180)).cos() # in degrees

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
pc = prc_OLI.arrayProject([0]).arrayFlatten([bands_OLI])

pcVisParams = {'bands': 'B4,B3,B2','min': 0,'max': 0.1}
pc = pc.multiply(mask)
pcImgID = pc.getMapId(pcVisParams)
# Calculate FAI
NIRprime = (pc.select('B4')).add((pc.select('B6').subtract(pc.select('B4'))).multiply((ee.Image(865).subtract(ee.Image(655))).divide((ee.Image(1609).subtract(ee.Image(655))))))
fai = (pc.select('B5').subtract(NIRprime))
fai = fai.multiply(mask)

faiVisParams = {'min': -0.05, 'max': 0.2, 'palette': '#000080,#0080FF,#7BFF7B,#FF9700,#800000'}
faiImgID = fai.getMapId(faiVisParams)

print pcImgID
print faiImgID












































