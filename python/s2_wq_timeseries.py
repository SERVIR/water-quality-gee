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

# This script processes and charts WQ parameter valuse from a collection of Sentinel-2 archive for a particular pixel throughout time.
# WQ parameters such chlor-a, secchi depth, trophic state index
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
# (5) Click Run
# (5) The "available imagery" ImageCollection within the console displays all available imagery. If you wish to investigate a single
#     image from this collection, highlight the FILE_ID and copy and paste it into the proper location within the "s2 Single Image" script.
# (6) The charts will appear on the right hand side of the interface, and will display the mean map in the map user interface.
# (6) Click "Run"
# (7) Export each chart as a csv by clicking on the export icon near the top-right of the chart.

# Author: Benjamin Page #
# Citations:
# Page, B.P. and Mishra, D.R., 2018. A modified atmospheric correction when coupling sentinel-2 and landsat-8 for inland water quality monitoring

#########################################################/
#########################################################/
#########################################################/
# Map.centerObject(geometry, 7)

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

# water mask
startMonth = 5
endMonth = 9
startYear = 2013
endYear = 2017

forMask = SRP.filterBounds(geometry).select('B6').filterMetadata('CLOUD_COVER', "less_than", 10).filter(
    ee.Filter.calendarRange(startMonth, endMonth, 'month')).filter(ee.Filter.calendarRange(startYear, endYear, 'year'))
mask = ee.Image(forMask.select('B6').median().lt(300))
mask = mask.updateMask(mask)

# filter sentinel 2 collection
FC = MSI.filterDate(iniDate, endDate).filterBounds(geometry).filterMetadata('CLOUDY_PIXEL_PERCENTAGE', "less_than",
                                                                            cloudPerc)


# #########################################################/
# #########################################################/
# #########################################################/
# Mapping functions #

def s2Correction(img):
    # msi bands
    bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']

    # rescale
    rescale = img.select(bands).divide(10000).multiply(mask)

    # tile footprint
    footprint = rescale.geometry()

    # dem
    DEM = ee.Image('USGS/SRTMGL1_003').clip(footprint)

    # ozone
    DU = ee.Image(ozone.filterDate(iniDate, endDate).filterBounds(footprint).mean())

    # Julian Day
    imgDate = ee.Date(img.get('system:time_start'))
    FOY = ee.Date.fromYMD(imgDate.get('year'), 1, 1)
    JD = imgDate.difference(FOY, 'day').int().add(1)

    # earth-sun distance
    myCos = ((ee.Image(0.0172).multiply(ee.Image(JD).subtract(ee.Image(2)))).cos()).pow(2)
    cosd = myCos.multiply(pi.divide(ee.Image(180))).cos()
    d = ee.Image(1).subtract(ee.Image(0.01673)).multiply(cosd).clip(footprint)

    # sun azimuth
    SunAz = ee.Image.constant(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')).clip(footprint)

    # sun zenith
    SunZe = ee.Image.constant(img.get('MEAN_SOLAR_ZENITH_ANGLE')).clip(footprint)
    cosdSunZe = SunZe.multiply(pi.divide(ee.Image(180))).cos()  # in degrees
    sindSunZe = SunZe.multiply(pi.divide(ee.Image(180))).sin()  # in degrees

    # sat zenith
    SatZe = ee.Image.constant(img.get('MEAN_INCIDENCE_ZENITH_ANGLE_B5')).clip(footprint)
    cosdSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).cos()
    sindSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).sin()

    # sat azimuth
    SatAz = ee.Image.constant(img.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B5')).clip(footprint)

    # relative azimuth
    RelAz = SatAz.subtract(SunAz)
    cosdRelAz = RelAz.multiply(pi.divide(ee.Image(180))).cos()

    # Pressure
    P = (ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM)).pow(5.25588)).multiply(
        0.01)).multiply(mask)
    Po = ee.Image(1013.25)

    # esun
    ESUN = ee.Image(ee.Array([ee.Image(img.get('SOLAR_IRRADIANCE_B1')),
                              ee.Image(img.get('SOLAR_IRRADIANCE_B2')),
                              ee.Image(img.get('SOLAR_IRRADIANCE_B3')),
                              ee.Image(img.get('SOLAR_IRRADIANCE_B4')),
                              ee.Image(img.get('SOLAR_IRRADIANCE_B5')),
                              ee.Image(img.get('SOLAR_IRRADIANCE_B6')),
                              ee.Image(img.get('SOLAR_IRRADIANCE_B7')),
                              ee.Image(img.get('SOLAR_IRRADIANCE_B8')),
                              ee.Image(img.get('SOLAR_IRRADIANCE_B8A')),
                              ee.Image(img.get('SOLAR_IRRADIANCE_B11')),
                              ee.Image(img.get('SOLAR_IRRADIANCE_B2'))]
                             )).toArray().toArray(1)

    ESUN = ESUN.multiply(ee.Image(1))

    ESUNImg = ESUN.arrayProject([0]).arrayFlatten([bands])

    # create empty array for the images
    imgArr = rescale.select(bands).toArray().toArray(1)

    # pTOA to Ltoa
    Ltoa = imgArr.multiply(ESUN).multiply(cosdSunZe).divide(pi.multiply(d.pow(2)))

    # band centers
    bandCenter = ee.Image(443).divide(1000).addBands(ee.Image(490).divide(1000)) \
        .addBands(ee.Image(560).divide(1000)) \
        .addBands(ee.Image(665).divide(1000)) \
        .addBands(ee.Image(705).divide(1000)) \
        .addBands(ee.Image(740).divide(1000)) \
        .addBands(ee.Image(783).divide(1000)) \
        .addBands(ee.Image(842).divide(1000)) \
        .addBands(ee.Image(865).divide(1000)) \
        .addBands(ee.Image(1610).divide(1000)) \
        .addBands(ee.Image(2190).divide(1000)) \
        .toArray().toArray(1)

    # ozone coefficients
    koz = ee.Image(0.0039).addBands(ee.Image(0.0213)) \
        .addBands(ee.Image(0.1052)) \
        .addBands(ee.Image(0.0505)) \
        .addBands(ee.Image(0.0205)) \
        .addBands(ee.Image(0.0112)) \
        .addBands(ee.Image(0.0075)) \
        .addBands(ee.Image(0.0021)) \
        .addBands(ee.Image(0.0019)) \
        .addBands(ee.Image(0)) \
        .addBands(ee.Image(0)) \
        .toArray().toArray(1)

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

    R_theta_SZ_p = (
        ((theta_SZ.multiply(pi.divide(180))).subtract(theta_j.multiply(pi.divide(180)))).tan().pow(2)).divide(
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
    bands_nm = ee.Image(443).addBands(ee.Image(490)) \
        .addBands(ee.Image(560)) \
        .addBands(ee.Image(665)) \
        .addBands(ee.Image(705)) \
        .addBands(ee.Image(740)) \
        .addBands(ee.Image(783)) \
        .addBands(ee.Image(842)) \
        .addBands(ee.Image(865)) \
        .addBands(ee.Image(0)) \
        .addBands(ee.Image(0)) \
        .toArray().toArray(1)

    # Lam in SWIR bands
    Lam_10 = LrcImg.select('B11')
    Lam_11 = LrcImg.select('B12')

    # Calculate aerosol type
    eps = (
        (((Lam_11).divide(ESUNImg.select('B12'))).log()).subtract(
            ((Lam_10).divide(ESUNImg.select('B11'))).log())).divide(
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

    # remote sensing reflectance
    Rrs_coll = (pw.divide(pi).arrayProject([0]).arrayFlatten([bands]).slice(0, 9))

    return (Rrs_coll.set('system:time_start', img.get('system:time_start')))


def chlorophyll(img):
    NDCI_coll = (img.select('B5').subtract(img.select('B4'))).divide(img.select('B5').add(img.select('B4')))
    chlor_a_coll = ee.Image(14.039).add(ee.Image(86.115).multiply(NDCI_coll)).add(
        ee.Image(194.325).multiply(NDCI_coll.pow(ee.Image(2))))
    return (chlor_a_coll.updateMask(chlor_a_coll.lt(100)).set('system:time_start', img.get('system:time_start')))


def secchi(img):

    blueRed_coll = (img.select('B2').divide(img.select('B4'))).log()
    lnMOSD_coll = (ee.Image(1.4856).multiply(blueRed_coll)).add(
        ee.Image(0.2734))  # R2 = 0.8748 with Anthony's in-situ data
    MOSD_coll = ee.Image(10).pow(lnMOSD_coll)
    sd_coll = (ee.Image(0.1777).multiply(MOSD_coll)).add(ee.Image(1.0813))
    return (sd_coll.updateMask(sd_coll.lt(10)).set('system:time_start', img.get('system:time_start')))


def trophicState(img):
    tsi_coll = ee.Image(30.6).add(ee.Image(9.81).multiply(img.log()))
    return (tsi_coll.updateMask(tsi_coll.lt(200)).set('system:time_start', img.get('system:time_start')))


def reclassify(img):
    # Create conditions
    mask1 = img.lt(30)  # (1)
    mask2 = img.gte(30).And(img.lt(40))  # (2)
    mask3 = img.gte(40).And(img.lt(50))  # (3)
    mask4 = img.gte(50).And(img.lt(60))  # (4)
    mask5 = img.gte(60).And(img.lt(70))  # (5)
    mask6 = img.gte(70).And(img.lt(80))  # (6)
    mask7 = img.gte(80)  # (7)

    # Reclassify conditions into new values
    img1 = img.where(mask1.eq(1), 1).mask(mask1)
    img2 = img.where(mask2.eq(1), 2).mask(mask2)
    img3 = img.where(mask3.eq(1), 3).mask(mask3)
    img4 = img.where(mask4.eq(1), 4).mask(mask4)
    img5 = img.where(mask5.eq(1), 5).mask(mask5)
    img6 = img.where(mask6.eq(1), 6).mask(mask6)
    img7 = img.where(mask7.eq(1), 7).mask(mask7)

    # Ouput of reclassified image
    tsi_collR = img1.unmask(img2).unmask(img3).unmask(img4).unmask(img5).unmask(img6).unmask(img7)
    return (tsi_collR.updateMask(tsi_collR.set('system:time_start', img.get('system:time_start'))))


#########################################################/
#########################################################/
#########################################################/
# Collection Processing #

# atmospheric correction
Rrs_coll = FC.map(s2Correction)

# chlorophyll-a
chlor_a_coll = Rrs_coll.map(chlorophyll)

# sd
sd_coll = Rrs_coll.map(secchi)

# tsi
tsi_coll = chlor_a_coll.map(trophicState)

# tsi reclass
tsi_collR = tsi_coll.map(reclassify)

Rsr = Rrs_coll.mean()
Rsr_vis = {'min': 0, 'max': 0.03, 'bands': 'B4,B3,B2'}
Rsr_mapid = Rsr.getMapId(Rsr_vis)

SD = sd_coll.mean()
SD_vis = {'min': 0, 'max': 2, 'palette': '#800000,#FF9700,#7BFF7B,#0080FF,#000080'}
SD_mapid = SD.getMapId(SD_vis)

TSI = tsi_coll.mean()
TSI_vis = {'min': 30, 'max': 80,'palette':'darkblue,blue,cyan,limegreen,yellow,orange,orangered,darkred'}
TSI_mapid = TSI.getMapId(TSI_vis)

TSI_R = tsi_collR.mean()
TSI_R_vis = {'min': 1, 'max': 7, 'palette': 'purple,blue,limegreen,yellow,orange,orangered,darkred'}
TSI_R_mapid = TSI_R.getMapId(TSI_R_vis)

print 'RSR',Rsr_mapid
print 'SD',SD_mapid
print 'TSI',TSI_mapid
print 'TSI R',TSI_R_mapid

# #########################################################/
# #########################################################/
# #########################################################/
# Map Layers #
# Map.addLayer(mask, {}, 'mask', false)
# Map.addLayer(Rrs_coll.mean(), {min: 0, max: 0.03, bands: ['B4', 'B3', 'B2']}, 'Mean RGB', false)
# Map.addLayer(chlor_a_coll.mean(), {min: 0, max: 40,
#                                    palette: ['darkblue', 'blue', 'cyan', 'limegreen', 'yellow', 'orange', 'orangered',
#                                              'darkred']}, 'Mean chlor-a', false)
# Map.addLayer(sd_coll.mean(), {min: 0, max: 2, palette: ['800000', 'FF9700', '7BFF7B', '0080FF', '000080']}, 'Mean Zsd',
#              false)
# Map.addLayer(tsi_coll.mean(), {min: 30, max: 80,
#                                palette: ['darkblue', 'blue', 'cyan', 'limegreen', 'yellow', 'orange', 'orangered',
#                                          'darkred']}, 'Mean TSI', false)
# Map.addLayer(tsi_collR.mode(),
#              {min: 1, max: 7, palette: ['purple', 'blue', 'limegreen', 'yellow', 'orange', 'orangered', 'darkred']},
#              'Mode TSI Class', true)

#########################################################/
#########################################################/
#########################################################/
# Time Series #

# # Chlorophyll-a time series
# chlorTimeSeries = ui.Chart.image.seriesByRegion(
# chlor_a_coll, geometry, ee.Reducer.mean())
# .setChartType('ScatterChart')
#     .setOptions({
#     title: 'Mean Chlorphyll-a',
#     vAxis: {title: 'Chlor-a [micrograms/L]'},
#     lineWidth: 1,
#     pointSize: 4,
# })
#
# # SD time series
# sdTimeSeries = ui.Chart.image.seriesByRegion(
# sd_coll, geometry, ee.Reducer.mean())
# .setChartType('ScatterChart')
#     .setOptions({
#     title: 'Mean Secchi Depth',
#     vAxis: {title: 'Zsd [m]'},
#     lineWidth: 1,
#     pointSize: 4,
# })
#
# # TSI time series
# tsiTimeSeries = ui.Chart.image.seriesByRegion(
# tsi_coll, geometry, ee.Reducer.mean())
# .setChartType('ScatterChart')
#     .setOptions({
#     title: 'Mean Trophic State Index',
#     vAxis: {title: 'TSI [1-100]'},
#     lineWidth: 1,
#     pointSize: 4,
# })
#
# # TSI Reclass time series
# tsiRTimeSeries = ui.Chart.image.seriesByRegion(
# tsi_collR, geometry, ee.Reducer.mode())
# .setChartType('ScatterChart')
#     .setOptions({
#     title: 'Mode Trophic State Index Class',
#     vAxis: {title: 'TSI Class'},
#     lineWidth: 1,
#     pointSize: 4,
# })
#
# print(chlorTimeSeries)
# print(sdTimeSeries)
# print(tsiTimeSeries)
# print(tsiRTimeSeries)


