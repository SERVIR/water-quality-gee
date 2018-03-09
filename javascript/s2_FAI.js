// This script processes a single Sentinel-2 TOA to Rayleigh corrected reflectances.
// The floating algal index (FAI) is calculated from Rrc over water bodies.
// How to use:
  // (1) If a "geometry" variable exists in the imports window, delete it.
  // (2) Within the map, select the point button near the top left side.
  // (3) Create a new point by clicking on a location of interest and click run.
  // (4) The "available imagery" ImageCollection within the console displays all available imagery
  //     over that location with the necessary filters described by the user below.
  // (5) Expand the features within the "available imagery" image collection and select an image
  //     by highlighting the FILE_ID and pasting it here:
         var s2 = ee.Image('COPERNICUS/S2/20180216T075011_20180216T080852_T36MXE')
  // (6) Click "Run"
  // (7) Export the image to your Google Drive by clicking on the "tasks" tab and clicking "RUN", be sure to specify
  //     the proper folder.
  
// Author: Benjamin Page //
// Citations: 
// Page, B.P., Kumar, A. and Mishra, D.R., 2018. A novel cross-satellite based assessment of the spatio-temporal development of a cyanobacterial harmful algal bloom. International Journal of Applied Earth Observation and Geoinformation, 66, pp.69-81.
  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// User Input //

// begin date
var iniDate = '2015-01-01';

// end date
var endDate = '2018-02-28';

var cloudPerc = 5

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Import Collections //

// sentinel-2
var MSI = ee.ImageCollection('COPERNICUS/S2');

// toms / omi
var ozone = ee.ImageCollection('TOMS/MERGED');

var SRP = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR');

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// water mask
var startMonth = 5;
var endMonth = 9;
var startYear = 2013;
var endYear = 2017;

var forMask = SRP.filterBounds(geometry).select('B6').filterMetadata('CLOUD_COVER', "less_than", 10).filter(ee.Filter.calendarRange(startMonth, endMonth, 'month')).filter(ee.Filter.calendarRange(startYear, endYear, 'year'));
var mask = ee.Image(forMask.select('B6').median().lt(300)) 
mask = mask.updateMask(mask)

// constants
var pi = ee.Image(3.141592);
var imDate = s2.date();

// filter sentinel 2 collection
var FC = MSI.filterDate(iniDate, endDate).filterBounds(geometry).filterMetadata('CLOUDY_PIXEL_PERCENTAGE', "less_than", cloudPerc);
print(FC, 'Sentinel 2 Collection')
print(imDate, 'image date')

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MSI Atmospheric Correction //

var bands = ['B1','B2','B3','B4','B5','B6','B7', 'B8', 'B8A', 'B11', 'B12'];

// rescale
var rescale = ee.Image(s2.divide(10000).copyProperties(s2)).select(bands);

// tile footprint
var footprint = s2.geometry()
Map.centerObject(footprint)

// dem
var DEM = ee.Image('USGS/SRTMGL1_003').clip(footprint); 

// ozone
var DU = ee.Image(ozone.filterDate(iniDate,endDate).filterBounds(geometry).mean());

//Julian Day
var imgDate = ee.Date(s2.get('system:time_start'));
var FOY = ee.Date.fromYMD(imgDate.get('year'),1,1);
var JD = imgDate.difference(FOY,'day').int().add(1);

// earth-sun distance
var myCos = ((ee.Image(0.0172).multiply(ee.Image(JD).subtract(ee.Image(2)))).cos()).pow(2)
var cosd = myCos.multiply(pi.divide(ee.Image(180))).cos();
var d = ee.Image(1).subtract(ee.Image(0.01673)).multiply(cosd).clip(footprint)

// sun azimuth
var SunAz = ee.Image.constant(s2.get('MEAN_SOLAR_AZIMUTH_ANGLE')).clip(footprint);

// sun zenith
var SunZe = ee.Image.constant(s2.get('MEAN_SOLAR_ZENITH_ANGLE')).clip(footprint);
var cosdSunZe = SunZe.multiply(pi.divide(ee.Image(180))).cos(); // in degrees
var sindSunZe = SunZe.multiply(pi.divide(ee.Image(180))).sin(); // in degrees

// sat zenith
var SatZe = ee.Image.constant(s2.get('MEAN_INCIDENCE_ZENITH_ANGLE_B5')).clip(footprint);
var cosdSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).cos();
var sindSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).sin();
  
// sat azimuth
var SatAz = ee.Image.constant(s2.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B5')).clip(footprint);

// relative azimuth
var RelAz = SatAz.subtract(SunAz);
var cosdRelAz = RelAz.multiply(pi.divide(ee.Image(180))).cos();

// Pressure
var P = ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM)).pow(5.25588)).multiply(0.01)
var Po = ee.Image(1013.25);

// esun
var ESUN = ee.Image(ee.Array([ee.Image(s2.get('SOLAR_IRRADIANCE_B1')),
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
                  )).toArray().toArray(1);

ESUN = ESUN.multiply(ee.Image(1))

var ESUNImg = ESUN.arrayProject([0]).arrayFlatten([bands]);

// create empty array for the images
var imgArr = rescale.select(bands).toArray().toArray(1);

// pTOA to Ltoa
var Ltoa = imgArr.multiply(ESUN).multiply(cosdSunZe).divide(pi.multiply(d.pow(2)));

// band centers
var bandCenter = ee.Image(443).divide(1000).addBands(ee.Image(490).divide(1000))
                                        .addBands(ee.Image(560).divide(1000))
                                        .addBands(ee.Image(665).divide(1000))
                                        .addBands(ee.Image(705).divide(1000))
                                        .addBands(ee.Image(740).divide(1000))
                                        .addBands(ee.Number(783).divide(1000))
                                        .addBands(ee.Number(842).divide(1000))
                                        .addBands(ee.Number(865).divide(1000))
                                        .addBands(ee.Number(1610).divide(1000))
                                        .addBands(ee.Number(2190).divide(1000))
                                        .toArray().toArray(1);

// ozone coefficients
var koz = ee.Image(0.0039).addBands(ee.Image(0.0213))
                        .addBands(ee.Image(0.1052))
                        .addBands(ee.Image(0.0505))
                        .addBands(ee.Image(0.0205))
                        .addBands(ee.Image(0.0112))
                        .addBands(ee.Image(0.0075))
                        .addBands(ee.Image(0.0021))
                        .addBands(ee.Image(0.0019))                          
                        .addBands(ee.Image(0))
                        .addBands(ee.Image(0))
                        .toArray().toArray(1);

//////////////////////////////////////

// Calculate ozone optical thickness
var Toz = koz.multiply(DU).divide(ee.Image(1000));

// Calculate TOA radiance in the absense of ozone
var Lt = Ltoa.multiply(((Toz)).multiply((ee.Image(1).divide(cosdSunZe)).add(ee.Image(1).divide(cosdSatZe))).exp());

// Rayleigh optical thickness
var Tr = (P.divide(Po)).multiply(ee.Image(0.008569).multiply(bandCenter.pow(-4))).multiply((ee.Image(1).add(ee.Image(0.0113).multiply(bandCenter.pow(-2))).add(ee.Image(0.00013).multiply(bandCenter.pow(-4)))));

// Specular reflection (s- and p- polarization states)
var theta_V = ee.Image(0.0000000001);
var sin_theta_j = sindSunZe.divide(ee.Image(1.333));

var theta_j = sin_theta_j.asin().multiply(ee.Image(180).divide(pi));

var theta_SZ = SunZe;

var R_theta_SZ_s = (((theta_SZ.multiply(pi.divide(ee.Image(180)))).subtract(theta_j.multiply(pi.divide(ee.Image(180))))).sin().pow(2)).divide((((theta_SZ.multiply(pi.divide(ee.Image(180)))).add(theta_j.multiply(pi.divide(ee.Image(180))))).sin().pow(2)));

var R_theta_V_s = ee.Image(0.0000000001);

var R_theta_SZ_p = (((theta_SZ.multiply(pi.divide(180))).subtract(theta_j.multiply(pi.divide(180)))).tan().pow(2)).divide((((theta_SZ.multiply(pi.divide(180))).add(theta_j.multiply(pi.divide(180)))).tan().pow(2)));

var R_theta_V_p = ee.Image(0.0000000001);

var R_theta_SZ = ee.Image(0.5).multiply(R_theta_SZ_s.add(R_theta_SZ_p));

var R_theta_V = ee.Image(0.5).multiply(R_theta_V_s.add(R_theta_V_p));
  
// Sun-sensor geometry
var theta_neg = ((cosdSunZe.multiply(ee.Image(-1))).multiply(cosdSatZe)).subtract((sindSunZe).multiply(sindSatZe).multiply(cosdRelAz));

var theta_neg_inv = theta_neg.acos().multiply(ee.Image(180).divide(pi));

var theta_pos = (cosdSunZe.multiply(cosdSatZe)).subtract(sindSunZe.multiply(sindSatZe).multiply(cosdRelAz));

var theta_pos_inv = theta_pos.acos().multiply(ee.Image(180).divide(pi));

var cosd_tni = theta_neg_inv.multiply(pi.divide(180)).cos(); // in degrees

var cosd_tpi = theta_pos_inv.multiply(pi.divide(180)).cos(); // in degrees

var Pr_neg = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tni.pow(2))));

var Pr_pos = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tpi.pow(2))));
  
// Rayleigh scattering phase function
var Pr = Pr_neg.add((R_theta_SZ.add(R_theta_V)).multiply(Pr_pos)); // for water Rayleigh correction
//var Pr = ee.Image(1); // for terrestrial Rayleigh correction

// rayleigh radiance contribution
var denom = ee.Image(4).multiply(pi).multiply(cosdSatZe);
var Lr = (ESUN.multiply(Tr)).multiply(Pr.divide(denom));

// rayleigh corrected radiance
var Lrc = Lt.subtract(Lr);
var LrcImg = Lrc.arrayProject([0]).arrayFlatten([bands]);

// rayleigh corrected reflectance
var prc = (Lrc.multiply(pi).multiply(d.pow(2)).divide(ESUN.multiply(cosdSunZe)));
var prcImg = prc.arrayProject([0]).arrayFlatten([bands]);
print(prcImg, 'prc')

////////////////////////////////////////////////
// models //

// Calculate FAI
var NIRprime = (prcImg.select('B4')).add((prcImg.select('B11').subtract(prcImg.select('B4'))).multiply((ee.Image(865).subtract(ee.Image(665))).divide((ee.Image(1610).subtract(ee.Image(665))))));
var FAI = ((prcImg.select('B8A').subtract(NIRprime))).multiply(mask);

// Map Layers
Map.addLayer(footprint, {}, 'footprint', true);
Map.addLayer(prcImg.multiply(mask), {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.4}, 'prc rgb', false)
Map.addLayer(FAI, {min: -0.05, max: 0.2, palette: ['000080','0080FF','7BFF7B','FF9700','800000']}, 'FAI', true);

// export image //

Export.image.toDrive({
  image: FAI,
  description: 's2_FAI',
  scale: 30,
  region: footprint
});

