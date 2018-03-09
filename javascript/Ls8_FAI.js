// This script processes a single Landsat 8 Tier 1 Raw image to Rayleigh corrected reflectances.
// The floating algal index (FAI) is calculated from Rrc over water bodies.
// How to use:
  // (1) If a "geometry" variable exists in the imports window, delete it.
  // (2) Within the map, select the point button near the top left side.
  // (3) Create a new point by clicking on a location of interest and click run.
  // (4) The "available imagery" ImageCollection within the console displays all available imagery
  //     over that location with the necessary filters described by the user below.
  // (5) Expand the features within the "available imagery" image collection and select an image
  //     by highlighting the FILE_ID and pasting it here:
         var l8 = ee.Image('LANDSAT/LC08/C01/T1/LC08_170060_20171226');
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
var iniDate = '2015-05-01'; 

// end date
var endDate = '2018-03-31';

var oliCloudPerc = 5

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Import Collections //

// landsat 8 raw dn
var OLI_DN = ee.ImageCollection('LANDSAT/LC08/C01/T1');

// landsat-8 surface reflactance product (for masking purposes)
var SRP = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR');

// toms / omi
var ozone = ee.ImageCollection('TOMS/MERGED');

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filtering Collection and Masking //

var pi = ee.Image(3.141592);

// water mask
var startMonth = 5;
var endMonth = 9;
var startYear = 2013;
var endYear = 2017;

var forMask = SRP.filterBounds(geometry).select('B6').filterMetadata('CLOUD_COVER', "less_than", 10).filter(ee.Filter.calendarRange(startMonth, endMonth, 'month')).filter(ee.Filter.calendarRange(startYear, endYear, 'year'));
var mask = ee.Image(forMask.select('B6').median().lt(600)) 
mask = mask.updateMask(mask)
Map.addLayer(mask, {}, 'mask', false)

// filter landsat 8 collection
var FC_OLI = OLI_DN.filterDate(iniDate, endDate).filterBounds(geometry).filterMetadata('CLOUD_COVER', "less_than", oliCloudPerc);
print(FC_OLI, 'Available Imagery')

print(l8, 'l8 image info')

// oli image date
var oliDate = l8.date();
print(oliDate, 'l8 Date')

var footprint = l8.geometry()
Map.centerObject(footprint)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// single Ls8 image Rayleigh Correction // 

// dem
var DEM_OLI = ee.Image('USGS/SRTMGL1_003').clip(footprint);

// ozone
var DU_OLI = ee.Image(ozone.filterDate(iniDate,endDate).filterBounds(footprint).mean());

//Julian Day
var imgDate_OLI = ee.Date(l8.get('system:time_start'));
var FOY_OLI = ee.Date.fromYMD(imgDate_OLI.get('year'),1,1);
var JD_OLI = imgDate_OLI.difference(FOY_OLI,'day').int().add(1); 

// Earth-Sun distance
var d_OLI = ee.Image.constant(l8.get('EARTH_SUN_DISTANCE'));

//Sun elevation
var SunEl_OLI = ee.Image.constant(l8.get('SUN_ELEVATION'));

//Sun azimuth
var SunAz_OLI = ee.Image.constant(l8.get('SUN_AZIMUTH'));

//Satellite zenith
var SatZe_OLI = ee.Image(0.0)
var cosdSatZe_OLI = (SatZe_OLI).multiply(pi.divide(ee.Image(180))).cos();
var sindSatZe_OLI = (SatZe_OLI).multiply(pi.divide(ee.Image(180))).sin();

//Satellite azimuth
var SatAz_OLI = ee.Image(0.0)

//Sun zenith
var SunZe_OLI = ee.Image(90).subtract(SunEl_OLI)
var cosdSunZe_OLI = SunZe_OLI.multiply(pi.divide(ee.Image.constant(180))).cos(); // in degrees
var sindSunZe_OLI = SunZe_OLI.multiply(pi.divide(ee.Image(180))).sin(); // in degrees

//Relative azimuth
var RelAz_OLI = ee.Image(SunAz_OLI);
var cosdRelAz_OLI = RelAz_OLI.multiply(pi.divide(ee.Image(180))).cos();

//Pressure calculation
var P_OLI = ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM_OLI)).pow(5.25588)).multiply(0.01);
var Po_OLI = ee.Image(1013.25);

// Radiometric Calibration //
//define bands to be converted to radiance
var bands_OLI = ['B1','B2','B3','B4','B5','B6','B7'];

// radiance_mult_bands
var rad_mult_OLI = ee.Image(ee.Array([ee.Image(l8.get('RADIANCE_MULT_BAND_1')),
                        ee.Image(l8.get('RADIANCE_MULT_BAND_2')),
                        ee.Image(l8.get('RADIANCE_MULT_BAND_3')),
                        ee.Image(l8.get('RADIANCE_MULT_BAND_4')),
                        ee.Image(l8.get('RADIANCE_MULT_BAND_5')),
                        ee.Image(l8.get('RADIANCE_MULT_BAND_6')),
                        ee.Image(l8.get('RADIANCE_MULT_BAND_7'))]
                        )).toArray(1);

// radiance add band                         
var rad_add_OLI = ee.Image(ee.Array([ee.Image(l8.get('RADIANCE_ADD_BAND_1')),
                        ee.Image(l8.get('RADIANCE_ADD_BAND_2')),
                        ee.Image(l8.get('RADIANCE_ADD_BAND_3')),
                        ee.Image(l8.get('RADIANCE_ADD_BAND_4')),
                        ee.Image(l8.get('RADIANCE_ADD_BAND_5')),
                        ee.Image(l8.get('RADIANCE_ADD_BAND_6')),
                        ee.Image(l8.get('RADIANCE_ADD_BAND_7'))]
                        )).toArray(1);

//create an empty image to save new radiance bands to
var imgArr_OLI = l8.select(bands_OLI).toArray().toArray(1);
var Ltoa_OLI = imgArr_OLI.multiply(rad_mult_OLI).add(rad_add_OLI)

// esun
var ESUN_OLI = ee.Image.constant(197.24790954589844)
                .addBands(ee.Image.constant(201.98426818847656))
                .addBands(ee.Image.constant(186.12677001953125))
                .addBands(ee.Image.constant(156.95257568359375))
                .addBands(ee.Image.constant(96.04714965820312))
                .addBands(ee.Image.constant(23.8833221450863))
                .addBands(ee.Image.constant(8.04995873449635)).toArray().toArray(1);
ESUN_OLI = ESUN_OLI.multiply(ee.Image(1))

var ESUNImg_OLI = ESUN_OLI.arrayProject([0]).arrayFlatten([bands_OLI]);

// Ozone Correction //
// Ozone coefficients
var koz_OLI = ee.Image.constant(0.0039).addBands(ee.Image.constant(0.0218))
                          .addBands(ee.Image.constant(0.1078))
                          .addBands(ee.Image.constant(0.0608))
                          .addBands(ee.Image.constant(0.0019))
                          .addBands(ee.Image.constant(0))
                          .addBands(ee.Image.constant(0))
                          .toArray().toArray(1);

// Calculate ozone optical thickness
var Toz_OLI = koz_OLI.multiply(DU_OLI).divide(ee.Image.constant(1000));

// Calculate TOA radiance in the absense of ozone
var Lt_OLI = Ltoa_OLI.multiply(((Toz_OLI)).multiply((ee.Image.constant(1).divide(cosdSunZe_OLI)).add(ee.Image.constant(1).divide(cosdSatZe_OLI))).exp());

// Rayleigh optical thickness
var bandCenter_OLI = ee.Image(443).divide(1000).addBands(ee.Image(483).divide(1000))
                                          .addBands(ee.Image(561).divide(1000))
                                          .addBands(ee.Image(655).divide(1000))
                                          .addBands(ee.Image(865).divide(1000))
                                          .addBands(ee.Image(1609).divide(1000))
                                          .addBands(ee.Number(2201).divide(1000))
                                          .toArray().toArray(1);

// create an empty image to save new Tr values to
var Tr_OLI = (P_OLI.divide(Po_OLI)).multiply(ee.Image(0.008569).multiply(bandCenter_OLI.pow(-4))).multiply((ee.Image(1).add(ee.Image(0.0113).multiply(bandCenter_OLI.pow(-2))).add(ee.Image(0.00013).multiply(bandCenter_OLI.pow(-4)))));

// Fresnel Reflection //
// Specular reflection (s- and p- polarization states)
var theta_V_OLI = ee.Image(0.0000000001);
var sin_theta_j_OLI = sindSunZe_OLI.divide(ee.Image(1.333));

var theta_j_OLI = sin_theta_j_OLI.asin().multiply(ee.Image(180).divide(pi));

var theta_SZ_OLI = SunZe_OLI;

var R_theta_SZ_s_OLI = (((theta_SZ_OLI.multiply(pi.divide(ee.Image(180)))).subtract(theta_j_OLI.multiply(pi.divide(ee.Image(180))))).sin().pow(2)).divide((((theta_SZ_OLI.multiply(pi.divide(ee.Image(180)))).add(theta_j_OLI.multiply(pi.divide(ee.Image(180))))).sin().pow(2)));

var R_theta_V_s_OLI = ee.Image(0.0000000001);

var R_theta_SZ_p_OLI = (((theta_SZ_OLI.multiply(pi.divide(180))).subtract(theta_j_OLI.multiply(pi.divide(180)))).tan().pow(2)).divide((((theta_SZ_OLI.multiply(pi.divide(180))).add(theta_j_OLI.multiply(pi.divide(180)))).tan().pow(2)));

var R_theta_V_p_OLI = ee.Image(0.0000000001);

var R_theta_SZ_OLI = ee.Image(0.5).multiply(R_theta_SZ_s_OLI.add(R_theta_SZ_p_OLI));

var R_theta_V_OLI = ee.Image(0.5).multiply(R_theta_V_s_OLI.add(R_theta_V_p_OLI));

// Rayleigh scattering phase function //
// Sun-sensor geometry
var theta_neg_OLI = ((cosdSunZe_OLI.multiply(ee.Image(-1))).multiply(cosdSatZe_OLI)).subtract((sindSunZe_OLI).multiply(sindSatZe_OLI).multiply(cosdRelAz_OLI));

var theta_neg_inv_OLI = theta_neg_OLI.acos().multiply(ee.Image(180).divide(pi));

var theta_pos_OLI = (cosdSunZe_OLI.multiply(cosdSatZe_OLI)).subtract(sindSunZe_OLI.multiply(sindSatZe_OLI).multiply(cosdRelAz_OLI));

var theta_pos_inv_OLI = theta_pos_OLI.acos().multiply(ee.Image(180).divide(pi));

var cosd_tni_OLI = theta_neg_inv_OLI.multiply(pi.divide(180)).cos(); // in degrees

var cosd_tpi_OLI = theta_pos_inv_OLI.multiply(pi.divide(180)).cos(); // in degrees

var Pr_neg_OLI = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tni_OLI.pow(2))));

var Pr_pos_OLI = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tpi_OLI.pow(2))));

// Rayleigh scattering phase function
var Pr_OLI = Pr_neg_OLI.add((R_theta_SZ_OLI.add(R_theta_V_OLI)).multiply(Pr_pos_OLI));

// Calulate Lr,
var denom_OLI = ee.Image(4).multiply(pi).multiply(cosdSatZe_OLI);
var Lr_OLI = (ESUN_OLI.multiply(Tr_OLI)).multiply(Pr_OLI.divide(denom_OLI));

// Rayleigh corrected radiance
var Lrc_OLI = (Lt_OLI.divide(ee.Image(10))).subtract(Lr_OLI);
var LrcImg_OLI = Lrc_OLI.arrayProject([0]).arrayFlatten([bands_OLI]);

// Rayleigh corrected reflectance
var prc_OLI = Lrc_OLI.multiply(pi).multiply(d_OLI.pow(2)).divide(ESUN_OLI.multiply(cosdSunZe_OLI));
var pc = prc_OLI.arrayProject([0]).arrayFlatten([bands_OLI]);

// Calculate FAI
var NIRprime = (pc.select('B4')).add((pc.select('B6').subtract(pc.select('B4'))).multiply((ee.Image(865).subtract(ee.Image(655))).divide((ee.Image(1609).subtract(ee.Image(655))))));
var fai = ((pc.select('B5').subtract(NIRprime))).multiply(mask);

// Map Layers //
// Map Layers
Map.addLayer(footprint, {}, 'footprint', false);
Map.addLayer(pc.multiply(mask), {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.1}, 'pc rgb', false)
Map.addLayer(fai, {min: -0.05, max: 0.2, palette: ['000080','0080FF','7BFF7B','FF9700','800000']}, 'FAI', true);

// Export Image to Drive
Export.image.toDrive({
  image: fai,
  description: 'l8_FAI',
  scale: 30,
  region: footprint
});