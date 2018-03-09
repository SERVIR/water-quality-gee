// This script processes a single Sentinel-2 image of interest.
// WQ parameters such chlor-a, secchi depth, trophic state index are calculated
// How to use:
  // (1) If a "geometry" variable exists in the imports window, delete it.
  // (2) Within the map, select the point button near the top left side.
  // (3) Create a new point by clicking on a location of interest.
  // (4) Adjust the time frame you wish to browse by adding here:
        // begin date
        var iniDate = '2015-05-01';
        // end date
        var endDate = '2018-03-31';
  // (5) Adjust a cloud % threshold here:
        var cloudPerc = 5
  // (6) Click Run
  // (7) The "available imagery" ImageCollection within the console displays all available imagery. If you wish to investigate a single
  //     image from this collection, highlight the FILE_ID in the available imagery features and copy and paste it here:
        var s2 = ee.Image('COPERNICUS/S2/20171128T075251_20171128T080644_T36MXE')
  // (8) The charts will appear on the right hand side of the interface, and will display the mean map in the map user interface.
  // (9) Click "Run"
  // (10) Export each image by clicking on the run button within the "tasks" tab.
  
// Author: Benjamin Page //
// Citations: 
// Page, B.P. and Mishra, D.R., 2018. A modified atmospheric correction when coupling sentinel-2 and landsat-8 for inland water quality monitoring

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Import Collections //

// sentinel-2
var MSI = ee.ImageCollection('COPERNICUS/S2');

// landsat-8 surface reflactance product (for masking purposes)
var SRP = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR');

// toms / omi
var ozone = ee.ImageCollection('TOMS/MERGED');

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

var pi = ee.Image(3.141592);
var imDate = s2.date();

var footprint = s2.geometry()

// water mask
var startMonth = 5;
var endMonth = 9;
var startYear = 2013;
var endYear = 2017;

var forMask = SRP.filterBounds(geometry).select('B6').filterMetadata('CLOUD_COVER', "less_than", 10).filter(ee.Filter.calendarRange(startMonth, endMonth, 'month')).filter(ee.Filter.calendarRange(startYear, endYear, 'year'));
var mask = ee.Image(forMask.select('B6').median().lt(300)) 
mask = mask.updateMask(mask).clip(footprint)

// filter sentinel 2 collection
var FC = MSI.filterDate(iniDate, endDate).filterBounds(geometry).filterMetadata('CLOUDY_PIXEL_PERCENTAGE', "less_than", cloudPerc);
print(FC, 'Sentinel 2 Collection')
print(imDate, 'image date')

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MSI Atmospheric Correction //

var bands = ['B1','B2','B3','B4','B5','B6','B7', 'B8', 'B8A', 'B11', 'B12'];

// rescale*
var rescale = ee.Image(s2.divide(10000).multiply(mask).copyProperties(s2)).select(bands);

// dem*
var DEM = ee.Image('USGS/SRTMGL1_003').clip(footprint); 

// ozone*
var DU = ee.Image(ozone.filterDate(iniDate,endDate).filterBounds(footprint).mean());

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
var P = (ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM)).pow(5.25588)).multiply(0.01)).multiply(mask);
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
var Pr = Pr_neg.add((R_theta_SZ.add(R_theta_V)).multiply(Pr_pos));

// rayleigh radiance contribution
var denom = ee.Image(4).multiply(pi).multiply(cosdSatZe);
var Lr = (ESUN.multiply(Tr)).multiply(Pr.divide(denom));

// rayleigh corrected radiance
var Lrc = Lt.subtract(Lr);
var LrcImg = Lrc.arrayProject([0]).arrayFlatten([bands]);

// Aerosol Correction //

// Bands in nm
var bands_nm = ee.Image(443).addBands(ee.Image(490))
                            .addBands(ee.Image(560))
                            .addBands(ee.Image(665))
                            .addBands(ee.Image(705))
                            .addBands(ee.Image(740))
                            .addBands(ee.Image(783))
                            .addBands(ee.Image(842))
                            .addBands(ee.Image(865))                            
                            .addBands(ee.Image(0))
                            .addBands(ee.Image(0))
                            .toArray().toArray(1);

// Lam in SWIR bands
var Lam_10 = LrcImg.select('B11');
var Lam_11 = LrcImg.select('B12');

// Calculate aerosol type
var eps = ((((Lam_11).divide(ESUNImg.select('B12'))).log()).subtract(((Lam_10).divide(ESUNImg.select('B11'))).log())).divide(ee.Image(2190).subtract(ee.Image(1610)));

// Calculate multiple scattering of aerosols for each band
var Lam = (Lam_11).multiply(((ESUN).divide(ESUNImg.select('B12')))).multiply((eps.multiply(ee.Image(-1))).multiply((bands_nm.divide(ee.Image(2190)))).exp());

// diffuse transmittance
var trans = Tr.multiply(ee.Image(-1)).divide(ee.Image(2)).multiply(ee.Image(1).divide(cosdSatZe)).exp();

// Compute water-leaving radiance
var Lw = Lrc.subtract(Lam).divide(trans);

// water-leaving reflectance
var pw = (Lw.multiply(pi).multiply(d.pow(2)).divide(ESUN.multiply(cosdSunZe)));
var pwImg = pw.arrayProject([0]).arrayFlatten([bands]);

// remote sensing reflectance
var Rrs = (pw.divide(pi).arrayProject([0]).arrayFlatten([bands]).slice(0,9));
print(Rrs, 'Rrs')

/// Bio optical Models ///
// chlor_a
var NDCI = (Rrs.select('B5').subtract(Rrs.select('B4'))).divide(Rrs.select('B5').add(Rrs.select('B4'))) // ndci
var chlor_a = ee.Image(14.039).add(ee.Image(86.115).multiply(NDCI)).add(ee.Image(194.325).multiply(NDCI.pow(ee.Image(2)))); // chlor_a

// SD
var ln_BlueRed = (Rrs.select('B2').divide(Rrs.select('B4'))).log()
var lnMOSD = (ee.Image(1.4856).multiply(ln_BlueRed)).add(ee.Image(0.2734)); // R2 = 0.8748 with in-situ
var MOSD = ee.Image(10).pow(lnMOSD); // log space to (m)
var SD = (ee.Image(0.1777).multiply(MOSD)).add(ee.Image(1.0813));

// tsi
var TSI_c = ee.Image(30.6).add(ee.Image(9.81)).multiply(chlor_a.log())
var TSI_s = ee.Image(60.0).subtract(ee.Image(14.41)).multiply(SD.log())
var TSI = (TSI_c.add(TSI_s)).divide(ee.Image(2));

// tsi reclassified
// Create conditions
  var mask1 = TSI.lt(30); // classical oligitrophy (1)
  var mask2 = TSI.gte(30).and(TSI.lt(40));// (2)
  var mask3 = TSI.gte(40).and(TSI.lt(50)); // (3)
  var mask4 = TSI.gte(50).and(TSI.lt(60)); // (4)
  var mask5 = TSI.gte(60).and(TSI.lt(70)); // (5)
  var mask6 = TSI.gte(70).and(TSI.lt(80)); // (6)
  var mask7 = TSI.gte(80); // (7)


  // Reclassify conditions into new values
  var img1 = TSI.where(mask1.eq(1), 1).mask(mask1);
  var img2 = TSI.where(mask2.eq(1), 2).mask(mask2);
  var img3 = TSI.where(mask3.eq(1), 3).mask(mask3);
  var img4 = TSI.where(mask4.eq(1), 4).mask(mask4);
  var img5 = TSI.where(mask5.eq(1), 5).mask(mask5);
  var img6 = TSI.where(mask6.eq(1), 6).mask(mask6);
  var img7 = TSI.where(mask7.eq(1), 7).mask(mask7);
  

  // Ouput of reclassified image
  var TSI_R = img1.unmask(img2).unmask(img3).unmask(img4).unmask(img5).unmask(img6).unmask(img7);

// Map Layers //
Map.centerObject(footprint, 7);
Map.addLayer(footprint, {}, 'footprint', true) // footprint
Map.addLayer(s2.multiply(mask), {bands: ['B4', 'B3', 'B2'], min: 0, max: 1000}, 'rgb', false) // rgb
Map.addLayer(SD, {min: 0, max: 3, palette: ['darkred', 'orange', 'yellow', 'limegreen', 'cyan' , 'blue']}, 'SD', false);
Map.addLayer(chlor_a, {min: 0, max: 40, palette: ['blue', 'cyan', 'limegreen', 'yellow', 'orange' , 'darkred']}, 'chlor_a', false);
Map.addLayer(TSI, {min: 30, max: 70, palette: ['blue', 'cyan', 'limegreen', 'yellow', 'orange' , 'darkred']}, 'TSI', false);
Map.addLayer(TSI_R, {min: 1, max: 7, palette: ['purple', 'blue', 'limegreen', 'yellow', 'orange', 'orangered', 'darkred']}, 'TSI_R', true);

// EXPORT IMAGES //

// Export SD
Export.image.toDrive({
  image: SD,
  description: 'SD',
  scale: 20,
  region: footprint
});

// Export chlor_a
Export.image.toDrive({
  image: chlor_a,
  description: 'chlor_a',
  scale: 20,
  region: footprint
});

// Export TSI
Export.image.toDrive({
  image: TSI,
  description: 'TSI',
  scale: 20,
  region: footprint
});

// Export TSI_R
Export.image.toDrive({
  image: TSI_R,
  description: 'TSI_R',
  scale: 20,
  region: footprint
});