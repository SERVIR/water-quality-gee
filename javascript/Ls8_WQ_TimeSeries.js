// This script processes and charts WQ parameter valuse from a collection of Lansdat 8 archive for a particular pixel throughout time.
// WQ parameters such secchi depth, tropshic state index, and lake temperature
// How to use:
  // (1) If a "geometry" variable exists in the imports window, delete it.
  // (2) Within the map, select the point button near the top left side.
  // (3) Create a new point by clicking on a location of interest.
  // (4) Adjust the time frame you wish to browse by adding here:
        // begin date
        var iniDate = '2016-04-01';
        // end date
        var endDate = '2016-10-31';
  // (5) Adjust a cloud % threshold here:
        var oliCloudPerc = 10
  // (5) Click Run
  // (5) The "available imagery" ImageCollection within the console displays all available imagery. If you wish to investigate a single
  //     image from this collection, highlight the FILE_ID and copy and paste it into the proper location within the "Ls8 Single Image" script.
  // (6) The charts will appear on the right hand side of the interface, and will display the mean map in the map user interface.
  // (6) Click "Run"
  // (7) Export each chart as a csv by clicking on the export icon near the top-right of the chart.
  
// Author: Benjamin Page //
// Citations: 
// Page, B.P. and Mishra, D.R., 2018. A modified atmospheric correction when coupling sentinel-2 and landsat-8 for inland water quality monitoring

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Map.centerObject(geometry, 7)

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

// filter landsat 8 collection
var FC_OLI = OLI_DN.filterDate(iniDate, endDate).filterBounds(geometry).filterMetadata('CLOUD_COVER', "less_than", oliCloudPerc);
print(FC_OLI, 'Available Imagery')

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Processing Collection //
var Rrs_coll = FC_OLI.map(l8Correction);

var sd_coll = Rrs_coll.map(secchi);

var tsi_coll = sd_coll.map(trophicState)

var tsi_collR = tsi_coll.map(reclassify);

var lst_coll = FC_OLI.map(LST)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mapping Functions //

function l8Correction(img){

// tile geometry

var l8Footprint = img.geometry()

// dem
var DEM_OLI = ee.Image('USGS/SRTMGL1_003').clip(l8Footprint); 

// ozone
var DU_OLI = ee.Image(ozone.filterDate(iniDate,endDate).filterBounds(l8Footprint).mean());

//Julian Day
var imgDate_OLI = ee.Date(img.get('system:time_start'));
var FOY_OLI = ee.Date.fromYMD(imgDate_OLI.get('year'),1,1);
var JD_OLI = imgDate_OLI.difference(FOY_OLI,'day').int().add(1); 

// Earth-Sun distance
var d_OLI = ee.Image.constant(img.get('EARTH_SUN_DISTANCE'));

//Sun elevation
var SunEl_OLI = ee.Image.constant(img.get('SUN_ELEVATION'));

//Sun azimuth
var SunAz_OLI = ee.Image.constant(img.get('SUN_AZIMUTH'));

//Satellite zenith
var SatZe_OLI = ee.Image(0.0)
var cosdSatZe_OLI = (SatZe_OLI).multiply(pi.divide(ee.Image(180))).cos();
var sindSatZe_OLI = (SatZe_OLI).multiply(pi.divide(ee.Image(180))).sin();

//Satellite azimuth
var SatAz_OLI = ee.Image(0.0).clip(l8Footprint);

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
var rad_mult_OLI = ee.Image(ee.Array([ee.Image(img.get('RADIANCE_MULT_BAND_1')),
                        ee.Image(img.get('RADIANCE_MULT_BAND_2')),
                        ee.Image(img.get('RADIANCE_MULT_BAND_3')),
                        ee.Image(img.get('RADIANCE_MULT_BAND_4')),
                        ee.Image(img.get('RADIANCE_MULT_BAND_5')),
                        ee.Image(img.get('RADIANCE_MULT_BAND_6')),
                        ee.Image(img.get('RADIANCE_MULT_BAND_7'))]
                        )).toArray(1);

// radiance add band                         
var rad_add_OLI = ee.Image(ee.Array([ee.Image(img.get('RADIANCE_ADD_BAND_1')),
                        ee.Image(img.get('RADIANCE_ADD_BAND_2')),
                        ee.Image(img.get('RADIANCE_ADD_BAND_3')),
                        ee.Image(img.get('RADIANCE_ADD_BAND_4')),
                        ee.Image(img.get('RADIANCE_ADD_BAND_5')),
                        ee.Image(img.get('RADIANCE_ADD_BAND_6')),
                        ee.Image(img.get('RADIANCE_ADD_BAND_7'))]
                        )).toArray(1);

//create an empty image to save new radiance bands to
var imgArr_OLI = img.select(bands_OLI).toArray().toArray(1);
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
var prcImg_OLI = prc_OLI.arrayProject([0]).arrayFlatten([bands_OLI]);

// Aerosol Correction // 
// Bands in nm
var bands_nm_OLI = ee.Image(443).addBands(ee.Image(483))
                            .addBands(ee.Image(561))
                            .addBands(ee.Image(655))
                            .addBands(ee.Image(865))
                            .addBands(ee.Image(0))
                            .addBands(ee.Image(0))
                            .toArray().toArray(1);

// Lam in SWIR bands
var Lam_6_OLI = LrcImg_OLI.select('B6');
var Lam_7_OLI = LrcImg_OLI.select('B7');

// Calculate aerosol type
var eps_OLI = (((((Lam_7_OLI).divide(ESUNImg_OLI.select('B7'))).log()).subtract(((Lam_6_OLI).divide(ESUNImg_OLI.select('B6'))).log())).divide(ee.Image(2201).subtract(ee.Image(1609)))).multiply(mask);

// Calculate multiple scattering of aerosols for each band
var Lam_OLI = (Lam_7_OLI).multiply(((ESUN_OLI).divide(ESUNImg_OLI.select('B7')))).multiply((eps_OLI.multiply(ee.Image(-1))).multiply((bands_nm_OLI.divide(ee.Image(2201)))).exp());

// diffuse transmittance
var trans_OLI = Tr_OLI.multiply(ee.Image(-1)).divide(ee.Image(2)).multiply(ee.Image(1).divide(cosdSatZe_OLI)).exp();

// Compute water-leaving radiance
var Lw_OLI = Lrc_OLI.subtract(Lam_OLI).divide(trans_OLI);

// water-leaving reflectance
var pw_OLI = (Lw_OLI.multiply(pi).multiply(d_OLI.pow(2)).divide(ESUN_OLI.multiply(cosdSunZe_OLI)));
var pwImg_OLI = pw_OLI.arrayProject([0]).arrayFlatten([bands_OLI]);

// Rrs
var Rrs_coll = (pw_OLI.divide(pi).arrayProject([0]).arrayFlatten([bands_OLI]).slice(0,5));

return(Rrs_coll.set('system:time_start',img.get('system:time_start')));

}
function LST(img){
  var TIRS_1 = img.select('B10');
  
  var b10_add_band = ee.Number(img.get('RADIANCE_ADD_BAND_10'));
  var b10_mult_band = ee.Number(img.get('RADIANCE_MULT_BAND_10'));
  
  var Oi = ee.Number(0.29);
  
  var TIRS_cal = (TIRS_1.multiply(b10_mult_band).add(b10_add_band).subtract(Oi)); //.multiply(lakes);
  
  var K1_b10 = ee.Number(TIRS_1.get('K1_CONSTANT_BAND_10'));
  var K2_b10 = ee.Number(TIRS_1.get('K2_CONSTANT_BAND_10'));
  
  var lst_coll = ((ee.Image(K2_b10).divide(((ee.Image(K1_b10).divide(ee.Image(TIRS_cal))).add(ee.Image(1))).log())).subtract(ee.Image(273))).multiply(mask); // Celsius
  lst_coll = (ee.Image(0.7745).multiply(lst_coll)).add(ee.Image(9.6502));
  
  // Define a square kernel with radius 1.5
  var sq_kernel = ee.Kernel.square(1.5,'pixels');

  // focal convolution 
  lst_coll = lst_coll.focal_median({kernel: sq_kernel, iterations: 1});
  
  return(lst_coll.set('system:time_start',img.get('system:time_start')))

}
function secchi(img){
  var blueRed_coll = (img.select('B2').divide(img.select('B4'))).log()
  var lnMOSD_coll = (ee.Image(1.4856).multiply(blueRed_coll)).add(ee.Image(0.2734)); // R2 = 0.8748 with Anthony's in-situ data
  var MOSD_coll = ee.Image(10).pow(lnMOSD_coll);
  var sd_coll = (ee.Image(0.1777).multiply(MOSD_coll)).add(ee.Image(1.0813));
  return(sd_coll.updateMask(sd_coll.lt(10)).set('system:time_start',img.get('system:time_start')))
}
function trophicState(img){
  var tsi_coll =  ee.Image(60).subtract(ee.Image(14.41).multiply(img.log()));
  return(tsi_coll.updateMask(tsi_coll.lt(200)).set('system:time_start',img.get('system:time_start')))

}
function reclassify(img){
  
  // Create conditions
  var mask1 = img.lt(30); // (1)
  var mask2 = img.gte(30).and(img.lt(40));// (2)
  var mask3 = img.gte(40).and(img.lt(50)); // (3)
  var mask4 = img.gte(50).and(img.lt(60)); // (4)
  var mask5 = img.gte(60).and(img.lt(70)); // (5)
  var mask6 = img.gte(70).and(img.lt(80)); // (6)
  var mask7 = img.gte(70); // (7)

  // Reclassify conditions into new values
  var img1 = img.where(mask1.eq(1), 1).mask(mask1);
  var img2 = img.where(mask2.eq(1), 2).mask(mask2);
  var img3 = img.where(mask3.eq(1), 3).mask(mask3);
  var img4 = img.where(mask4.eq(1), 4).mask(mask4);
  var img5 = img.where(mask5.eq(1), 5).mask(mask5);
  var img6 = img.where(mask6.eq(1), 6).mask(mask6);
  var img7 = img.where(mask7.eq(1), 7).mask(mask7);

  // Ouput of reclassified image
  var tsi_collR = img1.unmask(img2).unmask(img3).unmask(img4).unmask(img5).unmask(img6).unmask(img7);
  return(tsi_collR.updateMask(tsi_collR.set('system:time_start',img.get('system:time_start'))))
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Map Layers //
Map.addLayer(mask, {}, 'mask', false)

// l8
Map.addLayer(Rrs_coll.mean(), {bands: ['B4', 'B3', 'B1'], min: 0, max: 0.04}, 'l8 Rrs RGB', false);
Map.addLayer(lst_coll.mean(), {min: 23, max: 27, palette: ['darkblue', 'blue', 'white', 'red', 'darkred']}, 'LST [c]', false)
Map.addLayer(sd_coll.mean(), {min: 0, max: 3, palette: ['800000', 'FF9700', '7BFF7B', '0080FF', '000080']}, 'Zsd [m]', false);
Map.addLayer(tsi_coll.mean(), {min: 30, max: 70, palette: ['blue', 'cyan', 'limegreen', 'yellow', 'darkred']}, 'TSI [0-100]', false);
Map.addLayer(tsi_collR.mode(), {min: 1, max: 7, palette: ['purple', 'blue', 'limegreen', 'yellow', 'orange', 'orangered', 'darkred']}, 'TSI Reclass', true);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Time Series //

// LST time series
var lstTimeSeries = ui.Chart.image.seriesByRegion(
    lst_coll, geometry, ee.Reducer.mean())
        .setChartType('ScatterChart')
        .setOptions({
          title: 'Mean LST',
          vAxis: {title: 'LST [c]'},
          lineWidth: 1,
          pointSize: 4,
});
  
// SD time series
var sdTimeSeries = ui.Chart.image.seriesByRegion(
    sd_coll, geometry, ee.Reducer.mean())
        .setChartType('ScatterChart')
        .setOptions({
          title: 'Mean Secchi Depth',
          vAxis: {title: 'Zsd [m]'},
          lineWidth: 1,
          pointSize: 4,
});


// TSI time series
var tsiTimeSeries = ui.Chart.image.seriesByRegion(
    tsi_coll, geometry, ee.Reducer.mean())
        .setChartType('ScatterChart')
        .setOptions({
          title: 'Mean Trophic State Index',
          vAxis: {title: 'TSI [1-100]'},
          lineWidth: 1,
          pointSize: 4,
});

// TSI Reclass time series
var tsiRTimeSeries = ui.Chart.image.seriesByRegion(
    tsi_collR, geometry, ee.Reducer.mean())
        .setChartType('ScatterChart')
        .setOptions({
          title: 'Mode Trophic State Index Class',
          vAxis: {title: 'TSI Class'},
          lineWidth: 1,
          pointSize: 4,
});

print(lstTimeSeries)
print(sdTimeSeries)
print(tsiTimeSeries)
print(tsiRTimeSeries)