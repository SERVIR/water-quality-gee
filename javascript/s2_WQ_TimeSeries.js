// This script processes and charts WQ parameter valuse from a collection of Sentinel-2 archive for a particular pixel throughout time.
// WQ parameters such chlor-a, secchi depth, trophic state index
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
  // (5) Click Run
  // (5) The "available imagery" ImageCollection within the console displays all available imagery. If you wish to investigate a single
  //     image from this collection, highlight the FILE_ID and copy and paste it into the proper location within the "s2 Single Image" script.
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

// water mask
var startMonth = 5;
var endMonth = 9;
var startYear = 2013;
var endYear = 2017;

var forMask = SRP.filterBounds(geometry).select('B6').filterMetadata('CLOUD_COVER', "less_than", 10).filter(ee.Filter.calendarRange(startMonth, endMonth, 'month')).filter(ee.Filter.calendarRange(startYear, endYear, 'year'));
var mask = ee.Image(forMask.select('B6').median().lt(300)) 
mask = mask.updateMask(mask)

// filter sentinel 2 collection
var FC = MSI.filterDate(iniDate, endDate).filterBounds(geometry).filterMetadata('CLOUDY_PIXEL_PERCENTAGE', "less_than", cloudPerc);
print(FC, 'Sentinel 2 Collection')

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Collection Processing //

// atmospheric correction
var Rrs_coll = FC.map(s2Correction);


// chlorophyll-a
var chlor_a_coll = Rrs_coll.map(chlorophyll);

// sd
var sd_coll = Rrs_coll.map(secchi);

// tsi
var tsi_coll = chlor_a_coll.map(trophicState);

// tsi reclass
var tsi_collR = tsi_coll.map(reclassify);

// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mapping functions //

function s2Correction(img){

// msi bands 
var bands = ['B1','B2','B3','B4','B5','B6','B7', 'B8', 'B8A', 'B11', 'B12'];

// rescale
var rescale = img.select(bands).divide(10000).multiply(mask)

// tile footprint
var footprint = rescale.geometry()

// dem
var DEM = ee.Image('USGS/SRTMGL1_003').clip(footprint); 

// ozone
var DU = ee.Image(ozone.filterDate(iniDate,endDate).filterBounds(footprint).mean());

//Julian Day
var imgDate = ee.Date(img.get('system:time_start'));
var FOY = ee.Date.fromYMD(imgDate.get('year'),1,1);
var JD = imgDate.difference(FOY,'day').int().add(1);

// earth-sun distance
var myCos = ((ee.Image(0.0172).multiply(ee.Image(JD).subtract(ee.Image(2)))).cos()).pow(2)
var cosd = myCos.multiply(pi.divide(ee.Image(180))).cos();
var d = ee.Image(1).subtract(ee.Image(0.01673)).multiply(cosd).clip(footprint)

// sun azimuth
var SunAz = ee.Image.constant(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')).clip(footprint);

// sun zenith
var SunZe = ee.Image.constant(img.get('MEAN_SOLAR_ZENITH_ANGLE')).clip(footprint);
var cosdSunZe = SunZe.multiply(pi.divide(ee.Image(180))).cos(); // in degrees
var sindSunZe = SunZe.multiply(pi.divide(ee.Image(180))).sin(); // in degrees

// sat zenith
var SatZe = ee.Image.constant(img.get('MEAN_INCIDENCE_ZENITH_ANGLE_B5')).clip(footprint);
var cosdSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).cos();
var sindSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).sin();
  
// sat azimuth
var SatAz = ee.Image.constant(img.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B5')).clip(footprint);

// relative azimuth
var RelAz = SatAz.subtract(SunAz);
var cosdRelAz = RelAz.multiply(pi.divide(ee.Image(180))).cos();

// Pressure
var P = (ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM)).pow(5.25588)).multiply(0.01)).multiply(mask);
var Po = ee.Image(1013.25);

// esun
var ESUN = ee.Image(ee.Array([ee.Image(img.get('SOLAR_IRRADIANCE_B1')),
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
                                        .addBands(ee.Image(783).divide(1000))
                                        .addBands(ee.Image(842).divide(1000))
                                        .addBands(ee.Image(865).divide(1000))
                                        .addBands(ee.Image(1610).divide(1000))
                                        .addBands(ee.Image(2190).divide(1000))
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

// remote sensing reflectance
var Rrs_coll = (pw.divide(pi).arrayProject([0]).arrayFlatten([bands]).slice(0,9));

return(Rrs_coll.set('system:time_start',img.get('system:time_start')));

}
function chlorophyll(img){
  var NDCI_coll = (img.select('B5').subtract(img.select('B4'))).divide(img.select('B5').add(img.select('B4')));
  var chlor_a_coll = ee.Image(14.039).add(ee.Image(86.115).multiply(NDCI_coll)).add(ee.Image(194.325).multiply(NDCI_coll.pow(ee.Image(2))));
  return(chlor_a_coll.updateMask(chlor_a_coll.lt(100)).set('system:time_start',img.get('system:time_start')))
}
function secchi(img){
  var blueRed_coll = (img.select('B2').divide(img.select('B4'))).log()
  var lnMOSD_coll = (ee.Image(1.4856).multiply(blueRed_coll)).add(ee.Image(0.2734)); // R2 = 0.8748 with Anthony's in-situ data
  var MOSD_coll = ee.Image(10).pow(lnMOSD_coll);
  var sd_coll = (ee.Image(0.1777).multiply(MOSD_coll)).add(ee.Image(1.0813));
  return(sd_coll.updateMask(sd_coll.lt(10)).set('system:time_start',img.get('system:time_start')))
}
function trophicState(img){
  var tsi_coll =  ee.Image(30.6).add(ee.Image(9.81).multiply(img.log()));
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
  var mask7 = img.gte(80); // (7)

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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Time Series //

// Chlorophyll-a time series
var chlorTimeSeries = ui.Chart.image.seriesByRegion(
    chlor_a_coll, geometry, ee.Reducer.mean())
        .setChartType('ScatterChart')
        .setOptions({
          title: 'Mean Chlorphyll-a',
          vAxis: {title: 'Chlor-a [micrograms/L]'},
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
    tsi_collR, geometry, ee.Reducer.mode())
        .setChartType('ScatterChart')
        .setOptions({
          title: 'Mode Trophic State Index Class',
          vAxis: {title: 'TSI Class'},
          lineWidth: 1,
          pointSize: 4,
});


print(chlorTimeSeries)
print(sdTimeSeries)
print(tsiTimeSeries)
print(tsiRTimeSeries)

// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Map Layers //
Map.addLayer(mask, {}, 'mask', false)
Map.addLayer(Rrs_coll.mean(), {min: 0, max: 0.03, bands: ['B4', 'B3', 'B2']}, 'Mean RGB', false);
Map.addLayer(chlor_a_coll.mean(), {min: 0, max: 40, palette: ['darkblue','blue','cyan','limegreen','yellow', 'orange', 'orangered', 'darkred']}, 'Mean chlor-a', false);
Map.addLayer(sd_coll.mean(), {min: 0, max: 2, palette: ['800000', 'FF9700', '7BFF7B', '0080FF', '000080']}, 'Mean Zsd', false);
Map.addLayer(tsi_coll.mean(), {min: 30, max: 80, palette: ['darkblue','blue','cyan','limegreen','yellow', 'orange', 'orangered', 'darkred']}, 'Mean TSI', false);
Map.addLayer(tsi_collR.mode(), {min: 1, max: 7, palette: ['purple', 'blue', 'limegreen', 'yellow', 'orange', 'orangered', 'darkred']}, 'Mode TSI Class', true);
