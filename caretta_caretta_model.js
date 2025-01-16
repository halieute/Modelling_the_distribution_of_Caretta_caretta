////////////////////////////////////////////////////////////////////////////////////////////////////////
// Modelling the distribution of Caretta caretta (Linnaeus, 1758) species in the North Atlantic "2000 - 2010"
////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load Area of Interest (North Atlantic region defined using a polygon geometry)
var northAtlantic = ee.Geometry.Polygon([
    [[-80, 0], [-80, 60], [10, 60], [10, 0], [-80, 0]]
  ]);
  
  
  Map.centerObject(northAtlantic, 3);
  
  // Define text annotations
  var title = ui.Label('Modelling the distribution of Caretta caretta (Linnaeus, 1758) species in the North Atlantic "2000 - 2010"');
  title.style().set({
    position: 'top-center',
    fontWeight: 'bold',
    fontSize: '24px',
    backgroundColor: 'rgba(255, 255, 255, 0.5)',
    margin: '500px 0 0 10' // Move title down by increasing the top margin
  });
  Map.add(title);
  
  var authorDate = ui.Label('Created by Souleymane Maman Nouri Souley, Date: ' + (new Date()).toDateString());
  authorDate.style().set({
    position: 'bottom-right',
    fontWeight: 'normal',
    fontSize: '12px',
    backgroundColor: 'rgba(255, 255, 255, 0.7)',
    textAlign: 'right'
  });
  Map.add(authorDate);
  
  // Define colors for legend
  var legendColors = [
    {label: 'Low Elevation', color: '#006600'}, // Corrected single #
    {label: 'High Elevation', color: '#CC9966'},
    {label: 'Low Presence', color: '#FFFFFF'}, // White
    {label: 'High Presence', color: '#008000'}, // Green
    {label: 'Low PC1', color: '#FFFFFF'}, // White
    {label: 'High PC1', color: '#FF0000'}, // Red
    {label: 'Low PC2', color: '#FFFFFF'}, // White
    {label: 'High PC2', color: '#0000FF'}  // Blue
  ];
  
  
  // Create legend panel
  var legend = ui.Panel({
    style: {
      position: 'bottom-left',
      padding: '8px 15px'
    }
  });
  legend.add(ui.Label({
    value: 'Legend',
    style: {fontWeight: 'bold'}
  }));
  
  // Add legend colors and labels
  legendColors.forEach(function(element) {
    var colorBox = ui.Label({
      style: {
        backgroundColor: '' + element.color,
        padding: '10px',
        margin: '4px'
      }
    });
    var description = ui.Label(element.label, {margin: '5px'});
    legend.add(ui.Panel([colorBox, description], ui.Panel.Layout.Flow('horizontal')));
  });
  Map.add(legend);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Section 1 - Data Loading and Preprocessing
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Import the species occurrence dataset
  var speciesData = ee.FeatureCollection("projects/ee-ssouleyniger/assets/new_dataset_date_num"); // Presence CSV file of Caretta caretta species in the North Atlantic Ocean
  
  // Define spatial resolution for analysis (in meters, where 10 km = 10000 meters)
  var grainSize = 10000; // 10 km
  
  // Function to remove spatial duplicates based on grid size
  function removeDuplicates(data) {
    var randomRaster = ee.Image.random().multiply(1000000).reproject('EPSG:4326', null, grainSize);
    var randPointVals = randomRaster.sampleRegions({
      collection: ee.FeatureCollection(data), 
      scale: grainSize, 
      geometries: true
    });
    return randPointVals.distinct('random');
  }
  
  // Apply duplicate removal and print results
  var processedData = removeDuplicates(speciesData);
  print('Number of unique occurrence points:', processedData.size());
  
  // Load Area of Interest (North Atlantic region defined using a polygon geometry)
  var northAtlantic = ee.Geometry.Polygon([
    [[-80, 0], [-80, 60], [10, 60], [10, 0], [-80, 0]]
  ]);
  
  Map.addLayer(northAtlantic, {}, 'North Atlantic', false);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Section 2 - Environmental Variables Processing
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // Load bioclimatic variables and terrain data
  var bioclimatic = ee.Image("WORLDCLIM/V1/BIO");
  var terrain = ee.Algorithms.Terrain(ee.Image("USGS/SRTMGL1_003"));
  
  // Calculate median tree cover from MODIS collection
  var modisCollection = ee.ImageCollection("MODIS/006/MOD44B");
  var medianTreeCover = modisCollection
    .filterDate('2000-01-01', '2010-12-31')
    .select(['Percent_Tree_Cover'])
    .median();
  
  // Combine all environmental predictors
  var environmentalLayers = bioclimatic.addBands(terrain).addBands(medianTreeCover);
  
  // Create and apply land mask
  var landMask = terrain.select('elevation').gt(0);
  environmentalLayers = environmentalLayers.updateMask(landMask).clip(northAtlantic);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Section 3 - Principal Component Analysis (PCA)
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // Function to generate new band names for PCA components
  var generatePCABandNames = function(prefix, count) {
    return ee.List.sequence(1, count).map(function(i) {
      return ee.String(prefix).cat(ee.Number(i).int());
    });
  };
  
  // Enhanced PCA calculation function
  var calculatePCA = function(centeredImage, scale, region) {
    var bands = centeredImage.bandNames().length();
    var arrays = centeredImage.toArray();
  
    var covariance = arrays.reduceRegion({
      reducer: ee.Reducer.centeredCovariance(),
      geometry: region,
      scale: scale,
      maxPixels: 1e9
    });
  
    var covarArray = ee.Array(covariance.get('array'));
    var eigens = covarArray.eigen();
  
    var eigenValues = eigens.slice(1, 0, 1);
    var eigenVectors = eigens.slice(1, 1);
  
    var arrayImage = arrays.toArray(1);
    var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage);
  
    var sdImage = ee.Image(eigenValues.sqrt())
      .arrayProject([0])
      .arrayFlatten([generatePCABandNames('sd', bands)]);
  
    return principalComponents
      .arrayProject([0])
      .arrayFlatten([generatePCABandNames('pc', bands)])
      .divide(sdImage);
  };
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Section 4 - PCA Implementation and Final Variable Selection
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // Select bioclimatic variables for PCA
  var bioVariables = bioclimatic.select([
    'bio01', 'bio02', 'bio03', 'bio04', 'bio05', 'bio06', 'bio07', 'bio08',
    'bio09', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16',
    'bio17', 'bio18', 'bio19'
  ]);
  
  // Center the bioclimatic variables
  var meanStats = bioVariables.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: processedData,
    scale: 1000,
    maxPixels: 1e9
  });
  
  var meanImage = ee.Image.constant(ee.List(meanStats.values(bioVariables.bandNames())));
  var centeredVariables = bioVariables.subtract(meanImage);
  
  // Calculate PCA and select components
  var pcaResults = calculatePCA(centeredVariables, grainSize, processedData);
  var selectedPCs = pcaResults.select(['pc1', 'pc2']);
  
  // Print PCA results
  print('Selected Principal Components (PC1 and PC2):', selectedPCs);
  
  // Combine final predictor variables
  var finalPredictors = environmentalLayers
    .addBands(selectedPCs)
    .select(['pc1', 'pc2', 'elevation', 'Percent_Tree_Cover']);
  
  // Print final predictor information
  print('Final SDM Predictors:', finalPredictors);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Section 5 - Visualization and Correlation Analysis
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // Add layers to the map
  Map.addLayer(finalPredictors, {}, 'All Predictors', false);
  Map.addLayer(finalPredictors, {
    bands: ['elevation'],
    min: 0,
    max: 5000,
    palette: ['000000', '006600', '009900', '33CC00', '996600', 'CC9900', 'CC9966', 'FFFFFF']
  }, 'Elevation (m)', false);
  Map.addLayer(selectedPCs, {
    bands: ['pc1'],
    min: 190,
    max: 400,
    palette: 'white,red'
  }, 'PC1', false);
  Map.addLayer(selectedPCs, {
    bands: ['pc2'],
    min: 0,
    max: 4000,
    palette: 'white,blue'
  }, 'PC2', false);
  Map.addLayer(finalPredictors, {
    bands: ['Percent_Tree_Cover'],
    min: 1,
    max: 100,
    palette: 'white,yellow,green'
  }, 'Tree Cover (%)', false);
  
  // Calculate Spearman correlation between predictors
  var pixelValues = finalPredictors.sampleRegions({
    collection: processedData,
    scale: grainSize,
    tileScale: 16
  }).filter(ee.Filter.notNull(finalPredictors.bandNames()));
  
  var correlationMatrix = ee.List(finalPredictors.bandNames()).map(function(band1) {
    return ee.List(finalPredictors.bandNames()).map(function(band2) {
      var correlation = pixelValues.reduceColumns({
        reducer: ee.Reducer.spearmansCorrelation(),
        selectors: [band1, band2]
      }).get('correlation');
      return correlation;
    })
  });
  
  print('Predictor Variable Correlation Matrix:', correlationMatrix);
  