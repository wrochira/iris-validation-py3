let selectedVersion = 1;
let selectedChain = 0;
let selectedResidue = 0;

let residueSummary = null;
let switchMovementAnimations = [ ];
let switchColorAnimations = [ ];

let chainViews = [ ];
let residueSelectors = [ ];
let interactionSegmentSets = [ ];
let chainSelectors = [ ];
let shadeGroups = [ ];
let discreteGroupSets = [ ];
let lineAnimationSets = [ ];
let residueSelectorDragging = false;

let residueView = null;
let barChartsContainer = null;
let boxes = [ ];
let boxTexts = [ ];
let boxplots = [ ];
let barMainlines = [ ];
let barLabels = [ ];
let boxplotBoxes = [ ];
let boxplotLines = [ ];
let barLineYs = [ ];
let barOffsetY = null;
let barMultiplierY = null;


//
// Utility functions
//
function mean(values) {
  let sum = values.reduce(function(sum, value) {
    return sum + value;
  }, 0);
  let avg = sum / values.length;
  return avg;
};


function standardDeviation(values) {
  avg = mean(values);
  let squareDiffs = values.map(function(value) {
    let diff = value - avg;
    let sqrDiff = diff * diff;
    return sqrDiff;
  });
  let variance = mean(squareDiffs);
  let stdDev = Math.sqrt(variance);
  return stdDev;
};


function percentile(values, p) {
  values = values.sort();
  let pos = ((values.length) - 1) * p;
  let base = Math.floor(pos);
  let rest = pos - base;
  if( (values[base+1] !== undefined) ) {
    return values[base] + rest * (values[base+1] - values[base]);
  } else {
    return values[base];
  };
};


function coordsFromAngle(centre, angle, pointRadius) {
  let xc = centre[0];
  let yc = centre[1];
  let x1 = xc + pointRadius * Math.sin(angle);
  let y1 = yc - pointRadius * Math.cos(angle);
  return [Math.round(x1, 2), Math.round(y1, 2)];
};


function setPointer() {
    document.body.style.cursor = 'pointer';
};


function unsetPointer() {
    document.body.style.cursor = '';
};


//
// One-off functions
//
function getResidueViewData() {
  // Dimensions
  let bccPoints = [ ];
  for (var pointID = 0; pointID < barChartsContainer.points.numberOfItems; ++pointID) {
    let x = barChartsContainer.points.getItem(pointID).x;
    let y = barChartsContainer.points.getItem(pointID).y;
    bccPoints.push([x, y]);
  };
  barOffsetY = bccPoints[2][1];
  barMultiplierY = -(bccPoints[2][1]-bccPoints[0][1]) / 100;

  // Boxplot ranges
  for (var versionID = 0; versionID < modelData[selectedChain]['num_versions']; ++versionID) {
    barLineYs.push([ ]); // Model-version holder
    for (var barID = 0; barID < barMetricIDs.length; ++barID) {
      let allPercentileValues = [ ];
      let metricID = barMetricIDs[barID];
      for (var chainID = 0; chainID < numChains; ++chainID) {
        for (var residueID = 0; residueID < modelData[chainID]['aligned_length']; ++residueID) {
          if (modelData[chainID]['residue_validities'][versionID][residueID]) {
            let percentileValue = modelData[chainID]['percentile_values'][metricID][versionID][residueID];
            if (percentileValue !== null) {
              allPercentileValues.push(percentileValue);
            };
          };
        };
      };
      let metricMin = Math.min.apply(null, allPercentileValues);
      let metricMax = Math.max.apply(null, allPercentileValues);
      let metricMean = mean(allPercentileValues);
      let metricStd = standardDeviation(allPercentileValues);
      let metricLow = Math.max(0, metricMean-metricStd);
      let metricHigh = Math.min(100, metricMean+metricStd);
      let distributionValues = [ metricMin, metricMax, metricLow, metricMean, metricHigh ];
      let versionLineYs = [ ];
      for (var valueID = 0; valueID < 5; ++valueID) {
        let lineY = parseFloat((barOffsetY + barMultiplierY * distributionValues[valueID]).toFixed(1));
        versionLineYs.push(lineY);
      };
      barLineYs[versionID].push(versionLineYs);
    };
  };
};


//
// Interaction functions
//
function toggleVersion() {
  selectedVersion = (selectedVersion + 1) % modelData[selectedChain]['num_versions'];
  switchMovementAnimations[selectedVersion].beginElement();
  switchColorAnimations[selectedVersion].beginElement();
  updateSelectedVersion();
};


function setChain(chainID) {
  selectedChain = chainID;
  selectedResidue = 0;
  residueSelectorDragging = false;
  updateSelectedChain();
};


function handleSegment(actionID, residueID) {
  if (actionID === 1 || (actionID === 2 && residueSelectorDragging)) {
    residueSelectorDragging = true;
    selectedResidue = residueID;
    if (modelData[selectedChain]['residue_validities'][selectedVersion][selectedResidue]) {
      updateSelectedResidue();
    };
  } else if (actionID === 3) {
    residueSelectorDragging = false;
  };
};


function updateSelectedVersion() {
  // Chain view
  for (var chainID = 0; chainID < numChains; ++chainID) {
    for (var versionID = 0; versionID < modelData[selectedChain]['num_versions']; ++versionID) {
      let opacity = versionID === selectedVersion ? 1 : 0;
      shadeGroups[chainID][versionID].setAttribute('opacity', opacity);
      let discreteGroups = discreteGroupSets[chainID][versionID];
      for (var groupID = 0; groupID < discreteGroups.length; ++groupID) {
        discreteGroups[groupID].setAttribute('opacity', opacity);
      };
      if (versionID === selectedVersion) {
        let lineAnimations = lineAnimationSets[chainID][versionID];
        for (var animationID = 0; animationID < lineAnimations.length; ++animationID) {
          lineAnimations[animationID].beginElement();
        };
      };
    };
  };

  // Residue view
  for (var barID = 0; barID < barMetricIDs.length; ++barID) {
    // Update boxplots
    let boxplotPoints = [ ];
    for (var pointID = 0; pointID < boxplotBoxes[barID].points.numberOfItems; ++pointID) {
      let x = boxplotBoxes[barID].points.getItem(pointID).x;
      let y = boxplotBoxes[barID].points.getItem(pointID).y;
      boxplotPoints.push([x, y]);
    };
    boxplotPoints[1][1] = boxplotPoints[2][1] = barLineYs[selectedVersion][barID][0];
    boxplotPoints[0][1] = boxplotPoints[3][1] = barLineYs[selectedVersion][barID][1];
    boxplotBoxes[barID].setAttribute('points', boxplotPoints);

    // Update bar lines
    boxplotLines[barID][0].setAttribute('y1', barLineYs[selectedVersion][barID][2]);
    boxplotLines[barID][0].setAttribute('y2', barLineYs[selectedVersion][barID][2]);
    boxplotLines[barID][1].setAttribute('y1', barLineYs[selectedVersion][barID][3]);
    boxplotLines[barID][1].setAttribute('y2', barLineYs[selectedVersion][barID][3]);
    boxplotLines[barID][2].setAttribute('y1', barLineYs[selectedVersion][barID][4]);
    boxplotLines[barID][2].setAttribute('y2', barLineYs[selectedVersion][barID][4]);

    // Only show SD lines if they fall within the min/max range
    boxplotLines[barID][0].setAttribute('opacity', 0);
    if (barLineYs[selectedVersion][barID][2] < barLineYs[selectedVersion][barID][0]) {
      boxplotLines[barID][0].setAttribute('opacity', 1);
    };
    boxplotLines[barID][2].setAttribute('opacity', 0);
    if (barLineYs[selectedVersion][barID][4] > barLineYs[selectedVersion][barID][1]) {
      boxplotLines[barID][2].setAttribute('opacity', 1);
    };
  };

  updateSelectedChain();
};


function updateSelectedChain() {
  // Chain view
  for (var chainID = 0; chainID < modelData.length; ++chainID) {
    if (chainID === selectedChain) {
      chainViews[chainID].style.display = '';
      chainSelectors[chainID].setAttribute('fill', chainSelectorColors[1]);
    } else {
      chainViews[chainID].style.display = 'none';
      chainSelectors[chainID].setAttribute('fill', chainSelectorColors[0]);
    };
  };

  updateSelectedResidue();
};


function updateSelectedResidue() {
  // If selected residue is null, cycle residues
  while (!modelData[selectedChain]['residue_validities'][selectedVersion][selectedResidue]) {
    selectedResidue = (selectedResidue + 1) % modelData[selectedChain]['aligned_length'];
  };

  // Chain view 
  residueSelectors[selectedChain].setAttribute('transform', 'rotate(' + (360-gapDegrees)/modelData[selectedChain]['aligned_length']*selectedResidue + ', 500, 500' + ')');
  for (var residueID = 0; residueID < interactionSegmentSets[selectedChain].length; ++residueID) {
    if (residueID === selectedResidue) {
      interactionSegmentSets[selectedChain][residueID].setAttribute('stroke-opacity', 0.25);
      interactionSegmentSets[selectedChain][residueID].setAttribute('fill-opacity', 0.25);
    } else {
      interactionSegmentSets[selectedChain][residueID].setAttribute('stroke-opacity', 0);
      interactionSegmentSets[selectedChain][residueID].setAttribute('fill-opacity', 0);
    };
  };

  // Residue view
  // Set discrete boxes
  for (var boxID = 0; boxID < boxMetricIDs.length; ++boxID) {
    let metricID = boxMetricIDs[boxID];
    let discreteIndex = modelData[selectedChain]['discrete_values'][metricID][selectedVersion][selectedResidue];
    let boxColor = 'rgb(200, 200, 200)';
    let boxText = 'N/A';
    if (discreteIndex !== null) {
      boxColor = boxColors[boxID][discreteIndex];
      boxText = boxLabels[boxID][discreteIndex];
    };
    boxes[boxID].setAttribute('fill', boxColor);
    boxTexts[boxID].textContent = boxText;
  };

  // Make boxplots visible
  for (var barID = 0; barID < barMetricIDs.length; ++barID) {
    let metricID = barMetricIDs[barID];
    // Get bar value
    boxplots[barID].setAttribute('opacity', 1);
    let barValue = modelData[selectedChain]['percentile_values'][metricID][selectedVersion][selectedResidue];
    // Set main line coordinates
    barY = parseFloat((barOffsetY + barMultiplierY * barValue).toFixed(1));
    barMainlines[barID].setAttribute('y1', barY);
    barMainlines[barID].setAttribute('y2', barY);
    // Set bar label text and position
    barLabels[barID].textContent = barValue;
    if (barValue < 10) {
      barLabels[barID].setAttribute('y', barY-10);
    } else {
      barLabels[barID].setAttribute('y', barY+25);
    };
  };

  // Set summary text
  let seqNum = modelData[selectedChain]['residue_seqnos'][selectedVersion][selectedResidue];
  let aaCode = modelData[selectedChain]['residue_codes'][selectedVersion][selectedResidue];
  residueSummary.textContent = 'Residue ' + seqNum + ' (' + aaCode + ')';
};


function loadElements() {
  // Panel
  residueSummary = document.getElementById('iris-panel-residue-summary');
  switchMovementAnimations = [ document.getElementById('iris-panel-switch-move-animation-0'),
                               document.getElementById('iris-panel-switch-move-animation-1') ]
  switchColorAnimations = [ document.getElementById('iris-panel-switch-color-animation-0'),
                            document.getElementById('iris-panel-switch-color-animation-1') ]
  for (var chainID = 0; chainID < modelData.length; ++chainID) {
    let chainSelector = document.getElementById('iris-panel-chain-selector-' + chainID);
    chainSelectors.push(chainSelector);
  };

  // Chain view
  for (var chainID = 0; chainID < modelData.length; ++chainID) {
    let chainViewID = 'iris-chain-view-' + chainID;
    let chainView = document.getElementById(chainViewID);
    chainViews.push(chainView);
    let residueSelector = document.getElementById(chainViewID + '-residue-selector');
    residueSelectors.push(residueSelector);
    let interactionSegmentSet = document.querySelectorAll('[id^=' + chainViewID + '-interaction-segment-]');
    interactionSegmentSets.push(interactionSegmentSet);
    shadeGroups.push([ ]);
    discreteGroupSets.push([ ]);
    lineAnimationSets.push([ ]);
    for (var versionID = 0; versionID < modelData[selectedChain]['num_versions']; ++versionID) {
      let shadeGroup = document.getElementById(chainViewID + '-shade-' + versionID);
      shadeGroups[chainID].push(shadeGroup);
      let discreteGroupSet = document.querySelectorAll('[id^=' + chainViewID + '-discrete-' +  versionID + '-]');
      discreteGroupSets[chainID].push(discreteGroupSet);
      let lineAnimationSet = document.querySelectorAll('[id^=' + chainViewID + '-animation-' +  versionID + '-]');
      lineAnimationSets[chainID].push(lineAnimationSet);
    };
  };

  // Residue view
  residueView = document.getElementById('iris-residue-view');
  barChartsContainer = document.getElementById('iris-residue-view-bar-charts-container');
  for (var boxID = 0; boxID < boxMetricIDs.length; ++boxID) {
    box = document.getElementById('iris-residue-view-box-' + boxID);
    boxes.push(box);
    boxText = document.getElementById('iris-residue-view-box-' + boxID + '-text');
    boxTexts.push(boxText);
  };
  for (var barID = 0; barID < barMetricIDs.length; ++barID) {
    boxplot = document.getElementById('iris-residue-view-boxplot-' + barID);
    boxplots.push(boxplot);
    barMainline = document.getElementById('iris-residue-view-bar-' + barID + '-mainline');
    barMainlines.push(barMainline);
    barLabel = document.getElementById('iris-residue-view-bar-' + barID + '-label');
    barLabels.push(barLabel);
    boxplotBox = document.getElementById('iris-residue-view-boxplot-' + barID + '-box');
    boxplotBoxes.push(boxplotBox);
    let boxplotLineLow = document.getElementById('iris-residue-view-boxplot-' + barID + '-line-low');
    let boxplotLineMid = document.getElementById('iris-residue-view-boxplot-' + barID + '-line-mid');
    let boxplotLineHigh = document.getElementById('iris-residue-view-boxplot-' + barID + '-line-high');
    boxplotLines.push([ boxplotLineLow, boxplotLineMid, boxplotLineHigh ]);
  };
};


function main() {
  loadElements();
  getResidueViewData();
  updateSelectedVersion();
};

window.addEventListener('load', main);
