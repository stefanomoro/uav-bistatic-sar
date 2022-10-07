clear variables
close all
clc

addpath(genpath('../lib' ));       % add path of lib
addpath("./utils")

%% CONSTANTS
param = load("./test_parameters/20220826/test4").param;
[const] = initializeConstants(param);

%% PRE-PROCESSING
[radar,drone,targets] = loadData(const);
% Cut RC to process only a subset
[radar] = cutProcessedRC(radar,param);
plotRC(radar,[],1)
% Align drone time and radar time and interpolate
[RX,TX,radar] = alignInterpolateDroneRadarTime(const,drone,radar);
% Compute all useful distances
[scenario] = computeAllDistances(const,radar,RX,TX,targets);
%Rotate path to keep orthogonal to x axis
[RX,TX,targets,drone] = rotateFlightPath(RX,TX,targets,drone);

%% RC Synch, correct and prepare for Focusing

[radar]= interpolateRC(const,radar);
[radar] = computeCorrCrossTalk(const,radar,scenario);
[radar] = timeShiftCorrection(radar);
% Now radar.RC is interpolated and time shifted corrected
[radar] = getCrossTalk(radar);
[radar] = freqShiftCorrection(const,radar);
% Now radar.RC is also phase corrected with navigation computed phase 

% Apply gaussian filter to remove sidelobe of sinc in range
[radar] = filterGaussRC(const,radar,.75*const.dt);

% PLOTTING
% plotDroneTrack(drone,RX,TX,targets);
% plotRC(radar,scenario,1);
% plotRC(radar,scenario,2);
% plotRC(radar,scenario,3);


%% FOCUSING
[scenario] = defineFocusingGrid(const,scenario,RX);
[focus] = focusingTDBP_CUDA(const,radar,scenario,RX,TX);
% EQUALIZE distance
[focus] = equalizeDistanceRC(scenario,TX,focus);

%% Convert cellarray to matrix (backward compatibility)
if(iscell(focus.Focused_vec))
    focus = convertCellToMatrix(focus);
end
%% SAVE Focus result
save( getResultsFileName(focus),"param","const","scenario","focus")

%% PLOTS
figure,plotFocusedWithTargets(scenario,RX,TX,targets,20*log10(abs(focus.Focused_vec(:,:,1))),...
        "0° squint, not equalized");
caxis([100 160])

plotAllFocusedSquints(focus,scenario,RX,TX,targets,1,[100 160])
makeGIF(const,scenario,RX,TX,targets,focus,1,[120 160]);
% saveAllSquintImages(const,scenario,RX,TX,targets,focus,2);


%% SUM all the focused angles
focus_all_sum = sumAllFocusedAngle(focus.Focused_vec,1:size(focus.Focused_vec,3));
figure,plotFocusedWithTargets(scenario,RX,TX,targets,focus_all_sum,...
        "Sum all squint,linear");      
figure,plotFocusedWithTargets(scenario,RX,TX,targets,20*log10(focus_all_sum),...
    "Sum all squint, dB" );   
caxis([150 190]) 

figure,plotFocusedWithTargets(scenario,RX,TX,targets,...
    20*log10(filterHammingFocus(focus_all_sum,3)),...
        "Sum all squint,dB - Hamming filtered");      
caxis([140 170]) 

