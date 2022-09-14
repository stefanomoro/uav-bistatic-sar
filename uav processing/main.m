clear variables
close all
clc

addpath(genpath('..\lib' ));       % add path of lib
addpath("./utils")

%% CONSTANTS
param = load("./test_parameters/20220826/test6").param;

[const] = initializeConstants(param);

%% PRE-PROCESSING
[radar,drone,targets] = loadData(const);
% Cut processed RC
[radar] = cutProcessedRC(radar,param);
plotRC(radar,[],1)
% Align drone time and radar time and interpolate
[RX,TX,radar] = alignInterpolateDroneRadarTime(const,drone,radar);
% Compute all useful distances
[scenario] = computeAllDistances(const,radar,RX,TX,targets);
%Rotate path
[RX,TX,targets,drone] = rotateFlightPath(RX,TX,targets,drone);

%% RC Synch, correct and prepare for Focusing
[radar]= interpolateRC(const,radar);
[radar] = computeCorrCrossTalk(const,radar,scenario);
[radar] = timeShiftCorrection(radar);
[radar] = getCrossTalk(radar);
[radar] = freqShiftCorrection(const,radar);
plotDroneTrack(drone,RX,TX,targets);
[radar] = filterGaussRC(const,radar,.75*const.dt);

plotRC(radar,scenario,1);

plotRC(radar,scenario,2);
plotRC(radar,scenario,3);


%% FOCUSING
[scenario] = defineFocusingGrid(const,scenario,RX);
[focus] = focusingTDBP_GPU(const,radar,scenario,RX,TX);
figure,plotFocusedWithTargets(scenario,RX,TX,targets,20*log10(focus.Focused_vec(:,:,8)),...
        "0° squint, not equalized");
caxis([100 160])
% EQUALIZE distance
[focus] = equalizeDistanceRC(scenario,TX,focus);

%% Convert cellarray to matrix
if(iscell(focus.Focused_vec))
    focus = convertCellToMatrix(focus);
end
%% SUM all the focused angles
focus_all_sum = sumAllFocusedAngle(focus.Focus_eq,1);

figure,plotFocusedWithTargets(scenario,RX,TX,targets,focus_all_sum,...
        "Sum all squint,linear");      
figure,plotFocusedWithTargets(scenario,RX,TX,targets,20*log10(focus_all_sum),...
    "Sum all squint, dB" );   
caxis([180 230]) 
figure,plotFocusedWithTargets(scenario,RX,TX,targets,...
    20*log10(filterHammingFocus(focus_all_sum,3)),...
        "Sum all squint,dB - Hamming filtered");      


caxis([180 230]) 

%% Show 
figure
for idx = 1:8
subplot(2,4,idx)
F = filterHammingFocus(focus.Focus_eq(:,:,idx),3);
% F = focus.Focused_vec(:,:,idx);
plotFocusedWithTargets(scenario,RX,TX,targets,20*log10(F),...
        strcat("Squint ",num2str(focus.angle_vec(idx)) ,"°, dB"));  
    caxis([100 200])
end
figure
for idx = 8:15
subplot(2,4,idx-7)
F = filterHammingFocus(focus.Focus_eq(:,:,idx),3);
% F = focus.Focused_vec(:,:,idx);
plotFocusedWithTargets(scenario,RX,TX,targets,20*log10(F),...
        strcat("Squint ",num2str(focus.angle_vec(idx)) ,"°, dB"));  
    caxis([100 200])
end

%%
makeGIF(const,scenario,RX,TX,targets,focus,2,[100 200]);
% saveAllSquintImages(const,scenario,RX,TX,targets,focus,2);
%% SAVE vec
save( getResultsFileName(focus),"param","const","scenario","focus")

