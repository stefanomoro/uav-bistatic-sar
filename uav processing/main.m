clear variables
close all
clc
%% =============================================================================  LOAD RC MATRIX CUT IN TARGET ZONE
%% CONSTANTS
param = load("./test_parameters/20220502/test1_pass2_param").param;
[const] = initializeConstants(param);

%% PRE-PROCESSING
addpath(genpath('..\lib' ));       % add path of lib
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
plotDroneTrack(drone,RX,TX,targets);

%% RC Synch, correct and prepare for Focusing
[radar]= interpolateRC(const,radar);
[radar] = computeCorrCrossTalk(const,radar,scenario);
[radar] = timeShiftCorrection(radar);
[radar] = getCrossTalk(radar);
[radar] = freqShiftCorrection(const,radar);
[radar] = filterGaussRC(const,radar,.75*const.dt);

plotRC(radar,scenario,2);
plotRC(radar,scenario,3);


%% FOCUSING
[scenario] = defineFocusingGrid(const,scenario);
[focus] = focusingTDBP(const,radar,scenario,RX,TX);
% EQUALIZE distance
[focus] = equalizeDistanceRC(scenario,TX,focus);

%% SUM all the focused angles

focus_all_sum = sumAllFocusedAngle(focus.Focus_eq,1:15);

figure,plotFocusedWithTargets(scenario,RX,TX,targets,focus_all_sum,...
        "Sum all squint,linear");      
figure,plotFocusedWithTargets(scenario,RX,TX,targets,20*log10(focus_all_sum),...
    "Sum all squint, dB" );   
caxis([180 230]) 
figure,plotFocusedWithTargets(scenario,RX,TX,targets,...
    20*log10(filterHammingFocus(focus_all_sum,3)),...
        "Sum all squint,dB - Hamming filtered");      


caxis([180 230]) 

%% save subplot squint
figure
for idx = 3:8
subplot(2,3,idx-2)
% F = filterHammingFocus(focus.Focus_eq{idx},3);
F = focus.Focus_eq{idx};
plotFocusedWithTargets(scenario,RX,TX,targets,20*log10(F),...
        strcat("Squint ",num2str(focus.angle_vec(idx)) ,"°, dB"));  
    caxis([160 220])
end

%%
makeGIF(const,scenario,RX,TX,targets,focus,2);
saveAllSquintImages(const,scenario,RX,TX,targets,focus,2);
%% SAVE vec
save( getResultsFileName(focus),"param","const","scenario","focus")

