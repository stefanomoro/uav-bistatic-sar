clear variables
close all
clc
%% =============================================================================  LOAD RC MATRIX CUT IN TARGET ZONE
%% CONSTANTS
[const] = initializeConstants();

%% PRE-PROCESSING
addpath(genpath('..\lib' ));       % add path of lib
[radar,drone,targets] = loadData(const);
% Cut processed RC
[radar] = cutProcessedRC(radar);
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

% %% SUM all the focused angles
% 
% if exist("Focused_vec","var")
%     focus_all_sum = sumAllFocusedAngle(Focused_vec);
%     focus_all_sum = equalizeDistanceRC(focus_all_sum,x_ax,y_ax,z0,TX_pos,RX_pos,0);
%     figure,plotWithTargets(x_ax,y_ax,focus_all_sum,"SUM all focus angles",RX_pos,TX_pos,cars,humans)
%         
% end    
% %% Detect targets
% Focused_vec_norm = cell(size(Focused_vec));
% 
% for i = 1:max(size(Focused_vec))
%     Focused_vec_norm{i} = abs(Focused_vec{i}) ./ max(abs(Focused_vec{i}));
% end

makeGIF(const,scenario,RX,TX,targets,focus);
%% SAVE vec
save( getResultsFileName(focus),"const","scenario","focus")

