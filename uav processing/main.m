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

focus_all_sum = sumAllFocusedAngle(focus.Focus_eq);

figure,plotFocusedWithTargets(scenario,RX,TX,targets,10*log10(focus.Focused_vec{1}./focus.not_coh_sum{1}),...
        "sum all angles" );        
%%
A = conv2(focus_all_sum,hamming(3),'same')./3;
At = conv2(A,hamming(3).','same')./3;
figure,imagesc(scenario.grid.x_ax,scenario.grid.y_ax,10*log10(At.'));colorbar,caxis([190 240]),axis('xy')
%%
idx = 15;
figure,plotFocusedWithTargets(scenario,RX,TX,targets,10*log10(focus.Focus_eq{idx}),...
        strcat("single image ",num2str(focus.angle_vec(idx) )));        
caxis

%%
figure
for idx = 8:13
subplot(2,3,idx-7)

imagesc(scenario.grid.x_ax,scenario.grid.y_ax,abs(10*log10(focus.Focus_eq{idx}).')), axis xy , 
%         title(title_txt)
title(focus.angle_vec(idx))
        xlabel('[m]'), ylabel('[m]')
        hold on,
        plot3(RX.pos(1,1),RX.pos(2,1),RX.pos(3,1),'ro'), 
        plot3(RX.pos(1,end),RX.pos(2,end),RX.pos(3,end),'go'),
        plot3(TX.pos(1,1),TX.pos(2,1),TX.pos(3,1),'kd'), 
        plot3(targets.cars(1,:),targets.cars(2,:),targets.cars(3,:),'rp'), 
        plot3(targets.humans(1,:),targets.humans(2,:),targets.humans(3,:),'yh')
        legend('start drone track','end drone track','TX','Cars','Humans');
        hold off
        colorbar
        caxis([90 120])
end


%% Detect targets
Focused_vec_norm = cell(size(Focused_vec));

for i = 1:max(size(Focused_vec))
    Focused_vec_norm{i} = abs(Focused_vec{i}) ./ max(abs(Focused_vec{i}));
end
%%
makeGIF(const,scenario,RX,TX,targets,focus,2);
%% SAVE vec
save( getResultsFileName(focus),"const","scenario","focus")

