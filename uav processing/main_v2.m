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
radar = cutProcessedRC(radar);
plotGenericRC(radar)
% Align drone time and radar time and interpolate
[RX,TX,radar] = alignInterpolateDroneRadarTime(const,drone,radar);
plotDroneTrack(drone,RX,TX,targets);
% Compute all useful distances
scenario = computeAllDistances(const,radar,RX,TX,targets);
%Rotate path
[RX,TX,targets,drone] = rotateFlightPath(RX,TX,targets,drone);
plotDroneTrack(drone,RX,TX,targets);
%% RC Synch and correction
% Interpolate RC
radar = interpolateRC(const,radar);
% Compute correct cross talk position
[radar] = computeCorrCrossTalk(const,radar,scenario);
[radar] = timeShiftCorrection(radar);
[radar] = getCrossTalk(radar);
radar = freqShiftCorrection(radar);
%% Synchronism correction algorythm
run('RC_SynchrCorrection.m')

%% GAUSS filter RC
RC = filterGaussRC(RC,dt,OSF,.75*dt);

%% Plot
figure
imagesc(tau_ax,R_ax,abs(RC))
title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")
hold on
plot(tau_ax,distance_tx_rx,'r','LineWidth',1.5),
plot(tau_ax,distance_cars(1,:),'g')
legend('crosstalk','car')%,'car 2','human1','human2','human3');

%% Plot
figure
imagesc(tau_ax,R_ax,abs(RC))
title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")
hold on

plot(tau_ax,distance_tx_rx),
plot(tau_ax,distance_cars(1,:),'m'), plot(tau_ax,distance_cars(2,:),'g'), 
plot(tau_ax,distance_humans(1,:),'y'), plot(tau_ax,distance_humans(2,:),'m'), 
plot(tau_ax,distance_humans(3,:),'b')
legend('crosstalk','car','car 2','human1','human2','human3');

%% FOCUSING
run('Focusing_WnWaveNumb.m');

%% EQUALIZE distance
Focus_eq = equalizeDistanceRC(Focus,x_ax,y_ax,z0,TX_pos,RX_pos,psi_foc);

%% SUM all the focused angles

if exist("Focused_vec","var")
    focus_all_sum = sumAllFocusedAngle(Focused_vec);
    focus_all_sum = equalizeDistanceRC(focus_all_sum,x_ax,y_ax,z0,TX_pos,RX_pos,0);
    figure,plotWithTargets(x_ax,y_ax,focus_all_sum,"SUM all focus angles",RX_pos,TX_pos,cars,humans)
        
end    
%% Detect targets
Focused_vec_norm = cell(size(Focused_vec));

for i = 1:max(size(Focused_vec))
    Focused_vec_norm{i} = abs(Focused_vec{i}) ./ max(abs(Focused_vec{i}));
end


%%
if exist("Focused_vec","var")
    fig = figure("WindowState","maximized");
    images = cell(size(Focused_vec));
    for i = 1:length(angle_vec)
        F = equalizeDistanceRC(Focused_vec{i},x_ax,y_ax,z0,TX_pos,RX_pos,angle_vec(i));        
        plotWithTargets(x_ax,y_ax,F,strcat("Focused image with angle ",num2str(angle_vec(i)),"°" ),...
            RX_pos,TX_pos,cars,humans);
        frame = getframe(fig);
        images{i} = frame2im(frame);
    end
%% MAKE GIF
    filename = strcat(experiment_name,'.gif'); % Specify the output file name
    for idx = 1:length(Focused_vec)
        [A,map] = rgb2ind(images{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.5);
        end
    end
end
%% SAVE vec
if exist("Focused_vec","var")
    i = 0;
    fname = strcat("focused_35_",num2str(i));
    while isfile(strcat("focused_result/",fname))
        i = i +1;
    end
    save(strcat("focused_result/",fname),"Focused_vec")
end

return

%%
figure,
imagesc(x_ax,y_ax,abs(Focus.')), axis xy , 
%title('Final Sensor Position') 
title('Focused image by TDBP ')
xlabel('azimuth [m]'), ylabel('ground range [m]')

hold on,
plot3(drone_start_pos(1),drone_start_pos(2),drone_start_pos(3),'ro'), 
plot3(drone_end_pos(1),drone_end_pos(2),drone_end_pos(3),'go'),
plot3(TX_pos(1,1),TX_pos(2,1),TX_pos(3,1),'kd'), plot3(cars(1,:),cars(2,:),cars(3,:),'rp'), plot3(humans(1,:),humans(2,:),humans(3,:),'yh')
legend('start drone track','end drone track','TX','Cars','Humans');
%%
figure,
imagesc(x_ax,y_ax,abs(Focus_eq.')), axis xy , 
%title('Final Sensor Position') 
title('Equalized Image')
xlabel('azimuth [m]'), ylabel('ground range [m]')

hold on,
plot3(drone_start_pos(1),drone_start_pos(2),drone_start_pos(3),'ro'), 
plot3(drone_end_pos(1),drone_end_pos(2),drone_end_pos(3),'go'),
plot3(TX_pos(1,1),TX_pos(2,1),TX_pos(3,1),'kd'), plot3(cars(1,:),cars(2,:),cars(3,:),'rp'), plot3(humans(1,:),humans(2,:),humans(3,:),'yh')
legend('start drone track','end drone track','TX','Cars','Humans');

%%
figure,
imagesc(x_ax,y_ax,10*log10(abs(Focus.'))), axis xy , %colormap jet
%title('Initial Sensor Position') 
title('Focused image by TDBP - dB')
xlabel('azimuth [m]'), ylabel('ground range [m]')

hold on,
plot3(drone_start_pos(1),drone_start_pos(2),drone_start_pos(3),'ro'), 
plot3(drone_end_pos(1),drone_end_pos(2),drone_end_pos(3),'go'),
plot3(TX_pos(1,1),TX_pos(2,1),TX_pos(3,1),'kd'), plot3(cars(1,:),cars(2,:),cars(3,:),'rp'), plot3(humans(1,:),humans(2,:),humans(3,:),'yh')
legend('start swath','end swath','TX','Cars','Humans');
