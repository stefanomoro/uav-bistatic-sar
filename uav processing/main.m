clear all
close all
clc
%% =============================================================================  LOAD RC MATRIX CUT IN TARGET ZONE
%% COSTANTS
chirp_bw = 56e6*0.9;                                % actual chirp bandwidth (B)
chirp_sr = 56e6;                                % SDR sample rate (fs)
norm_B = chirp_bw / chirp_sr;

experiment_name = "test1_cut1";

radar_folder_name = '..\mat_files\uav_test\20220502\radar\test1\';
drone_track_folder = '..\mat_files\uav_test\20220502\drone track\';

drone_data = load(strcat(drone_track_folder, 'flight1')).drone;
radar_last_mod_date = load(strcat(radar_folder_name,'test1_last_mod_date')).last_mod_date;
targets = load(strcat(drone_track_folder,'target.mat'));

tx_wave = load(strcat('..\tx_waveform/tx_waveform_S56M.mat')).s_pad;
tx_wave = single(tx_wave);
samples_per_chirp = length(tx_wave);
f0 = 1.65e9;
lambda = physconst('LightSpeed')/f0;
PRI = samples_per_chirp / chirp_sr;
PRF = 1 / PRI;
dt = 1/chirp_sr;
dR = (physconst('LightSpeed') * dt);

%% RC Load
addpath(genpath('..\lib' ));       % add path of lib

RC_file = strcat(radar_folder_name,experiment_name);

RC = load_bin(RC_file);
R_margin = size(RC,1)*dR;
R_ax = -R_margin/2:dR:R_margin/2;
tau_ax = 0:PRI:PRI*size(RC,2);
N_PRI_tot = length(tau_ax);

%% taglio sul cross talk
idx_wind_row = 150:300;%180:230;
idx_wind_col = 1:2.4e4;%1:1.4e4;

RC = RC(idx_wind_row,idx_wind_col);
R_ax = R_ax(idx_wind_row); 
tau_ax = tau_ax(idx_wind_col); 

N_PRI = length(tau_ax);
%%
figure
imagesc(abs(RC))
title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")
%% Plot
figure
imagesc(tau_ax,R_ax,abs(RC))
title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

figure
imagesc(tau_ax,R_ax,angle(RC))
title('RC - angle plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

%%  ======================================================================= SYNCHRONISM ERROR CORRECTION
% Distance tx-rx for each PRI (cross-talk distance)
run('DronePath')

%%
% Synchronism correction algorythm
run('RC_SynchrCorrection.m')

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

%% ===========================================================================  CROSSRTALK REMOVAL - BRUTAL
% cross_talk_size_idxs = 10;
% 
% for n = 1:N_PRI
%     RC(1:cross_talk_idxs_corr(n)+cross_talk_size_idxs,n) = zeros(cross_talk_idxs_corr(n)+cross_talk_size_idxs,1);
% end
% 
% %% Plot
% figure
% imagesc(tau_ax,R_ax,abs(RC))
% title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")
% 
% figure
% imagesc(tau_ax,R_ax,angle(RC))
% title('RC - angle plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

%% =========================================================================== FOCUSING
%run('Focusing_WnAngle.m');

run('Focusing_WnWaveNumb.m');
%%
if exist("Focused_vec","var")
    fig = figure;
    images = cell(size(Focused_vec));
    for i = 1:length(angle_vec)
        imagesc(x_ax,y_ax,abs(Focused_vec{i})), axis xy , 
        title(strcat("Focused image with angle ",num2str(angle_vec(i)),"°" ))
        xlabel('[m]'), ylabel('[m]')

        hold on,
        plot3(RX_pos(1,1),RX_pos(2,1),RX_pos(3,1),'ro'), 
        plot3(RX_pos(1,end),RX_pos(2,end),RX_pos(3,end),'go'),
        plot3(TX_pos(1,1),TX_pos(2,1),TX_pos(3,1),'kd'), plot3(cars(1,:),cars(2,:),cars(3,:),'rp'), plot3(humans(1,:),humans(2,:),humans(3,:),'yh')
        legend('start drone track','end drone track','TX','Cars','Humans');
        hold off
        frame = getframe(fig);
        images{i} = frame2im(frame);
    end
% MAKE GIF
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

return
%%
figure,
imagesc(x_ax,y_ax,abs(Focus)), axis xy , 
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
imagesc(x_ax,y_ax,10*log10(abs(Focus))), axis xy , %colormap jet
%title('Initial Sensor Position') 
title('Focused image by TDBP - dB')
xlabel('azimuth [m]'), ylabel('ground range [m]')

hold on,
plot3(drone_start_pos(1),drone_start_pos(2),drone_start_pos(3),'ro'), 
plot3(drone_end_pos(1),drone_end_pos(2),drone_end_pos(3),'go'),
plot3(TX_pos(1,1),TX_pos(2,1),TX_pos(3,1),'kd'), plot3(cars(1,:),cars(2,:),cars(3,:),'rp'), plot3(humans(1,:),humans(2,:),humans(3,:),'yh')
legend('start swath','end swath','TX','Cars','Humans');
