clear all
close all
clc
%% COSTANTS
chirp_sr = 56e6;                                % SDR sample rate
chirp_bw = .9*chirp_sr;                         % actual chirp bandwidth 

experiment_name = 'test';
folder_name = 'mat_files/uav/20220502/RC/cut/';

tx_wave = load(strcat('tx_waveform/tx_waveform_S56M.mat')).s_pad;
tx_wave = single(tx_wave);
samples_per_chirp = length(tx_wave);            % 2^15 mew, 33002 for 30MSps(old)

PRI = samples_per_chirp / chirp_sr;
PRF = 1 / PRI;
dt = 1/chirp_sr;
dR = (physconst('LightSpeed') * dt)/2;
R_margin = 500;                                 % meter of margin in RC                 
samp_margin = round(R_margin / dR);
%% LOAD chunk
addpath(genpath([pwd, filesep, 'lib' ]));       % add path of lib

exp_num = 3; cut_num = 1;
RC_file_name = strcat(folder_name,...
    experiment_name,num2str(exp_num),'_cut',num2str(cut_num));
RC = load_bin(RC_file_name);
% load last edited date
last_mod_date = load(strcat(folder_name,...
    experiment_name,num2str(exp_num),'_last_mod_date')).last_mod_date;

start_datetime = last_mod_date - seconds(PRI) * (size(RC,2)-1);
radar_slow_time = start_datetime + seconds(PRI) * (0:size(RC,2)-1);

save(strcat(folder_name,...
    experiment_name,num2str(exp_num),'_slow_time'),"radar_slow_time");
%% PLOT
R_ax = dR * (0:size(RC,1)-1);
slow_ax = PRI*(0:size(RC,2)-1);
figure,
imagesc(slow_ax,R_ax,abs(RC))