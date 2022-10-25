clear all
close all
clc
%% COSTANTS
chirp_sr = 40e6;                                % SDR sample rate
chirp_bw = .9*chirp_sr;                         % actual chirp bandwidth 
f0 = 1.65e9;
lambda = physconst('LightSpeed')/f0;
folder_name = '../mat_files/giuriati_test/20221020';

tx_wave = load(strcat('../tx_waveform/tx_waveform_S56M.mat')).s_pad;
tx_wave = single(tx_wave);
samples_per_chirp = length(tx_wave);            % 2^15 mew, 33002 for 30MSps(old)

PRI = samples_per_chirp / chirp_sr;
PRF = 1 / PRI;
dt = 1/chirp_sr;
dR = (physconst('LightSpeed') * dt);

ant_distance = 1.07;                            % distance between antennas
%% LOAD RAW
addpath("../lib");       % add path of lib
addpath("./utils")
allProcTimer = tic;
%% List Files
file_paths = listOnlyBinFiles(strcat(folder_name,'/RC/cut'));
%%

for exp_num =4 : length(file_paths)
file_path = file_paths(exp_num).complete_path;
disp(' '),disp(['Loading raw data ' num2str(exp_num)]),tic
RC =load_bin(file_path(1:end-3)); 
disp(['Loaded in ' num2str(toc) ' s']);


%%
RC_raw = RC;
OSF = 16;
RC = interpolateRC(RC_raw,OSF,.9);
t_ax = 0:dt/OSF:dt/OSF*size(RC,1);
[~,max_idx ] = max(RC(:,100));

zero_idx = floor(max_idx - ant_distance/(dR/OSF));
RC = RC(zero_idx:end,:);
R_ax = 0:dR/OSF/2:dR/OSF/2*(size(RC,1)-1);
% figure
% imagesc([],R_ax,abs(RC));


%% Windowed doppler (Normalized)
% dop_win_size = 2^nextpow2(round(1 / PRI));     % window time of 1s at least
% win_idx = 1:dop_win_size;
% N_cycle = floor(size(RC,2)/dop_win_size)*2-1;
% slow_freq = linspace(-PRF/2,PRF/2,dop_win_size);
% speed_ax = slow_freq * lambda/2 * 3.6;
% x_idx = and(speed_ax<10,speed_ax>-10);
% y_idx = and(R_ax>=-10,R_ax<100);
% 
% x = speed_ax(x_idx);
% y = R_ax(y_idx);
% 
% 
% fig = figure;
% images ={}; 
% for i = 1:N_cycle
%     RC_doppler = fftshift(fft(RC(:,win_idx),[],2),2);
%     AA = sum(abs(RC(:,win_idx)),2);
%     RC_doppler = RC_doppler./(AA*ones(1,length(win_idx)));
%     
% 
%     BB = abs(RC_doppler(y_idx,x_idx));
%    
%   
%     
%     win_idx = win_idx + dop_win_size/2;
%     imagesc(x,y,BB,[0 1]);
%    
%     title(num2str(i)),xlabel("speed [km/h]"),ylabel("Range [m]"), colormap jet, colorbar 
%     drawnow
%     frame = getframe(fig);
%     images{i} = frame2im(frame);
% end
% close
%%
%% Windowed doppler (Normalized dB)
dop_win_size = 2^nextpow2(round(1 / PRI));     % window time of 1s at least
win_idx = 1:dop_win_size;
N_cycle = floor(size(RC,2)/dop_win_size)*2-1;

slow_freq = linspace(-PRF/2,PRF/2,dop_win_size);
speed_ax = slow_freq * lambda/2 * 3.6;
x_idx = and(speed_ax<20,speed_ax>-20);
y_idx = and(R_ax>=-10,R_ax<200);

x = speed_ax(x_idx);
y = R_ax(y_idx);

fig = figure("WindowState","maximized");
images ={}; 
for i = 1:N_cycle
    RC_doppler = fftshift(fft(RC(:,win_idx),[],2),2);
    AA = 20*log10(abs(RC_doppler(y_idx,x_idx)));

    
    win_idx = win_idx + dop_win_size/2;
    imagesc(x,y,AA,[120 220]);
    title(num2str(i)),xlabel("speed [km/h]"),ylabel("Range [m]"), colormap jet, colorbar, 
    drawnow
    frame = getframe(fig);
    images{i} = frame2im(frame);
end
close

%% CREATE GIF
filename = strcat(file_paths(exp_num).exp_name,'.gif'); % Specify the output file name
for idx = 1:N_cycle
    [A,map] = rgb2ind(images{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',.1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.1);
    end
end
end