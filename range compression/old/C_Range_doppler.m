% clear all
close all
clc
%% COSTANTS
norm_B = .9;
chirp_sr = 56e6;                                % SDR sample rate
chirp_bw = norm_B*chirp_sr;                         % actual chirp bandwidth 

experiment_name = "test1_cut1";
% folder_name = 'mat_files/radar_window_test/';
folder_name = 'mat_files/giuriati_test/22_03_18/RC/cut/fixed/';

tx_wave = load(strcat('tx_waveform/tx_waveform_S56M.mat')).s_pad;tx_wave = single(tx_wave);
samples_per_chirp = length(tx_wave);            % 22002 for 20 MSps, 33002 for 30MSps
f0 = 1.65e9;
lambda = physconst('LightSpeed')/f0;
PRI = samples_per_chirp / chirp_sr;
PRF = 1 / PRI;
dt = 1/chirp_sr;
dR = (physconst('LightSpeed') * dt)/2;
R_margin = 500;                     % .5km of margin for RC
samp_margin = floor(R_margin/dR);
OSF = 8;
R_zero_shift = 10;                  % first peak real physical position
%% RC Load
addpath(genpath([pwd, filesep, 'lib' ]));       % add path of lib

RC_file = strcat(folder_name,experiment_name);

RC = load_bin(RC_file);

%% PLOTTING
R_ax = -R_margin:dR/OSF:R_margin - dR/OSF;
samp_shift = floor(R_zero_shift/dR * OSF);
tau_ax = linspace(0,PRI * size(RC,2),size(RC,2));

%% Windowed dopppler
dop_win_size = 2^nextpow2(round(1 / PRI));     % window time of 1s at least
win_idx = 1:dop_win_size;
N_cycle = floor(size(RC,2)/dop_win_size)*2-1;

slow_freq = linspace(-PRF/2,PRF/2,dop_win_size);
speed_ax = slow_freq * lambda/2 * 3.6;
x_idx = and(speed_ax<20,speed_ax>-20);
y_idx = and(R_ax>=-10,R_ax<150);

speed_ax_cut = speed_ax(x_idx);
R_ax_cut = R_ax(y_idx);


fig = figure;
images ={}; snapshots = {};
RC_shifted = circshift(RC,samp_shift,1);
% define idx where taget run around 10km/h
run1 = 13:60;
run2 = 82:100;
for i = run1%1:N_cycle
    win_idx = (1:dop_win_size) + (i-1) *(dop_win_size/2);
    
    RC_doppler = fftshift(fft(RC_shifted(:,win_idx),[],2),2);
    AA = abs(RC_doppler);
    %AA = 20*log10(AA);
    snapshots{end + 1} = RC_doppler(y_idx,x_idx);
    pre_dop = repmat(sum(abs(RC_shifted(:,win_idx)),2),1,size(AA,2));
    AA = AA ./ pre_dop;
    imagesc(speed_ax_cut,R_ax_cut,AA(y_idx,x_idx),[.1 1]);colormap jet, colorbar
    title(['Range-speed' num2str(i)]),xlabel("speed [km/h]"),ylabel("Range [m]")
    drawnow
    frame = getframe(fig);
    images{end + 1} = frame2im(frame);
    pause(.5)
end
%% FFT of moving target
N_fft = 4*2^nextpow2(size(snapshots{1},1));
figure
target_idx = [];
peak_idx_wind = 24:44;          % empirical idx where target is
for i = 1:length(snapshots)
    RC_fft = fftshift(fft(snapshots{i},N_fft,1),1);
    freq_ax = (-N_fft/2:N_fft/2-1)/N_fft * chirp_sr * OSF;
    temp = max(abs(RC_fft),[],1);
    temp(and(speed_ax_cut<5,speed_ax_cut>-5)) = 0;
    [~, target_idx(end + 1)] = max(temp(peak_idx_wind)); target_idx(end) = target_idx(end) +peak_idx_wind(1); 
     
    
    subplot(3,1,1)
    imagesc(speed_ax_cut,R_ax_cut,20*log10(abs(snapshots{i}))), colorbar, colormap jet
    title("Range - doppler")
    subplot(3,1,2)
    imagesc([],freq_ax,20*log10(abs(RC_fft))),colorbar,colormap jet, 
    title("Freq - doppler")
    subplot(3,1,3)
    
%     ff = abs(RC_fft(:,target_idx(end)));
    ff = mean(abs(RC_fft(:,peak_idx_wind)) .^2,2);
    plot(freq_ax,ff / max(ff)), hold on
    title(['FFT of Target in ' num2str(speed_ax_cut(target_idx(end))) ' km/h'])
    xlim([-chirp_sr chirp_sr]),grid on
    xlabel("Freq [Hz]")
    drawnow
    pause(.5)
end
%% PLOT peak descent
tgt_ampl_range = {};
figure
for i = 1:length(target_idx)
    tgt_ampl = abs(snapshots{i}(:,target_idx(i)));
    subplot(2,1,1)
    hold on 
    p1 = plot(R_ax_cut,tgt_ampl);
    [v,idx] = max(tgt_ampl);
    tgt_ampl_range{end +1,1} = R_ax_cut(idx); tgt_ampl_range{end,2} = v;
    plot(R_ax_cut(idx),v,'o')
    title("Target amplitude"), grid on
%     ylim([0 1e8]), xlabel("range")
    subplot(2,1,2)
    semilogy(R_ax_cut,20*log10(tgt_ampl)), grid on
    title("Target amplitude dB")
%     ylim([100 170]), xlabel("range")
    
    pause(.5)
    set(p1,'Visible','off')
end

%% Amplitude descent
tgt_range = 20 * log10(cell2mat(tgt_ampl_range(:,1)));
tgt_ampl = cell2mat(tgt_ampl_range(:,2));
[ampl_fitted,coeff] = poly_regression(tgt_range(:),20*log10(tgt_ampl(:)),1);

figure,

plot(tgt_range,20*log10(tgt_ampl) ),grid on
hold on
plot(tgt_range,ampl_fitted)
title(['Target Amplitude decay, exp: ' num2str(coeff(2)) ]),
xlabel("20log10(Range) [dB]"),ylabel ("Amplitude [dB]")
legend("Target amplitude","Fitted line")


return
%% CREATE GIF
filename = strcat(experiment_name,'.gif'); % Specify the output file name
for idx = 1:N_cycle
    [A,map] = rgb2ind(images{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.1);
    end
end