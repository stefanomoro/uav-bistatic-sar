clear all
close all
clc
%% COSTANTS
chirp_bw = 27e6;                                % actual chirp bandwidth 
chirp_sr = 30e6;                                % SDR sample rate
norm_B = chirp_bw / chirp_sr;
% experiment_name = strcat('MULTISTATIC_GPS1_B',...
%     num2str(chirp_bw/1e6),'M_S',num2str(chirp_sr/1e6),'M');
experiment_name = "data6_cut1";
% folder_name = 'mat_files/radar_window_test/';
folder_name = 'mat_files/giuriati_test_22_02_25/RC/cut/';

tx_wave = load(strcat('tx_waveform/tx_waveform_2pow_B',...
    num2str(chirp_bw/1e6),'M_S',num2str(chirp_sr/1e6),'M.mat')).s_pad;
tx_wave = single(tx_wave);
samples_per_chirp = length(tx_wave);            % 22002 for 20 MSps, 33002 for 30MSps
f0 = 1.65e9;
PRI = samples_per_chirp / chirp_sr;
PRF = 1 / PRI;
dt = 1/chirp_sr;
dR = (physconst('LightSpeed') * dt)/2;
R_margin = 500;                     % .5km of margin for RC
samp_margin = floor(R_margin/dR);

R_zero_shift = 10;                  % first peak real physical position
%% RC Load
addpath(genpath([pwd, filesep, 'lib' ]));       % add path of lib

RC_file = strcat(folder_name,experiment_name);

RC = load_bin(RC_file);
% RC = RC(:,1:98273);
%% PEAK time-shift tracking
[~,peak_idx] = max(abs(RC(:,1)));
peak_margin = 5;
RC_peak = RC(peak_idx-peak_margin:peak_idx+peak_margin,:);
OSF = 8;                                        %over sampling factor
%% Peak interpolation
x = 0:dt:dt*(size(RC_peak,1) -1);
xq = 0:dt/OSF:dt/OSF*(size(RC_peak,1) * OSF -1);

tf_men_t =xq(:) - x;
W = sinc(chirp_bw*tf_men_t);
RC_peak_OS = W*RC_peak;
%% Peak max tracking
[~,peak_idxs] = max(RC_peak_OS);
peak_idxs = movmean(peak_idxs,2000); 
peak_shifts = peak_idxs -peak_idxs(1);

%% ESTIMATE err1
% estimate only checking starting and ending idx of the peak
diff = abs(peak_shifts(end));
err1= diff / (samples_per_chirp * OSF *  size(RC,2));
clear diff

%% TIME SHIFT CORRECTION
t = 0:size(RC,1) -1;
tf = 0:1/OSF:(size(RC,1) -1);

tf_men_t =tf(:) - t;
W = sinc(norm_B*tf_men_t);
% interpolate whole RC signal
RC_OS = W*RC;

t =0:size(RC_OS,1)-1;               % centered around zero t axis
t = t- t(round(size(RC_OS,1)/2));

Nf = 2^nextpow2(size(RC_OS,1));
X = fftshift(fft(RC_OS,Nf,1),1);

f = (-Nf/2:Nf/2 - 1)/Nf;
H = exp(1i*2*pi*f(:)*peak_shifts);
H(1,:) = 0;                         % 

RC_Dt_fixed = ifft(ifftshift(X.* H,1),Nf,1);
RC_Dt_fixed = RC_Dt_fixed(1:size(RC_OS,1),:);

%% Df correction
time_shifts = peak_shifts * dt/OSF;
polinom = polyfit(0:length(time_shifts)-1,time_shifts,7);
fitted_shifts = polyval(polinom,0:length(time_shifts)-1);
%%
phase_fixer = repmat(exp(-2i* pi*f0*fitted_shifts),size(RC_Dt_fixed,1),1);
RC_Df_fixed = RC_Dt_fixed .* phase_fixer;

%% PLOTTING
R_ax = -R_margin:dR/OSF:R_margin;
samp_shift = floor(R_zero_shift/dR * OSF);
tau_ax = linspace(0,PRI * size(RC,2),size(RC,2));

figure,subplot(2,1,1)
imagesc(tau_ax,R_ax,abs(circshift(RC_OS,samp_shift,1)));
title('RC original' ),xlabel("Slow time [s]"),ylabel("Range [m]")
subplot(2,1,2)
imagesc(tau_ax,R_ax,abs(circshift(RC_Dt_fixed,samp_shift,1)));
title('RC time fixed' ),xlabel("Slow time [s]"),ylabel("Range [m]")
%% phase 
figure,subplot(2,1,1)
imagesc(tau_ax,R_ax,angle(circshift(RC_Dt_fixed,samp_shift,1)));
title('RC Dt fixed' ),xlabel("Slow time [s]"),ylabel("Range [m]")
subplot(2,1,2)
imagesc(tau_ax,R_ax,angle(circshift(RC_Df_fixed,samp_shift,1)));
title('RC freq fixed' ),xlabel("Slow time [s]"),ylabel("Range [m]")
%% Peak Phase
figure,subplot(2,1,1),
plot(tau_ax,angle(RC_OS(peak_idx(1),:)))
title('Peak phase original' ),xlabel("Slow time [s]"),ylabel("Phase [rad]")

subplot(2,1,2),
plot(tau_ax,angle(RC_Df_fixed(peak_idx(1),:)))
title('Peak phase fixed' ),xlabel("Slow time [s]"),ylabel("Phase [rad]")

return
%% Coherent SUM
coh_sum = sum(RC_Df_fixed,2);

figure,plot(R_ax,abs(coh_sum));grid on; 
title("Coherent sum of RC centered"),xlabel("Range [m]"),ylabel("Amplitude")

%% Compute Doppler 
RC_doppler = fftshift(fft(RC_Df_fixed,[],2),2);
slow_freq = linspace(-PRF/2,PRF/2,size(RC_doppler,2));
figure,imagesc(slow_freq,R_ax,abs(RC_doppler));title("Range-doppler (freq)"),xlabel("Doppler freq [Hz]"),ylabel("Range [m]")
speed_ax = slow_freq/f0 * physconst('LightSpeed') * 3.6;
figure,imagesc(speed_ax,R_ax,abs(RC_doppler));title("Range-speed (m/s)"),xlabel("speed [km/h]"),ylabel("Range [m]")

