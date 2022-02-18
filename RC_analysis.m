clear all
close all
clc
%% COSTANTS
chirp_bw = 27e6;                                % actual chirp bandwidth 
chirp_sr = 30e6;                                % SDR sample rate
experiment_name = strcat('MULTISTATIC_GPS1_B',...
    num2str(chirp_bw/1e6),'M_S',num2str(chirp_sr/1e6),'M');

mat_file_raw = strcat('mat_files/radar_window_test/',...
    experiment_name,'/chunk');

mat_file = strcat('mat_files/radar_window_test/',...
    experiment_name,'/RC/chunk');

tx_wave = load(strcat('tx_waveform/tx_waveform_B',...
    num2str(chirp_bw/1e6),'M_S',num2str(chirp_sr/1e6),'M.mat')).s_pad;
tx_wave = single(tx_wave);
samples_per_chirp = length(tx_wave);            % 22002 for 20 MSps, 33002 for 30MSps
f0 = 2.03e9;
PRI = samples_per_chirp / chirp_sr;
PRF = 1 / PRI;
dt = 1/chirp_sr;
dR = (physconst('LightSpeed') * dt)/2;
R_margin = 1e3;                     % 1km of margin for RC
samp_margin = floor(R_margin/dR);

N_input_lags = 2^nextpow2(samples_per_chirp);
R_zero_shift = 57;                  % first peak is at 57m
%% RC good up to : 98273
addpath(genpath([pwd, filesep, 'lib' ]));       % add path of lib

result_file = strcat('mat_files/radar_window_test/',...
    experiment_name,'/RC/RC_total');

RC = load_bin(result_file);

% WARNING : empirical value to cut RC signal
disp("!!! WARNING !!! Change empirical value")
RC = RC(:,1:98273);
%% ESTIMATE err1
% estimate only checking starting and ending idx of the peak

[~,st_pos] = max(abs(RC(:,1)));
[~,end_pos] = max(abs(RC(:,end)));
diff = abs(end_pos - st_pos);
err1= diff / (samples_per_chirp * size(RC,2));
clear diff st_pos end_pos
%% PEAK time-shift tracking
[~,peak_idx] = max(abs(RC(:,1)));
RC_peak = RC(peak_idx-10:peak_idx+10,:);
OSF = 8;                                        %over sampling factor
%% Peak interpolation
x = 0:dt:dt*(size(RC_peak,1) -1);
xq = 0:dt/OSF:dt/OSF*(size(RC_peak,1) * OSF -1);

tf_men_t =xq(:) - x;
W = sinc(chirp_bw*tf_men_t);
RC_peak_OS = W*RC_peak;
%%
[~,peak_idxs] = max(RC_peak_OS);
peak_idxs = movmean(peak_idxs,2000); 
peak_shifts = peak_idxs -peak_idxs(1);
% ord = 1;
% poly_fit=   polyfit(0:length(peak_idxs)-1,peak_idxs,ord);
% fit_y = polyval(poly_fit,0:length(peak_idxs)-1);
% figure,plot(0:length(peak_idxs)-1,fit_y,0:length(peak_idxs)-1,peak_idxs)
% figure,imagesc([],xq,abs(RC_peak_OS))

%% TIME SHIFT CORRECTION
norm_B = chirp_bw / chirp_sr;
t =single( 0:size(RC_peak_OS,1)-1);
tf = t;
tf_men_t =tf(:) - t;
tf_men_Dt = zeros(size(tf_men_t,1),size(tf_men_t,2),size(RC_peak_OS,2),'single');
for i = 1:size(RC_peak_OS,2)
    tf_men_Dt(:,:,i)= tf_men_t + peak_shifts(i);
end

W = sinc(norm_B*tf_men_Dt);

RC_peak_fix = zeros(size(RC_peak_OS),'like',RC_peak_OS);
for i = 1:size(RC_peak_OS,2)
    RC_peak_fix(:,i) = W(:,:,i)*RC_peak_OS(:,i);
end

figure,subplot(2,1,1);
imagesc(abs(RC_peak_OS))

subplot(2,1,2);
imagesc(abs(RC_peak_fix))

%% PLOTTING
R_ax = -R_margin:dR:R_margin;
samp_shift = floor(R_zero_shift/dR);
tau_ax = linspace(0,PRI * size(RC,2),size(RC,2));
plot_idx = find(and(R_ax >= 0,R_ax <= 200));

% figure,imagesc(tau_ax,R_ax,abs(RC_centered))
% title('RC time corrected' ),xlabel("Slow time [s]"),ylabel("Range [m]")

figure,subplot(2,1,1);
imagesc(tau_ax,R_ax,abs(circshift(RC,samp_shift,1)))
title('RC original' ),xlabel("Slow time [s]"),ylabel("Range [m]")

subplot(2,1,2);
imagesc(tau_ax,R_ax,abs(circshift(RC_fix,samp_shift,1)));
title('RC time fixed' ),xlabel("Slow time [s]"),ylabel("Range [m]")

%% Check window 
% check with which window size we obtain the best phase linearity
% [ ~, idx ] = max(RC_centered(:,1));
% 
% peak = RC_centered(idx,:);
% win_size = 1e1:.5e3:1e4;
% mm = zeros(length(win_size),1);
% for j = 1:length(win_size)
%     pipp = [];
% for i = 1:floor(length(peak)/win_size(j))
%     aa = peak(win_size(j) * (i-1) + 1:win_size(j) * i );
%     pipp(i) = max(abs(fft(aa))/ sum(abs(aa)));
% end
% mm(j) = mean(pipp);
% end
% figure,plot(win_size,mm)
%% VARIABLE Df computation
win_size = 1024;
RC_centered = circshift(RC_fix,samp_shift,1);

[ ~, idx ] = max(RC_centered(:,1));
peak = RC_centered(idx,:);
% compare magnitude of peak with closer values in t to see its stability
figure,plot(abs(peak)),hold on, plot(abs(RC_centered(idx + 1,:))),plot(abs(RC_centered(idx - 1,:)))
legend("Peak at start","Peak + 1 samples","Peak - 1 samples")
%%
% number of shifts of the window over the peak
N_cycle = floor(length(peak)/win_size) * 2 - 1;
Df = zeros(N_cycle,1);
fft_max = Df;   fft_angle = Df;
N_fft = 8*win_size;
freq_ax = linspace(-.5,.5,N_fft) * PRF;
win_idx = 1 : win_size;
for i = 1:N_cycle
    aa = peak(win_idx);
    ff = fftshift(fft(aa,N_fft));
    [fft_max(i),ii] = max(abs(ff));
    fft_angle(i) = angle(ff(ii));
    Df(i) = freq_ax(ii);
    win_idx = win_idx + win_size/2;
% 	figure,plot(freq_ax,abs(ff));
end

%%
figure,plot(Df);title("Df values"),xlabel("Sliding window iteration"),ylabel("Frequency [Hz]"),grid on
figure,plot(fft_max);title("FFT peak amplitude");xlabel("Sliding Window iteration"), ylabel("Amplitude"),grid on
figure,plot(fft_angle);title("FFT peak angle ");xlabel("Sliding Window iteration"),ylabel("Phase [rad]"),grid on

%% FIT Df curve
Df_x = linspace(1,length(peak),length(Df));
Df_x1 = 1:length(peak);
ord = 20;
Df_poly=   polyfit(Df_x,Df,ord);
Df_y1 = polyval(Df_poly,Df_x1);
figure,plot(Df_x,Df,Df_x1,Df_y1);
title(strcat("Order ",num2str(ord)," poly fitting "))

%% Df fix
% phasor = exp(-1i* 2*pi * Df_y1 .* tau_ax);
phasor = exp(-1i*angle(peak));
% peak_fixed = peak.* phasor;
% figure,plot(unwrap(angle(peak))),hold on,plot(unwrap(angle(peak_fixed)));

phasor_mat = repmat(phasor,size(RC_centered,1),1);
RC_Df_fixed = RC_centered .* phasor_mat;
%% Plot
figure,
subplot(2,1,1);imagesc(tau_ax,R_ax,angle(RC_centered));title("Original")
subplot(2,1,2);imagesc(tau_ax,R_ax,angle(RC_Df_fixed));title("Df fixed")
figure,
subplot(2,1,1);imagesc(tau_ax,R_ax,abs(RC_centered));title("Original")
subplot(2,1,2);imagesc(tau_ax,R_ax,abs(RC_Df_fixed));title("Df fixed")
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

