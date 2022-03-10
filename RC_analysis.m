% clear all
close all
clc
%% COSTANTS
chirp_bw = 27e6;                                % actual chirp bandwidth 
chirp_sr = 30e6;                                % SDR sample rate
norm_B = chirp_bw / chirp_sr;
% experiment_name = strcat('MULTISTATIC_GPS1_B',...
%     num2str(chirp_bw/1e6),'M_S',num2str(chirp_sr/1e6),'M');
experiment_name = "data5_cut1";
% folder_name = 'mat_files/radar_window_test/';
folder_name = 'mat_files/giuriati_test_22_02_25/RC/cut/';

tx_wave = load(strcat('tx_waveform/tx_waveform_2pow_B',...
    num2str(chirp_bw/1e6),'M_S',num2str(chirp_sr/1e6),'M.mat')).s_pad;
tx_wave = single(tx_wave);
samples_per_chirp = length(tx_wave);            % 22002 for 20 MSps, 33002 for 30MSps
f0 = 1.65e9;
lambda = physconst('LightSpeed')/f0;
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
peak_idxs = movmean(peak_idxs,1000); 
peak_shifts = peak_idxs -peak_idxs(1);
%% Remove noise from peak
deriv = [diff([diff(peak_shifts) 0]) 0];
deriv = deriv/max(deriv);
peak_shifts_filt = zeros(size(deriv));
safe_counter = -1;
for i = 1:length(deriv)
    if abs(deriv(i)) < 0.3 
        if safe_counter <0
            good_val = peak_shifts(i);
        end
    else
        safe_counter = 5e3;                   % Noisy zone, don't use this value
    end
    peak_shifts_filt(i) = good_val;
    safe_counter =safe_counter - 1;
end
peak_shifts_filt = movmean(peak_shifts_filt,5e3);

clear safe_counter good_val deriv
%% ESTIMATE err1
% estimate only checking starting and ending idx of the peak
diff = abs(peak_shifts_filt(end));
err1= diff / (samples_per_chirp * OSF *  size(RC,2));
clear diff

%% TIME SHIFT CORRECTION
t = 0:size(RC,1) -1;
tf = 0:1/OSF:(size(RC,1) -1);

tf_men_t =tf(:) - t;
W = sinc(norm_B*tf_men_t);
% interpolate whole RC signal
RC_OS = W*RC;

Nf = 2^nextpow2(size(RC_OS,1));
X = fftshift(fft(RC_OS,Nf,1),1);

f = (-Nf/2:Nf/2 - 1)/Nf;
H = exp(1i*2*pi*f(:)*peak_shifts_filt);
H(1,:) = 0;                         % 

RC_Dt_fixed = ifft(ifftshift(X.* H,1),Nf,1);
RC_Dt_fixed = RC_Dt_fixed(1:size(RC_OS,1),:);

%% VARIABLE Df computation
win_size = 1024;

[ ~, idx ] = max(RC_Dt_fixed(:,1));
peak = RC_Dt_fixed(idx,:);
% compare magnitude of peak with closer values in t to see its stability
% figure,plot(abs(peak)),hold on, plot(abs(RC_Dt_fixed(idx + OSF,:))),plot(abs(RC_Dt_fixed(idx - OSF,:)))
% legend("Peak at start","Peak + 1 samples","Peak - 1 samples")
%% Cycle 
% number of shifts of the window over the peak
N_cycle = floor(length(peak)/win_size) * 2 - 1;
Df = zeros(N_cycle,1); fft_max = Df; fft_angle = Df;
N_fft = 8*win_size;
freq_ax = (-N_fft/2:N_fft/2-1)/N_fft * PRF;
win_idx = 1 : win_size;

for i = 1:N_cycle
    aa = peak(win_idx);
    ff = fftshift(fft(aa,N_fft));
    [fft_max(i),ii] = max(abs(ff));
    fft_angle(i) = angle(ff(ii));
    Df(i) = freq_ax(ii);
    ang = unwrap(angle(aa));
    Df_(i) = (ang(end)-ang(1)) / win_size/PRI/2/pi;         %% Compute Df as phase difference over the window
    win_idx = win_idx + win_size/2;
% figure,plot(freq_ax,abs(ff));
% pause
end
% figure,plot(Df_),hold on,plot(Df)

%% INTEGRATE Df to get phase
computed_phase=[];pp = angle(peak(1));
t = 0:PRI:(win_size/2-1)*PRI;
for i = 1:length(Df)
    pp = pp(end) + 2*pi*Df(i)*t(:);
    computed_phase = [computed_phase; pp];
end
if length(computed_phase)<length(peak)
    temp = repmat(computed_phase(end),length(peak)-length(computed_phase),1);
    computed_phase = [computed_phase(:); temp]; 
end
% figure,plot(unwrap(angle(peak))),hold on, plot(computed_phase)

%% FFT PLOT for peak,angle,Df
% figure,plot(Df);title("Df values"),xlabel("Sliding window iteration"),ylabel("Frequency [Hz]"),grid on
% figure,plot(fft_max);title("FFT peak amplitude");xlabel("Sliding Window iteration"), ylabel("Amplitude"),grid on
% figure,plot(fft_angle);title("FFT peak angle ");xlabel("Sliding Window iteration"),ylabel("Phase [rad]"),grid on

%% Df FIX - FFT windowed based
phasor = exp(-1i*computed_phase);
peak_fixed = peak(:).* phasor(:);
%% FULL Df fix
phasor_mat = repmat(phasor(:).',size(RC_Dt_fixed,1),1);
RC_Df_fixed = RC_Dt_fixed .* phasor_mat;

%% Final Correction
peak_phase = angle(peak_fixed);
peak_phase_avg = movmean(peak_phase,2e3);

phasor = exp(-1i*peak_phase_avg);
peak_fixed = peak_fixed(:).* phasor(:);

phasor_mat = repmat(phasor(:).',size(RC_Df_fixed,1),1);
RC_Df_fixed = RC_Df_fixed .* phasor_mat;

%% PLOTTING
R_ax = -R_margin:dR/OSF:R_margin - dR/OSF;
samp_shift = floor(R_zero_shift/dR * OSF);
tau_ax = linspace(0,PRI * size(RC,2),size(RC,2));
%% Amplitude
figure,subplot(2,1,1)
imagesc(tau_ax,R_ax,abs(circshift(RC_OS,samp_shift,1)));
title('RC original' ),xlabel("Slow time [s]"),ylabel("Range [m]")
subplot(2,1,2)
imagesc(tau_ax,R_ax,abs(circshift(RC_Df_fixed,samp_shift,1)));
title('RC fixed' ),xlabel("Slow time [s]"),ylabel("Range [m]")
%% Phase 
figure,subplot(2,1,1)
imagesc(tau_ax,R_ax,angle(circshift(RC_Dt_fixed,samp_shift,1)));
title('RC Dt fixed' ),xlabel("Slow time [s]"),ylabel("Range [m]")
subplot(2,1,2)
imagesc(tau_ax,R_ax,angle(circshift(RC_Df_fixed,samp_shift,1)));
title('RC freq fixed' ),xlabel("Slow time [s]"),ylabel("Range [m]")
%% Peak Phase
figure,plot(angle(peak)),hold on,plot(angle(peak_fixed)); legend("Original","fixed")

%% Coherent SUM
coh_sum = sum(RC_Df_fixed,2);
[ ~, zero_idx ] = max(coh_sum);

coh_sum = circshift(coh_sum,samp_shift - (zero_idx-round(length(coh_sum)/2)));
figure,plot(R_ax,abs(coh_sum));grid on; 
title("Coherent sum of RC centered"),xlabel("Range [m]"),ylabel("Amplitude")

%% Windowed dopppler
dop_win_size = 2^nextpow2(round(1 / PRI));     % window time of 1s at least
win_idx = 1:dop_win_size;
N_cycle = floor(size(RC_Df_fixed,2)/dop_win_size)*2-1;

slow_freq = linspace(-PRF/2,PRF/2,dop_win_size);
speed_ax = slow_freq * lambda/2 * 3.6;
x_idx = and(speed_ax<20,speed_ax>-20);
y_idx = and(R_ax>=-10,R_ax<100);

x = speed_ax(x_idx);
y = R_ax(y_idx);


fig = figure;
images ={}; 
RC_shifted = circshift(RC_Df_fixed,samp_shift,1);
for i = 1:N_cycle
    RC_doppler = fftshift(fft(RC_shifted(:,win_idx),[],2),2);
    RC_doppler = 10*log10(RC_doppler);
    AA = abs(RC_doppler(y_idx,x_idx));

    
    win_idx = win_idx + dop_win_size/2;
    imagesc(x,y,AA);
    title("Range-speed"),xlabel("speed [km/h]"),ylabel("Range [m]")
    drawnow
    frame = getframe(fig);
    images{i} = frame2im(frame);
end
close
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