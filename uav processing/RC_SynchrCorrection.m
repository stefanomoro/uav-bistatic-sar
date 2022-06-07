%% ======================================================================== INTERPOLATION
OSF = 8; 

t = 0:size(RC,1) -1;
tf = 0:1/OSF:(size(RC,1) -1);

tf_men_t =tf(:) - t;
W = sinc(norm_B*tf_men_t);

% interpolate whole RC signal
RC_OS = W*RC;
R_ax_OS = R_ax(1):dR/OSF:R_ax(end);

% figure
% imagesc(tau_ax,R_ax_OS,abs(RC_OS))
% title('RC interp - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

clear t tf tf_men_t W

%% Plot
% figure
% subplot(2,1,1)
% imagesc(tau_ax,R_ax,abs(RC_OS))
% title('Range compressed matrix - absolute value' )
% xlabel("\tau [s]"),ylabel("Range [m]")
% subplot(2,1,2)
% imagesc(tau_ax,R_ax,angle(RC_OS))
% title('Range compressed matrix - phase' )
% xlabel("\tau [s]"),ylabel("Range [m]")
% saveas(gcf,'Figures/RC_OS.png')

%%  =======================================================================  TIME SHIFT CORRECTION
%% compute correct cross-talk idx from distance
% define the idxs where the cross talk should be

cross_talk_idxs_corr = zeros(size(tau_ax));
for n = 1:length(tau_ax)
    idxs = find(R_ax_OS > distance_tx_rx(n));
    cross_talk_idxs_corr(n) = idxs(1);
end

cross_talk_idxs_corr = movmean(cross_talk_idxs_corr,1e3);
clear idxs

%% Cross-talk tracking
[~,cross_talk_idxs] = max(RC_OS,[],1);
% apply median filter and then moving average
cross_talk_idxs = movmean(medfilt1(cross_talk_idxs,1e3) ,1000); 

%% TIME SHIFT CORRECTION
cross_talk_shifts = cross_talk_idxs - cross_talk_idxs_corr;

Nf = 2^nextpow2(size(RC_OS,1));
X = fftshift(fft(RC_OS,Nf,1),1);

f = (-Nf/2:Nf/2 - 1)/Nf;
H = exp(1i*2*pi*f(:)*cross_talk_shifts);
H(1,:) = 0;                         

RC_Dt_fixed = ifft(ifftshift(X.* H,1),Nf,1);
RC_Dt_fixed = RC_Dt_fixed(1:size(RC_OS,1),:);

% figure
% imagesc(tau_ax,R_ax,abs(RC_Dt_fixed))
% title('RC Dt fixed - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

clear Nf X f H
%%  ======================================================================= PHASE CORRECTION
%% Find Crosstalk
% get crosstalk
[~,cross_talk_idxs] = max(RC_Dt_fixed,[],1);
% apply median filter and then moving average
cross_talk_idxs = round(movmean(medfilt1(cross_talk_idxs,1e3) ,1e3)); 

cross_talk = zeros(size(RC_Dt_fixed,2),1);
for i = 1:length(cross_talk)
    cross_talk(i) = RC_Dt_fixed(cross_talk_idxs(i),i);
end

%% Compute correct cross talk phase from distance
% Correct phase
c = physconst('Lightspeed');
phase_corr = -2*pi*f0*(distance_tx_rx/c); % - 2*pi*f0*(2*speed_distance_tx_rx/lambda);
phase_corr = phase_corr(:).';

% move first phase in 2pi interval
phase_0 = mod(phase_corr(1),2*pi)-pi;

phase_corr = phase_corr - phase_corr(1) + phase_0;

% figure, plot(tau_ax, unwrap(angle(cross_talk))),
% hold on, plot(tau_ax, phase_corr),  xlabel('\tau [s]'), ylabel('\phi [rad]')
% legend('Crosstalk phase','Correct phase')
% saveas(gcf,'Figures/phase_cross_corr.png')

clear phase_0
%%  ===================== PHASE CORRECTION

wind_size = 2^8;
windowed_phase = windowedPhase(cross_talk,PRI,wind_size);

% figure,plot(tau_ax,unwrap(angle(cross_talk))),hold on, plot(tau_ax,windowed_phase),  
% legend('Unwrap angle peak','Computed phase'), 
% title('Check if phase value is correct'), xlabel('Slow time [samples]')
% xlabel('tau axis [s]'), ylabel('Phase axis [Hz]')

%% First pass phase correction 
phase_shift = windowed_phase(:);
phasor = exp(-1i*phase_shift);
peak_wind_fixed = cross_talk(:).* phasor(:);

% figure
% plot(tau_ax,angle(cross_talk)),title("Crosstalk phase")
% xlabel('\tau [s]'), ylabel('\phi [rad]')
% saveas(gcf,'Figures/phase_orig.png')

% figure
% plot(tau_ax,angle(peak_wind_fixed)),title("Compensated crosstalk phase (FFT window approach)")
% xlabel('\tau [s]'), ylabel('\phi [rad]')
% saveas(gcf,'Figures/phase_comp.png')

%% FULL matrix correction
phasor_mat = repmat(phasor(:).',size(RC_Dt_fixed,1),1);
RC_Df_fixed = RC_Dt_fixed .* phasor_mat;

clear phasor phasor_mat

%% Second step phase correction
% average the phase and leave only noise
peak_mean = movmean(peak_wind_fixed,500);
phase_shift = angle(peak_mean); 
phasor = exp(-1i*phase_shift);
peak_fixed = peak_wind_fixed(:).* phasor(:);

% figure
% plot(tau_ax,angle(peak_fixed)),title("Averaged crosstalk phase")
% xlabel('\tau [s]'), ylabel('\phi [rad]')
% saveas(gcf,'Figures/phase_avg.png')

%% Plot
% figure,
% subplot(3,1,1)
% plot(tau_ax,angle(cross_talk)),title("Crosstalk phase")
%  xlabel('\tau [s]'), ylabel('\phi [rad]')
% subplot(3,1,2)
% plot(tau_ax,angle(peak_wind_fixed)),title("Compensated crosstalk phase (FFT window approach)")
% xlabel('\tau [s]'), ylabel('\phi [rad]')
% subplot(3,1,3)
% plot(tau_ax,angle(peak_fixed)),title("Averaged crosstalk phase")
% xlabel('\tau [s]'), ylabel('\phi [rad]')
% saveas(gcf,'Figures/phase_corr_passage.png')

%% FULL matrix correction
phasor_mat = repmat(phasor(:).',size(RC_Dt_fixed,1),1);
RC_Df_fixed = RC_Df_fixed .* phasor_mat;
clear phasor phasor_mat

%% FINAL phase correction with distance computed
phase_shift = -phase_corr;
phasor = exp(-1i*phase_shift);
peak_final = peak_fixed(:).* phasor(:);

% figure
% plot(tau_ax,angle(peak_final)),title("Fixed crosstalk phase")
% xlabel('\tau [s]'), ylabel('\phi [rad]')
% saveas(gcf,'Figures/phase_fixed.png')
% 
% figure,
% plot(tau_ax,unwrap(angle(peak_final))),hold on, plot(tau_ax,phase_corr)
% legend("Crosstalk phase from distance value","Fixed crosstalk phase")
% xlabel('\tau [s]'), ylabel('\phi [rad]')
% saveas(gcf,'Figures/phase_fixed_check.png')

%% FULL matrix correction
phasor_mat = repmat(phasor(:).',size(RC_Dt_fixed,1),1);
RC_Df_fixed = RC_Df_fixed .* phasor_mat;
clear phasor phasor_mat

%%  ======================================================================= FINAL CORRECTION RESULTS
%% Amplitude
% figure,subplot(2,1,1)
% imagesc(tau_ax,R_ax_OS,abs(RC_OS));
% title('Range compressed matrix - absolute value'),xlabel("\tau [s]"),ylabel("Range [m]")
% subplot(2,1,2)
% imagesc(tau_ax,R_ax_OS,abs(RC_Df_fixed));
% title('Range compressed matrix synchronism correction - absolute value'),xlabel("\tau [s]"),ylabel("Range [m]")
% saveas(gcf,'Figures/RC_fixed_abs.png')

%% Phase 
% figure
% subplot(2,1,1)
% imagesc(tau_ax,R_ax_OS,angle(RC_OS));
% title('Range compressed matrix - phase' ),xlabel("\tau [s]"),ylabel("Range [m]")
% subplot(2,1,2)
% imagesc(tau_ax,R_ax_OS,angle(RC_Df_fixed));
% title('Range compressed matrix synchronism correction - phase' ),xlabel("\tau  [s]"),ylabel("Range [m]")
% saveas(gcf,'Figures/RC_fixed_angle.png')

%% ======================================================================== OUTPUT
RC = RC_Df_fixed;
R_ax = R_ax_OS;

%% Plot
% figure
% subplot(2,1,1)
% imagesc(tau_ax,R_ax,abs(RC))
% title('Range compressed matrix synchronism correction - absolute value' )
% xlabel("\tau [s]"),ylabel("Range [m]")
% subplot(2,1,2)
% imagesc(tau_ax,R_ax,angle(RC))
% title('Range compressed matrix synchronism correction - phase' )
% xlabel("\tau [s]"),ylabel("Range [m]")
% saveas(gcf,'Figures/RC_fixed.png')



