function [radar] = freqShiftCorrection(const,radar)
%FREQSHIFTCORRECTION correct the frequency shift by multi step approach
%   

wind_size = 2^8;
windowed_phase = windowedPhase(const,radar,wind_size);

% figure,plot(tau_ax,unwrap(angle(cross_talk))),hold on, plot(tau_ax,windowed_phase),  
% legend('Unwrap angle peak','Computed phase'), 
% title('Check if phase value is correct'), xlabel('Slow time [samples]')
% xlabel('tau axis [s]'), ylabel('Phase axis [Hz]')

%% First pass phase correction 
phase_shift = windowed_phase(:);
phasor = exp(-1i*phase_shift);
peak_wind_fixed = radar.cross_talk.value(:).* phasor(:);

% figure
% plot(tau_ax,angle(cross_talk)),title("Crosstalk phase")
% xlabel('\tau [s]'), ylabel('\phi [rad]')
% saveas(gcf,'Figures/phase_orig.png')

% figure
% plot(tau_ax,angle(peak_wind_fixed)),title("Compensated crosstalk phase (FFT window approach)")
% xlabel('\tau [s]'), ylabel('\phi [rad]')
% saveas(gcf,'Figures/phase_comp.png')

%% FULL matrix correction
phasor_mat = repmat(phasor(:).',size(radar.RC,1),1);
radar.RC = radar.RC .* phasor_mat;

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
figure,
subplot(3,1,1)
plot(radar.tau_ax,angle(radar.cross_talk.value)),title("Crosstalk phase")
 xlabel('\tau [s]'), ylabel('\phi [rad]')
subplot(3,1,2)
plot(radar.tau_ax,angle(peak_wind_fixed)),title("Compensated crosstalk phase (FFT window approach)")
xlabel('\tau [s]'), ylabel('\phi [rad]')
subplot(3,1,3)
plot(radar.tau_ax,angle(peak_fixed)),title("Averaged crosstalk phase")
xlabel('\tau [s]'), ylabel('\phi [rad]')
% saveas(gcf,'Figures/phase_corr_passage.png')

%% FULL matrix correction
phasor_mat = repmat(phasor(:).',size(radar.RC,1),1);
radar.RC = radar.RC .* phasor_mat;
clear phasor phasor_mat

%% FINAL phase correction with distance computed
phase_shift = -radar.cross_talk.phase_corr;
phasor = exp(-1i*phase_shift);
peak_final = peak_fixed(:).* phasor(:);

% figure
% plot(radar.tau_ax,angle(peak_final)),title("Fixed crosstalk phase")
% xlabel('\tau [s]'), ylabel('\phi [rad]')
% saveas(gcf,'Figures/phase_fixed.png')

% figure,
% plot(tau_ax,unwrap(angle(peak_final))),hold on, plot(radar.tau_ax,radar.cross_talk.phase_corr)
% legend("Crosstalk phase from distance value","Fixed crosstalk phase")
% xlabel('\tau [s]'), ylabel('\phi [rad]')
% saveas(gcf,'Figures/phase_fixed_check.png')

%% FULL matrix correction
phasor_mat = repmat(phasor(:).',size(radar.RC,1),1);
radar.RC = radar.RC .* phasor_mat;

end

