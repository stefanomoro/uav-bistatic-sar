function [outputArg1,outputArg2] = freqShiftCorrection(inputArg1,inputArg2)
%FREQSHIFTCORRECTION Summary of this function goes here
%   Detailed explanation goes here

wind_size = 2^8;
windowed_phase = windowedPhase_v(cross_talk,PRI,wind_size);

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
end

