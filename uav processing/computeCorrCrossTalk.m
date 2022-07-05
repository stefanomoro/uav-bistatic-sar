function [radar] = computeCorrCrossTalk(const,radar,scenario)
%UNTITLED compute correct cross-talk idx and phase from distance
%   [radar] = computeCorrCrossTalk(const,radar,scenario)

% define the idxs where the cross talk should be

idxs_corr = zeros(size(radar.tau_ax));
for n = 1:length(radar.tau_ax)
    idxs = find(radar.R_ax > scenario.distance.tx_rx(n));
    idxs_corr(n) = idxs(1);
end

radar.cross_talk.idxs_corr = movmean(idxs_corr,1e3);

% Compute correct cross talk phase from distance

c = physconst('Lightspeed');
phase_corr = -2*pi*const.f0*(scenario.distance.tx_rx/c); % - 2*pi*f0*(2*speed_distance_tx_rx/lambda);
phase_corr = phase_corr(:).';

% move first phase in 2pi interval
phase_0 = mod(phase_corr(1),2*pi)-pi;

phase_corr = phase_corr - phase_corr(1) + phase_0;

% figure, plot(tau_ax, unwrap(angle(cross_talk))),
% hold on, plot(tau_ax, phase_corr),  xlabel('\tau [s]'), ylabel('\phi [rad]')
% legend('Crosstalk phase','Correct phase')
% saveas(gcf,'Figures/phase_cross_corr.png')

radar.cross_talk.phase_corr = phase_corr;

end

