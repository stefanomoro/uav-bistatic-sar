function [windowed_phase] = windowedPhase(cross_talk,PRI,win_size)
%WINDOWED_PHASE Summary of this function goes here
%   Detailed explanation goes here

%number of shifts of the window over the peak
N_cycle = floor(length(cross_talk)/win_size) * 2 - 1;
Df = zeros(N_cycle,1);
N_fft = 8*win_size;
freq_ax = (-N_fft/2:N_fft/2-1)/N_fft * 1/PRI;
win_idx = 1 : win_size;

% check_win_size = zeros(N_cycle,1);
for i = 1:N_cycle
    aa = cross_talk(win_idx);
    ff = fftshift(fft(aa,N_fft));
    [~,ii] = max(abs(ff));
    Df(i) = freq_ax(ii);
    win_idx = win_idx + win_size/2;
    
%     check_win_size(i) = max(abs(ff))/sum(abs(aa));
end

% figure,plot(Df)
% title('Df from fft')
% xlabel('Slow time [samples]')


% figure,plot(check_win_size)
% title('check window size')
% xlabel('Slow time [samples]')

%% INTEGRATE Df to get phase
phase_0 = angle(cross_talk(1));

% computed_phase time axis is in the centroid of fft windows
t = (1: length(Df))*win_size/2*PRI;
windowed_phase = cumsum(2*pi*Df*win_size/2 * PRI);
windowed_phase = interp1(t,windowed_phase,(0:length(cross_talk)-1)*PRI,'pchip');
% start from same phase of crosstalk
windowed_phase = windowed_phase - windowed_phase(1) + phase_0;

end

