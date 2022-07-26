function [radar] = timeShiftCorrection(radar)
%TIME SHIFT CORRECTION Summary of this function goes here
%   [radar] = timeShiftCorrection(radar)

[~,cross_talk_idxs] = max(radar.RC,[],1);
% apply median filter and then moving average
cross_talk_idxs = movmean(medfilt1(cross_talk_idxs,1e3) ,1000); 

%% TIME SHIFT CORRECTION
cross_talk_shifts = cross_talk_idxs - radar.cross_talk.idxs_corr;

Nf = 2^nextpow2(size(radar.RC,1));
X = fftshift(fft(radar.RC,Nf,1),1);

f = (-Nf/2:Nf/2 - 1)/Nf;
H = exp(1i*2*pi*f(:)*cross_talk_shifts);
H(1,:) = 0;                         

RC_Dt_fixed = ifft(ifftshift(X.* H,1),Nf,1);
RC_Dt_fixed = RC_Dt_fixed(1:size(radar.RC,1),:);

radar.RC = RC_Dt_fixed; 
end

