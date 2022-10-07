function [radar] = interpolateRC(const,radar)
%INTERPOLATERC interpolate RC wit OSF 
%   [radar] = interpolateRC(const,radar)
 

t = 0:size(radar.RC,1) -1;
tf = 0:1/const.OSF:(size(radar.RC,1) -1);

tf_men_t =tf(:) - t;
W = sinc(const.norm_B*tf_men_t);

% interpolate whole RC signal
RC_OS = W*radar.RC;
R_ax_OS = radar.R_ax(1): const.dR/const.OSF :radar.R_ax(end);

% figure
% imagesc(tau_ax,R_ax_OS,abs(RC_OS))
% title('RC interp - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

radar.RC = RC_OS;
radar.R_ax = R_ax_OS;

end

