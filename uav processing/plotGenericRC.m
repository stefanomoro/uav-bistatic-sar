function [] = plotGenericRC(radar)
%PLOTGENERICRC generic plot for showing RC matrix

figure
imagesc(radar.tau_ax,radar.R_ax,abs(radar.RC))
title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

figure
imagesc(radar.tau_ax,radar.R_ax,angle(radar.RC))
title('RC - angle plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

end

