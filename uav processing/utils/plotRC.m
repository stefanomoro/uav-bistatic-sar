function [] = plotRC(radar,scenario,mode)
%PLOTGENERICRC generic plot for showing RC matrix
%   [] = plotRC(radar,scenario,mode)

switch mode
    case 1

figure
imagesc(radar.tau_ax,radar.R_ax,abs(radar.RC))
title('RC - abs' ),xlabel("Slow time [s]"),ylabel("Range [m]")

% figure
% imagesc(radar.tau_ax,radar.R_ax,angle(radar.RC))
% title('RC - angle plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

    case 2
        figure
imagesc(radar.tau_ax,radar.R_ax,abs(radar.RC))
title('RC - abs' ),xlabel("Slow time [s]"),ylabel("Range [m]")
hold on
plot(radar.tau_ax,scenario.distance.tx_rx,'r','LineWidth',1.5),
plot(radar.tau_ax,scenario.distance.targets(1,:),'g')
legend('crosstalk','target1')
    case 3
        
        figure
imagesc(radar.tau_ax,radar.R_ax,abs(radar.RC))
title('RC - abs' ),xlabel("Slow time [s]"),ylabel("Range [m]")
hold on

plot(radar.tau_ax,scenario.distance.tx_rx),
col_letter = ['m','g','y','b'];
for idx = 1:size(scenario.distance.targets,1)
    col = col_letter(mod(idx,length(col_letter)) + 1 );
    plot(radar.tau_ax,scenario.distance.targets(idx,:),col), 
end
legend('crosstalk','targets');

end
end

