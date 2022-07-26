function [] = plotRC(radar,scenario,mode)
%PLOTGENERICRC generic plot for showing RC matrix
%   [] = plotRC(radar,scenario,mode)

switch mode
    case 1

figure
imagesc(radar.tau_ax,radar.R_ax,abs(radar.RC))
title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

figure
imagesc(radar.tau_ax,radar.R_ax,angle(radar.RC))
title('RC - angle plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

    case 2
        figure
imagesc(radar.tau_ax,radar.R_ax,abs(radar.RC))
title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")
hold on
plot(radar.tau_ax,scenario.distance.tx_rx,'r','LineWidth',1.5),
plot(radar.tau_ax,scenario.distance.cars(1,:),'g')
legend('crosstalk','car')
    case 3
        
        figure
imagesc(radar.tau_ax,radar.R_ax,abs(radar.RC))
title('RC - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")
hold on

plot(radar.tau_ax,scenario.distance.tx_rx),
plot(radar.tau_ax,scenario.distance.cars(1,:),'m'), 
plot(radar.tau_ax,scenario.distance.cars(2,:),'g'), 
plot(radar.tau_ax,scenario.distance.humans(1,:),'y'), 
plot(radar.tau_ax,scenario.distance.humans(2,:),'m'), 
plot(radar.tau_ax,scenario.distance.humans(3,:),'b')
legend('crosstalk','car','car 2','human1','human2','human3');
end
end

