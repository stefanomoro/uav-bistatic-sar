function [] = plotWithTargets(x_ax,y_ax,F,title_txt,RX_pos,TX_pos,cars,humans)
%PLOTWITHTARGET useful function to wrap all the plotting commands
%   Detailed explanation goes here
imagesc(x_ax,y_ax,abs(F.')), axis xy , 
        title(title_txt)
        xlabel('[m]'), ylabel('[m]')
        hold on,
        plot3(RX_pos(1,1),RX_pos(2,1),RX_pos(3,1),'ro'), 
        plot3(RX_pos(1,end),RX_pos(2,end),RX_pos(3,end),'go'),
        plot3(TX_pos(1,1),TX_pos(2,1),TX_pos(3,1),'kd'), plot3(cars(1,:),cars(2,:),cars(3,:),'rp'), plot3(humans(1,:),humans(2,:),humans(3,:),'yh')
        legend('start drone track','end drone track','TX','Cars','Humans');
        hold off
end

