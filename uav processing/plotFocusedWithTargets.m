function [] = plotFocusedWithTargets(scenario,RX,TX,targets,F,title_txt)
%PLOTWITHTARGET useful function to wrap all the plotting commands
%   [] = plotFocusedWithTargets(scenario,RX,TX,targets,F,title_txt)
imagesc(scenario.grid.x_ax,scenario.grid.y_ax,abs(F.')), axis xy , 
        title(title_txt)
        xlabel('[m]'), ylabel('[m]')
        hold on,
        plot3(RX.pos(1,1),RX.pos(2,1),RX.pos(3,1),'ro'), 
        plot3(RX.pos(1,end),RX.pos(2,end),RX.pos(3,end),'go'),
        plot3(TX.pos(1,1),TX.pos(2,1),TX.pos(3,1),'kd'), 
        plot3(targets.cars(1,:),targets.cars(2,:),targets.cars(3,:),'rp'), 
        plot3(targets.humans(1,:),targets.humans(2,:),targets.humans(3,:),'yh')
        legend('start drone track','end drone track','TX','Cars','Humans');
        hold off
end

