function [] = plotFocusedWithTargets(scenario,RX,TX,targets,F,title_txt)
%PLOTWITHTARGET useful function to wrap all the plotting commands
%   [] = plotFocusedWithTargets(scenario,RX,TX,targets,F,title_txt)
F(or(F==Inf,F==-Inf)) = 0;

imagesc(scenario.grid.x_ax,scenario.grid.y_ax,abs(F.')), axis xy , 
        title(title_txt)
        xlabel('[m]'), ylabel('[m]')
        hold on,
        plot3(RX.pos(1,1),RX.pos(2,1),RX.pos(3,1),'ro'), 
        plot3(RX.pos(1,end),RX.pos(2,end),RX.pos(3,end),'go'),
        plot3(TX.pos(1,1),TX.pos(2,1),TX.pos(3,1),'kd'), 
        plot3(targets(1,:),targets(2,:),targets(3,:),'rp'), 
        legend('start drone track','end drone track','TX','targets');
        hold off
        colorbar
        colormap('jet')
end

