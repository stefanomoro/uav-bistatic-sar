function [] = plotDroneTrack(drone,RX,TX,targets)
%PLOTDRONETRACK plot the track of the drone with targets and TX
%   [] = plotDroneTrack(drone,RX,TX,targets)
figure, plot3(drone.rx.x,drone.rx.y,drone.rx.alt), title('Flight path'), xlabel('x'), ylabel('y'), zlabel('z')

hold on,
plot3(RX.pos(1,1),RX.pos(2,1),RX.pos(3,1),'ro'), 
plot3(RX.pos(1,end),RX.pos(2,end),RX.pos(3,end),'go'),
plot3(RX.pos(1,:),RX.pos(2,:),RX.pos(3,:),'b','LineWidth',1.5), 
plot3(TX.pos(1,:),TX.pos(2,:),TX.pos(3,:),'dk'), 
plot3(targets(1,:),targets(2,:),targets(3,:),'p'), 

legend('RX path complete','start','end','processed track','TX','Targets');grid on
end

