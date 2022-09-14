function [RX,TX,targets,drone] = rotateFlightPath(RX,TX,targets,drone)
%ROTATE FLIGHT PATH rotate processed flight path to make it parallel with x
%   [RX,TX,targets,drone] = rotateFlightPath(RX,TX,targets,drone)

% linear regression of path
x = RX.pos(1,1:1e3:end)';
X = [ones(length(x),1) x];
y = RX.pos(2,1:1e3:end)';
b = X\y;
ang = deg2rad(90) + atan(b(2));
% Compute the angle(to have y = const during swath)
% ang = deg2rad(90)+atan((RX.pos(2,end)-RX.pos(2,1)) /  (RX.pos(1,end)-RX.pos(1,1)));

% Center of rotation
Xc = RX.pos(1,1);  % X cener of rotation
Yc = RX.pos(2,1);  % Y center of rotation

%% All drone path
if (isfield(drone,"rx"))
    [pos] = rotateCoordinates(drone.rx.utm_x,drone.rx.utm_y,drone.rx.alt,Xc,Yc,ang);
else
% old setup mode doesn't have utm coordinates
    [pos] = rotateCoordinates(drone.x,drone.y,drone.alt,Xc,Yc,ang);
end
% omologate old and new steup mode
drone.rx.x = pos(1,:).';
drone.rx.y = pos(2,:).';
drone.rx.alt = pos(3,:).';

if(isfield(drone,"tx"))
    [pos] = rotateCoordinates(drone.tx.utm_x,drone.tx.utm_y,drone.tx.alt,Xc,Yc,ang);
    drone.tx.x = pos(1,:).';
    drone.tx.y = pos(2,:).';
    drone.tx.alt = pos(3,:).';
end

%% RX
RX.pos = rotateCoordinates(RX.pos(1,:),RX.pos(2,:),RX.pos(3,:),Xc,Yc,ang);
%% TX
TX.pos = rotateCoordinates(TX.pos(1,:),TX.pos(2,:),TX.pos(3,:),Xc,Yc,ang);

%% Targets
targets = rotateCoordinates(targets(1,:),targets(2,:),targets(3,:),Xc,Yc,ang);
end