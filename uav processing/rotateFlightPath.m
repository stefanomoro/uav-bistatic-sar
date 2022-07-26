function [RX,TX,targets,drone] = rotateFlightPath(RX,TX,targets,drone)
%ROTATE FLIGHT PATH rotate processed flight path to make it parallel with x
%   [RX,TX,targets,drone] = rotateFlightPath(RX,TX,targets,drone)

%% ========================================= Rotation of the path
% Compute the angle(to have y = const during swath)
ang = deg2rad(90)+atan((RX.pos(2,end)-RX.pos(2,1)) /  (RX.pos(1,end)-RX.pos(1,1)));

% Center of rotation
Xc = RX.pos(1,1);  % X cener of rotation
Yc = RX.pos(2,1);  % Y center of rotation

%% All drone path
X = drone.x;
Y = drone.y;
Z = drone.alt;

% Shift X/Y to the rotation center
Xshift = X - Xc;        
Yshift = Y - Yc;
% Rotate the coordinates
Xsrot =  Xshift*cos(ang) + Yshift*sin(ang);
Ysrot = -Xshift*sin(ang) + Yshift*cos(ang);
% Shift the rotated coordinates back to the original reference center
Xrot = Xsrot + Xc;
Yrot = Ysrot + Yc;
% Assemble position matrix
drone.x = Xrot;
drone.y = Yrot;
drone.alt = Z;


%% RX
X = RX.pos(1,:);
Y = RX.pos(2,:);
Z = RX.pos(3,:);

% Shift X/Y to the rotation center
Xshift = X - Xc;        
Yshift = Y - Yc;
% Rotate the coordinates
Xsrot =  Xshift*cos(ang) + Yshift*sin(ang);
Ysrot = -Xshift*sin(ang) + Yshift*cos(ang);
% Shift the rotated coordinates back to the original reference center
Xrot = Xsrot + Xc;
Yrot = Ysrot + Yc;
% Assemble position matrix
RX.pos = [Xrot;Yrot;Z];

%% TX
X = TX.pos(1,:);
Y = TX.pos(2,:);
Z = TX.pos(3,:);

% Shift X/Y to the rotation center
Xshift = X - Xc;        
Yshift = Y - Yc;
% Rotate the coordinates
Xsrot =  Xshift*cos(ang) + Yshift*sin(ang);
Ysrot = -Xshift*sin(ang) + Yshift*cos(ang);
% Shift the rotated coordinates back to the original reference center
Xrot = Xsrot + Xc;
Yrot = Ysrot + Yc;
% Assemble position matrix
TX.pos = [Xrot;Yrot;Z];

%% Cars
X = targets.cars(1,:);
Y = targets.cars(2,:);
Z = targets.cars(3,:);

% Shift X/Y to the rotation center
Xshift = X - Xc;        
Yshift = Y - Yc;
% Rotate the coordinates
Xsrot =  Xshift*cos(ang) + Yshift*sin(ang);
Ysrot = -Xshift*sin(ang) + Yshift*cos(ang);
% Shift the rotated coordinates back to the original reference center
Xrot = Xsrot + Xc;
Yrot = Ysrot + Yc;
% Assemble position matrix
targets.cars = [Xrot;Yrot;Z];

%% Humans
X = targets.humans(1,:);
Y = targets.humans(2,:);
Z = targets.humans(3,:);

% Shift X/Y to the rotation center
Xshift = X - Xc;        
Yshift = Y - Yc;
% Rotate the coordinates
Xsrot =  Xshift*cos(ang) + Yshift*sin(ang);
Ysrot = -Xshift*sin(ang) + Yshift*cos(ang);
% Shift the rotated coordinates back to the original reference center
Xrot = Xsrot + Xc;
Yrot = Ysrot + Yc;
% Assemble position matrix
targets.humans = [Xrot;Yrot;Z];


end