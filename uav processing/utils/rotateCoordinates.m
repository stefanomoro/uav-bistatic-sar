function [pos] = rotateCoordinates(X,Y,Z,Xc,Yc,ang)
%ROTATECOORDINATES rotate coordinates around a point of a given angle
%   [pos] = rotateCoordinates(X,Y,Z,Xc,Yc,ang)

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
pos = [Xrot(:) Yrot(:) Z(:)].';
end

