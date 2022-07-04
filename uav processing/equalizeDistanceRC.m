function [Focus_eq] = equalizeDistanceRC(Focus,x_ax,y_ax,z0,TX_pos,RX_pos,psi_foc)
%EQUALIZEDISTANCERC multiply focused image for an equalizer based on distance 
%   Detailed explanation goes here

[X,Y] = ndgrid(x_ax,y_ax);
R = zeros(size(X));

TX_pos = TX_pos(:,1);
for i = 1:size(Focus,1)
   for j = 1:size(Focus,2)
       R_tx = sqrt((TX_pos(1)-X(i,j)).^2 + (TX_pos(2)-Y(i,j)).^2  + (TX_pos(3)-z0).^2);
       R_rx = sqrt((TX_pos(1)-X(i,j)).^2 + (TX_pos(3)-z0).^2) / sind(90-psi_foc);
       R(i,j) = (R_tx + R_rx)/2;
    end
end

Focus_eq = Focus.*R;

end

