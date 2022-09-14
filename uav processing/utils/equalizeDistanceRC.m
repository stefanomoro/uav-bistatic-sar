function [focus] = equalizeDistanceRC(scenario,TX,focus)
%EQUALIZEDISTANCERC multiply focused image for an equalizer based on distance 
%   [Focus_eq] = equalizeDistanceRC(scenario,TX,focus)

X = scenario.grid.X; 
Y = scenario.grid.Y;
z0 = scenario.grid.z0;
R = zeros(size(X));
% assume transmitter static
TX_pos = TX.pos(:,1);
Focus_eq = zeros(size(focus.Focused_vec));
for ang_idx = 1 : length(focus.angle_vec)
    psi_foc = deg2rad(focus.angle_vec(ang_idx));
    Focus = focus.Focused_vec(:,:,ang_idx);
    for i = 1:size(Focus,1)
       for j = 1:size(Focus,2)
           R_tx = sqrt((TX_pos(1)-X(i,j)).^2 + (TX_pos(2)-Y(i,j)).^2  + (TX_pos(3)-z0).^2);
           R_rx = sqrt((TX_pos(1)-X(i,j)).^2 + (TX_pos(3)-z0).^2) / sind(90-psi_foc);
           R(i,j) = (R_tx + R_rx)/2;
        end
    end

Focus_eq(:,:,ang_idx) = Focus.*R;
end
focus.Focus_eq = Focus_eq;

end

