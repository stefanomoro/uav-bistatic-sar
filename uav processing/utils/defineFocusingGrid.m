function [scenario] = defineFocusingGrid(const,scenario,RX)
%DEFINEFOCUSINGGRID define the focusing grid for the TDBP
%   [scenario] = defineFocusingGrid(const,scenario,RX)

% In this case, we choose to focus on the plane x,y  with z = zero
scenario.grid.pho_az = 1;

% Set XY grid where to focus
scenario.grid.x_min = min(RX.pos(1,:)) - 5;
scenario.grid.x_max = max(RX.pos(1,:)) + 150;

scenario.grid.dx = scenario.grid.pho_az /4;
scenario.grid.x_ax =  scenario.grid.x_min: scenario.grid.dx: scenario.grid.x_max;

scenario.grid.y_min = min(RX.pos(2,:)) - 5;
scenario.grid.y_max = max(RX.pos(2,:)) + 5;
scenario.grid.dy = scenario.grid.dx;
scenario.grid.y_ax =  scenario.grid.y_min: scenario.grid.dy: scenario.grid.y_max;

scenario.grid.delta_psi_proc = const.lambda/scenario.grid.pho_az;

[scenario.grid.X,scenario.grid.Y] = ndgrid(scenario.grid.x_ax,scenario.grid.y_ax);

scenario.grid.z0 = 0;

end

