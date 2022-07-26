function [scenario] = defineFocusingGrid(const,scenario)
%DEFINEFOCUSINGGRID Summary of this function goes here
%   Detailed explanation goes here

% In this case, we choose to focus on the plane x,y  with z = c
% Set XY grid where to focus
scenario.grid.x_min = -20;
scenario.grid.x_max = 200;
scenario.grid.pho_az = 1;
scenario.grid.dx = scenario.grid.pho_az /4;
scenario.grid.x_ax =  scenario.grid.x_min: scenario.grid.dx: scenario.grid.x_max;

scenario.grid.y_min = -100;
scenario.grid.y_max = 50;
scenario.grid.dy = scenario.grid.dx;
scenario.grid.y_ax =  scenario.grid.y_min: scenario.grid.dy: scenario.grid.y_max;

scenario.grid.delta_psi_proc = const.lambda/scenario.grid.pho_az;

[scenario.grid.X,scenario.grid.Y] = ndgrid(scenario.grid.x_ax,scenario.grid.y_ax);

scenario.grid.z0 = 0;

end

