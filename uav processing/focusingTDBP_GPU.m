function [focus] = focusingTDBP_GPU(const,radar,scenario,RX,TX)
%FOCUSINGTDBP compute the focusing on the defined grid with TDBP
%   [focus] = focusingTDBP(const,radar,scenario,RX,TX)

Ny = length(scenario.grid.y_ax);
Nx = length(scenario.grid.x_ax);
t = radar.R_ax./physconst('LightSpeed');

% Processed Antenna aperture
focus.psi_proc = const.lambda / 2 / scenario.grid.pho_az;
focus.R_min = 50;
focus.synt_apert = 2 * tan(focus.psi_proc/2) * focus.R_min;
% Processed wavenumbers
Dk = 2*pi/scenario.grid.pho_az;

%Sqint angle vectors
focus.angle_vec = [0];%-35:5:35;
%Initialize vectors for the result
focus.Focused_vec = cell(size(focus.angle_vec));
focus.not_coh_sum = cell(size(focus.angle_vec));

wbar = waitbar(0,strcat('Backprojecting n 1/',num2str(length(focus.angle_vec))));

% copy variables for optimizing parfor
TX_pos_x = TX.pos(1,const.focus_PRI_cut);TX_pos_y = TX.pos(2,const.focus_PRI_cut);TX_pos_z = TX.pos(3,const.focus_PRI_cut); 
RX_pos_x = RX.pos(1,const.focus_PRI_cut);RX_pos_y = RX.pos(2,const.focus_PRI_cut);RX_pos_z = RX.pos(3,const.focus_PRI_cut); 
X = gpuarray(scenario.grid.X); gpuarray(Y = scenario.grid.Y); z0 = gpuarray(scenario.grid.z0);
lambda = const.lambda; f0 = const.f0;
RC = radar.RC(const.focus_PRI_cut);
x_ax = scenario.grid.x_ax;

tic
for ang_idx = 1:length(focus.angle_vec)
    waitbar(ang_idx/length(focus.angle_vec),wbar,strcat('Backprojecting n '...
        ,num2str(ang_idx),"/",num2str(length(focus.angle_vec))));
    psi_foc = deg2rad(focus.angle_vec(ang_idx));
    k_rx_0 = sin(psi_foc).*(2*pi/const.lambda); 
 
    S = zeros(Nx,Ny);
    A = zeros(Nx,Ny);
    
        Sn = arrayfun( @elementFuncTDBP, ...
            X,Y,z0,TX_pos_x,TX_pos_y,TX_pos_z,RX_pos_x,RX_pos_y,RX_pos_z,lambda,Dk,RC,t,f0,k_rx_0);
        Sn = gather(Sn);

        % Coherent sum over all positions along the trajectory 
        S = S + Sn;
        % Inchoerent sum over all positions along the trajectory
        A = A + abs(Sn);
    
    waitbar(ang_idx/length(focus.angle_vec),wbar);
    


% 	Focus = (S./A)';
    Focus = S;
%     Focus = A';
    focus.Focused_vec{ang_idx} = Focus;
    focus.not_coh_sum{ang_idx} = A; 
end

close(wbar)

disp (strcat("Total elaboration time: ",num2str(toc/3600)," h"))
end

