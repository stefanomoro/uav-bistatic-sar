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
Dk = single(2*pi/scenario.grid.pho_az);

% Sqint angle vectors
focus.angle_vec = (-35:5:35).';

% Copy variables for optimizing parfor
%reset(gpuDevice)
idxs = t >= 0;
t = t(idxs);
RC = radar.RC(idxs,:);

TX_pos = single(TX.pos);
TX_pos_x = gpuArray(TX_pos(1,:));TX_pos_y = gpuArray(TX_pos(2,:));TX_pos_z = gpuArray(TX_pos(3,:)); 
RX_pos = single(RX.pos);
RX_pos_x = gpuArray(RX_pos(1,:));RX_pos_y = gpuArray(RX_pos(2,:));RX_pos_z = gpuArray(RX_pos(3,:)); 
RX_speed = gpuArray(single(RX.speed));
X = gpuArray(single(scenario.grid.X)); Y = gpuArray(single(scenario.grid.Y)); z0 = single(scenario.grid.z0);
lambda = single(const.lambda); f0 = single(const.f0);
RC = gpuArray(single(RC));
t = gpuArray(single(t));
%max_speed = max(single(RX_speed));


% Initialize vectors for the result
% focus.Focused_vec = zeros(size(X,1),size(X,2),length(focus.angle_vec),'single');
% focus.not_coh_sum = zeros(size(focus.Focused_vec),'single');
% focus.SumCount = zeros(size(focus.Focused_vec),'single');

tic
% for ang_idx = 1:length(focus.angle_vec)   
    
    k_rx_0_vec = single(sin(deg2rad(focus.angle_vec)).*(2*pi/const.lambda)); 
 
      
    [S,SumCount] = cudaFocusing(X,Y,z0,TX_pos_x,TX_pos_y,TX_pos_z,RX_pos_x,...
        RX_pos_y,RX_pos_z,lambda,Dk,RC,t,f0,k_rx_0_vec);
    wait(gpuDevice);
    
    focus.SumCount = gather(SumCount);
    focus.Focused_vec = gather(S);


disp (strcat("Total elaboration time: ",num2str(toc/60)," min"))
end

