function [focus] = focusingTDBP_CUDA_double(const,radar,scenario,RX,TX,ang_vec)
%FOCUSINGTDBP compute the focusing on the defined grid with TDBP
%   [focus] = focusingTDBP_CUDA(const,radar,scenario,RX,TX,ang_vec)

t = radar.R_ax./physconst('LightSpeed');

% Processed Antenna aperture
focus.psi_proc = const.lambda / 2 / scenario.grid.pho_az;
focus.R_min = 50;
focus.synt_apert = 2 * tan(focus.psi_proc/2) * focus.R_min;

% Processed wavenumbers
Dk = 2*pi/scenario.grid.pho_az;

% Sqint angle vectors
focus.angle_vec = ang_vec(:);

% Copy variables for optimizing parfor
%reset(gpuDevice)
idxs = t >= 0;
t = t(idxs);
RC = radar.RC(idxs,:);

TX_pos = TX.pos;
TX_pos_x = gpuArray(TX_pos(1,:));TX_pos_y = gpuArray(TX_pos(2,:));TX_pos_z = gpuArray(TX_pos(3,:)); 
RX_pos = RX.pos;
RX_pos_x = gpuArray(RX_pos(1,:));RX_pos_y = gpuArray(RX_pos(2,:));RX_pos_z = gpuArray(RX_pos(3,:)); 
RX_speed = RX.speed;
X = gpuArray(scenario.grid.X); Y = gpuArray(scenario.grid.Y); z0 = scenario.grid.z0;
lambda = const.lambda; f0 = const.f0;
RC = gpuArray(RC);
t = gpuArray(t);
median_speed = median(RX_speed);



tic
    
k_rx_0_vec = sin(deg2rad(focus.angle_vec)).*(2*pi/const.lambda); 
  
[S,SumCount] = cudaFocusing_double(X,Y,z0,TX_pos_x,TX_pos_y,TX_pos_z,RX_pos_x,...
    RX_pos_y,RX_pos_z,lambda,Dk,RC,t,f0,k_rx_0_vec,RX_speed,median_speed);
wait(gpuDevice);

focus.SumCount = gather(SumCount);
focus.Focused_vec = gather(S);


disp (strcat("Total elaboration time: ",num2str(toc/60)," min"))
end

