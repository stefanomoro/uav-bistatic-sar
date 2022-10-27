function [focus] = focusingTDBP_CUDA(const,radar,scenario,RX,TX,ang_vec)
%FOCUSINGTDBP compute the focusing on the defined grid with TDBP
%   [focus] = focusingTDBP_CUDA(const,radar,scenario,RX,TX,ang_vec)


t = radar.R_ax./physconst('LightSpeed');

% Processed Antenna aperture
focus.psi_proc = const.lambda / 2 / scenario.grid.pho_az;
focus.R_min = 50;
focus.synt_apert = 2 * tan(focus.psi_proc/2) * focus.R_min;

% Processed wavenumbers
Dk = single(2*pi/scenario.grid.pho_az);

% Sqint angle vectors
focus.angle_vec = ang_vec(:);

% Copy variables for optimizing parfor
%reset(gpuDevice)
idxs = t >= 0;
t = t(idxs);
RC = radar.RC(idxs,:);


% remove mean value from grid 
X = scenario.grid.X;
Y = scenario.grid.Y;
reference_x = mean(X(:,1));
reference_y = mean(Y(1,:));
X = single(X - reference_x);
Y = single(Y - reference_y);
ref = [reference_x;reference_y;0];

TX_pos = single(TX.pos - ref);
TX_pos_x = gpuArray(TX_pos(1,:));TX_pos_y = gpuArray(TX_pos(2,:));TX_pos_z = gpuArray(TX_pos(3,:)); 
RX_pos = single(RX.pos - ref);
RX_pos_x = gpuArray(RX_pos(1,:));RX_pos_y = gpuArray(RX_pos(2,:));RX_pos_z = gpuArray(RX_pos(3,:)); 
RX_speed = single(RX.speed);
X = gpuArray(X); Y = gpuArray(Y); z0 = single(scenario.grid.z0);
lambda = single(const.lambda); f0 = single(const.f0);
RC = gpuArray(single(RC));
t = gpuArray(single(t));
median_speed = median(RX_speed);



tic
    
k_rx_0_vec = single(sin(deg2rad(focus.angle_vec)).*(2*pi/const.lambda)); 
  
[S,SumCount] = cudaFocusing(X,Y,z0,TX_pos_x,TX_pos_y,TX_pos_z,RX_pos_x,...
    RX_pos_y,RX_pos_z,lambda,Dk,RC,t,f0,k_rx_0_vec,RX_speed,median_speed);
wait(gpuDevice);

focus.SumCount = gather(SumCount);
focus.Focused_vec = gather(S);


disp (strcat("Total elaboration time: ",num2str(toc/60)," min"))
end

