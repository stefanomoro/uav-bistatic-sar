function [focus] = focusingTDBP_GPU_double(const,radar,scenario,RX,TX,ang_vec)
%FOCUSINGTDBP compute the focusing on the defined grid with TDBP
%   [focus] = focusingTDBP(const,radar,scenario,RX,TX,ang_vec)

Ny = length(scenario.grid.y_ax);
Nx = length(scenario.grid.x_ax);
t = radar.R_ax./physconst('LightSpeed');

% Processed Antenna aperture
focus.psi_proc = const.lambda / 2 / scenario.grid.pho_az;
focus.R_min = 50;
focus.synt_apert = 2 * tan(focus.psi_proc/2) * focus.R_min;

% Processed wavenumbers
Dk = 2*pi/scenario.grid.pho_az;

% Squint angle vectors
focus.angle_vec = ang_vec;

wbar = waitbar(0,strcat('Backprojecting n 1/',num2str(length(focus.angle_vec))));

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
X = X - reference_x;
Y = Y - reference_y;
ref = [reference_x;reference_y;0];

TX_pos = TX.pos - ref;
TX_pos_x = gpuArray(TX_pos(1,:));TX_pos_y = gpuArray(TX_pos(2,:));TX_pos_z = gpuArray(TX_pos(3,:)); 
RX_pos = RX.pos - ref;
RX_pos_x = gpuArray(RX_pos(1,:));RX_pos_y = gpuArray(RX_pos(2,:));RX_pos_z = gpuArray(RX_pos(3,:)); 
RX_speed = gpuArray(RX.speed);
X = gpuArray(X); Y = gpuArray(Y); z0 = scenario.grid.z0;
lambda = const.lambda; f0 = const.f0;
RC = gpuArray(RC);
x_ax = gpuArray(scenario.grid.x_ax);
t = gpuArray(t);
median_speed = median(RX_speed);


% Initialize vectors for the result
focus.Focused_vec = zeros(size(X,1),size(X,2),length(focus.angle_vec));
% focus.not_coh_sum = zeros(size(focus.Focused_vec),'single');
focus.SumCount = zeros(size(focus.Focused_vec));

tic
for ang_idx = 1:length(focus.angle_vec)
    waitbar(ang_idx/length(focus.angle_vec),wbar,strcat("Backprojecting n "...
        ,num2str(ang_idx),"/",num2str(length(focus.angle_vec))));
    
    psi_foc = deg2rad(focus.angle_vec(ang_idx));
    k_rx_0 = sin(psi_foc).*(2*pi/const.lambda); 
 
    S = gpuArray(zeros(Nx,Ny));
%     A = zeros(Nx,Ny,'gpuArray');
    SumCount = gpuArray(zeros(Nx,Ny));
    parfor n = 1 : size(RC,2)
        [Sn,Wn] = elementFuncTDBP(X,Y,z0,TX_pos_x(n),TX_pos_y(n),TX_pos_z(n),RX_pos_x(n),...
           RX_pos_y(n),RX_pos_z(n),lambda,Dk,RC(:,n),t,f0,k_rx_0,x_ax);
        
        
        % Give less weight to not moving positions
        speed_norm = RX_speed(n)/median_speed;
        % Count number of summations for each pixel
        SumCount = SumCount + speed_norm.*Wn;
        
        % Coherent sum over all positions along the trajectory 
        S = S + speed_norm .* Sn;
        % Inchoerent sum over all positions along the trajectory
%         A = A + abs(Sn);
    end
    waitbar(ang_idx/length(focus.angle_vec),wbar);
    
    focus.SumCount(:,:,ang_idx) = gather(SumCount);
    focus.Focused_vec(:,:,ang_idx) = gather(S);
%     focus.not_coh_sum(:,:,ang_idx) = gather(A); 
end

close(wbar)

disp (strcat("Total elaboration time: ",num2str(toc/60)," min"))
end

