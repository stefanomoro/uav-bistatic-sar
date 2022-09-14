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

% Sqint angle vectors
focus.angle_vec = -35:5:35;

wbar = waitbar(0,strcat('Backprojecting n 1/',num2str(length(focus.angle_vec))));

% Copy variables for optimizing parfor
TX_pos_x = gpuArray(TX.pos(1,:));TX_pos_y = gpuArray(TX.pos(2,:));TX_pos_z = gpuArray(TX.pos(3,:)); 
RX_pos_x = gpuArray(RX.pos(1,:));RX_pos_y = gpuArray(RX.pos(2,:));RX_pos_z = gpuArray(RX.pos(3,:)); 
RX_speed = gpuArray(RX.speed);
X = gpuArray(scenario.grid.X); Y = gpuArray(scenario.grid.Y); z0 = gpuArray(scenario.grid.z0);
lambda = const.lambda; f0 = const.f0;
RC = gpuArray(radar.RC);
x_ax = scenario.grid.x_ax;
max_speed = max(RX_speed);

% Initialize vectors for the result
focus.Focused_vec = zeros(size(X,1),size(X,2),length(focus.angle_vec),'gpuArray');
focus.not_coh_sum = zeros(size(focus.Focused_vec),'gpuArray');
focus.SumCount = zeros(size(focus.Focused_vec),'gpuArray');

tic
for ang_idx = 1:length(focus.angle_vec)
    waitbar(ang_idx/length(focus.angle_vec),wbar,strcat("Backprojecting n "...
        ,num2str(ang_idx),"/",num2str(length(focus.angle_vec))));
    
    psi_foc = deg2rad(focus.angle_vec(ang_idx));
    k_rx_0 = sin(psi_foc).*(2*pi/const.lambda); 
 
    S = zeros(Nx,Ny,'gpuArray');
    A = zeros(Nx,Ny,'gpuArray');
    SumCount = gpuArray(zeros(Nx,Ny));
    parfor n = 1 : size(RC,2)
        [Sn,Wn] = elementFuncTDBP(X,Y,z0,TX_pos_x(n),TX_pos_y(n),TX_pos_z(n),RX_pos_x(n),...
            RX_pos_y(n),RX_pos_z(n),lambda,Dk,RC(:,n),t,f0,k_rx_0,x_ax);
        
        %Sn = gather(Sn);
        % Give less weight to not moving positions
        speed_norm = RX_speed(n)/max_speed;
        % Count number of summations for each pixel
        SumCount = SumCount + speed_norm.*Wn;
        
        % Coherent sum over all positions along the trajectory 
        S = S + speed_norm .* Sn;
        % Inchoerent sum over all positions along the trajectory
        A = A + abs(Sn);
    end
    waitbar(ang_idx/length(focus.angle_vec),wbar);
    
    focus.SumCount(:,:,ang_idx) = SumCount;
    focus.Focused_vec(:,:,ang_idx) = S;
    focus.not_coh_sum(:,:,ang_idx) = A; 
end

close(wbar)

disp (strcat("Total elaboration time: ",num2str(toc/60)," min"))
end

