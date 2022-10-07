function [focus] = focusingTDBP(const,radar,scenario,RX,TX)
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

%Compute path angle
% psi_path = atan( (RX.pos(2,end)-RX.pos(2,1)) / (RX.pos(1,end)-RX.pos(1,1)) );

%Sqint angle vectors
focus.angle_vec = [0];%-35:5:35;
%Initialize vectors for the result
focus.Focused_vec = cell(size(focus.angle_vec));
focus.not_coh_sum = cell(size(focus.angle_vec));

wbar = waitbar(0,strcat('Backprojecting n 1/',num2str(length(focus.angle_vec))));

% copy variables for optimizing parfor
TX_pos_x = TX.pos(1,:);TX_pos_y = TX.pos(2,:);TX_pos_z = TX.pos(3,:); 
RX_pos_x = RX.pos(1,:);RX_pos_y = RX.pos(2,:);RX_pos_z = RX.pos(3,:); 
X = scenario.grid.X; Y = scenario.grid.Y; z0 = scenario.grid.z0;
lambda = const.lambda; f0 = const.f0;
RC = radar.RC;
x_ax = scenario.grid.x_ax;

tic
for ang_idx = 1:length(focus.angle_vec)
    waitbar(ang_idx/length(focus.angle_vec),wbar,strcat('Backprojecting n '...
        ,num2str(ang_idx),"/",num2str(length(focus.angle_vec))));
    psi_foc = deg2rad(focus.angle_vec(ang_idx));
    k_rx_0 = sin(psi_foc).*(2*pi/const.lambda); 
 
    S = zeros(Nx,Ny);
    A = zeros(Nx,Ny);
    
    parfor n = 1:size(RC,2)
        
        
        R_tx = sqrt((TX_pos_x(n)-X).^2 + (TX_pos_y(n)-Y).^2  + ...
            (TX_pos_z(n)-z0).^2);                                           %Range distances from the tx antenna [m]
        
        R_rx = sqrt((RX_pos_x(n)-X).^2 + (RX_pos_y(n)-Y).^2  + ...
            (RX_pos_z(n)-z0).^2);                                           %Range distances from the rx antenna [m]
        distance = R_tx+R_rx;                                               %Total Tx-target-Rx distance [m]
        delay = distance./physconst('LightSpeed');

        %Compute target wave number
        R = sqrt((RX_pos_x(n)-X).^2 + (RX_pos_y(n)-Y).^2);
        psi = asin((Y-RX_pos_y(n))./R);
        k_rx = sin(psi).*(2*pi/lambda);

        %Weight function
%         Wn = rectpuls((k_rx - k_rx_0)./psi_proc);
        sigma = Dk/2;
        Wn = gaussActivFunc(k_rx - k_rx_0,sigma);
        
        cut = find(x_ax>RX_pos_x(n));                                       %Cut the back-lobe
        cut = cut(1);
        Wn(1:cut,:) = zeros(size(Wn(1:cut,:)));

        % Backprojection of data from a single Radar position 
        Sn = Wn.*interp1(t,RC(:,n),delay).*exp(+1i*2*pi*f0*delay);

        % Coherent sum over all positions along the trajectory 
        S = S + Sn;
        % Inchoerent sum over all positions along the trajectory
        A = A + abs(Sn);
    end
    waitbar(ang_idx/length(focus.angle_vec),wbar);
    


% 	Focus = (S./A)';
    Focus = S;
%     Focus = A';
    focus.Focused_vec{ang_idx} = Focus;
    focus.not_coh_sum{ang_idx} = A; 
end

close(wbar)

disp (strcat("Total elaboration time: ",num2str(toc/60)," min"))
end

