%% =========================================================================== FOCUSING

%% FOCUSING BY TIME DOMAIN BACK PROJECTION (TDBP)
% In this case, we choose to focus on the plane x,y  with z = c
%% Set XY grid where to focus
x_min = -20;
x_max = 200;
pho_az = 1;
dx = pho_az*0.4;
x_ax =  x_min:dx:x_max;

y_min = -100;
y_max = 50;
dy = dx;
y_ax =  y_min:dy:y_max;

delta_psi_proc = lambda/pho_az;

[X,Y] = ndgrid(x_ax,y_ax);

z0 = 0;
%% Focusing (Wn form wave numbers)


Ny = length(y_ax);
Nx = length(x_ax);
t = R_ax./c;

% Wavenumber resolution
Dk = 2*pi/pho_az;

%Compute path angle
psi_path = atan( (RX_pos(2,end)-RX_pos(2,1)) / (RX_pos(1,end)-RX_pos(1,1)) );

%Compute pointing angle 
psi_point = psi_path+deg2rad(90); % for case with track not rotated

%Focalization wave number
angle_vec = [0];
Focused_vec = cell(size(angle_vec));
tic
for ang_idx = 1:length(angle_vec)
    wbar = waitbar(0,strcat('Backprojecting n',num2str(ang_idx),"/",num2str(length(angle_vec))));
    psi_foc = deg2rad(angle_vec(ang_idx));
    k_rx_0 = sin(psi_foc).*(2*pi/lambda); 
 
    S = zeros(Nx,Ny);
    Sn = S; A =S;
    figure
    for n = 1:N_PRI
        waitbar(n/N_PRI,wbar)

        % Distance 
        R_tx = sqrt((TX_pos(1,n)-X).^2 + (TX_pos(2,n)-Y).^2  + (TX_pos(3,n)-z0).^2); %  Range distances from the tx antenna [m]
        R_rx = sqrt((RX_pos(1,n)-X).^2 + (RX_pos(2,n)-Y).^2  + (RX_pos(3,n)-z0).^2); %  Range distances from the rx antenna [m]
        distance = R_tx+R_rx; %Total Tx-target-Rx distance [m]
        delay = distance./c;    %Delay

        %Compute target wave number
        R = sqrt((RX_pos(1,n)-X).^2 + (RX_pos(2,n)-Y).^2);
        psi = asin((Y-RX_pos(2,n))./R);
        k_rx = sin(psi).*(2*pi/lambda);

        %Weight function
%         Wn = rectpuls((k_rx - k_rx_0)./Dk);
        sigma = Dk/2;
        gauss = @(x) 1/(sigma*sqrt(2*pi)) * exp(-0.5*((x)./sigma).^2); 
        Wn = gauss(k_rx - k_rx_0);
        
        cut = find(x_ax>RX_pos(1,n)); %Cut the otherside lobe
        cut = cut(1);
        Wn(1:cut,:) = zeros(size(Wn(1:cut,:)));

        % Backprojection of data from a single Radar position 
        Sn = Wn.*interp1(t,RC(:,n),delay).*exp(+1i*2*pi*f0*delay);
        if mod(n,100) == 0 
            imagesc(abs(Sn))
            drawnow
            pause(.5)
        end
        % Coherent sum over all positions along the trajectory 
        S = S + Sn;
        % Inchoerent sum over all positions along the trajectory (choerent sum)
        A = A + abs(Sn);
    end
    close(wbar)


     Focus = (S./A)';
    %Focus = S';
    %Focus = A';
    Focused_vec{ang_idx} = Focus;
end
disp (strcat("Total elaboration time: ",num2str(toc/3600)," h"))