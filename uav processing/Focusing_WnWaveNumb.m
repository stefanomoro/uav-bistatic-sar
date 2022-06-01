%% =========================================================================== FOCUSING

%% FOCUSING BY TIME DOMAIN BACK PROJECTION (TDBP)
% In this case, we choose to focus on the plane x,y  with z = c

%% Focusing (Wn form wave numbers)

wbar = waitbar(0,'Backprojecting');
Ny = length(y_ax);
Nx = length(x_ax);
S = zeros(Nx,Ny);
A = zeros(Nx,Ny);
t = R_ax./c;

% Wavenumber resolution
Dk = 2*pi/pho_az;

%Compute path angle
psi_path = atan( (RX_pos(2,end)-RX_pos(2,1)) / (RX_pos(1,end)-RX_pos(1,1)) );

%Compute pointing angle 
psi_point = psi_path+deg2rad(90);

%Focalization wave number
psi_foc = deg2rad(0); % = psi_point is better
k_rx_0 = sin(psi_foc).*(2*pi/lambda); 

clear Sn    
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
    Wn = rectpuls((k_rx - k_rx_0)./Dk);
    cut = find(x_ax>RX_pos(1,n)); %Cut the otherside lobe
    cut = cut(1);
    Wn(1:cut,:) = zeros(size(Wn(1:cut,:)));
    
    % Backprojection of data from a single Radar position 
    Sn = Wn.*interp1(t,RC(:,n),delay).*exp(+1i*2*pi*f0*delay);
    
    % Coherent sum over all positions along the trajectory 
    S = S + Sn;
    % Inchoerent sum over all positions along the trajectory (choerent sum)
    A = A + abs(Sn);
end
close(wbar)


 Focus = (S./A)';
%Focus = S';
%Focus = A';


figure,
imagesc(x_ax,y_ax,abs(Focus)), axis xy , 
%title('Final Sensor Position') 
title('Focused image by TDBP')
xlabel('x [m]'), ylabel('y [m]')

hold on,
plot3(RX_pos(1,1),RX_pos(2,1),RX_pos(3,1),'ro'), plot3(RX_pos(1,end),RX_pos(2,end),RX_pos(3,end),'go'),
plot3(TX_pos(1,1),TX_pos(2,1),TX_pos(3,1),'kd'), plot3(cars(1,:),cars(2,:),cars(3,:),'rp'), plot3(humans(1,:),humans(2,:),humans(3,:),'yh')
legend('start swath','end swath','TX','Cars','Humans');