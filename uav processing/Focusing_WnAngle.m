%% =========================================================================== FOCUSING

%% FOCUSING BY TIME DOMAIN BACK PROJECTION (TDBP)
% In this case, we choose to focus on the plane x,y  with z = c
%% Set XY grid where to focus
x_min = -20;
x_max = 100;
pho_az = 1;
dx = pho_az*0.4;
x_ax =  x_min:dx:x_max;

y_min = -80;
y_max = 20;
dy = dx;
y_ax =  y_min:dy:y_max;

delta_psi_proc = lambda/pho_az;

[X,Y] = ndgrid(x_ax,y_ax);

z0 = 0;
%% Focusing (Wn form angle)

wbar = waitbar(0,'Backprojecting');
Ny = length(y_ax);
Nx = length(x_ax);
S = zeros(Nx,Ny);
A = zeros(Nx,Ny);
t = R_ax./c;

%Compute path angle
psi_path = atan( (RX_pos(2,end)-RX_pos(2,1)) / (RX_pos(1,end)-RX_pos(1,1)) );

%Compute pointing angle 
psi_point = psi_path+deg2rad(90);

%Focalization angle
% psi_foc = deg2rad(45); % = psi_point is better
psi_foc = psi_point;
clear Sn    
for n = 1:N_PRI
    waitbar(n/N_PRI,wbar)
   
    % Distance 
    R_tx = sqrt((TX_pos(1,n)-X).^2 + (TX_pos(2,n)-Y).^2  + (TX_pos(3,n)-z0).^2); %  Range distances from the tx antenna [m]
    R_rx = sqrt((RX_pos(1,n)-X).^2 + (RX_pos(2,n)-Y).^2  + (RX_pos(3,n)-z0).^2); %  Range distances from the rx antenna [m]
    distance = R_tx+R_rx; %Total Tx-target-Rx distance [m]
    delay = distance./c;    %Delay
    
    %Compute target angle
    R = sqrt((RX_pos(1,n)-X).^2 + (RX_pos(2,n)-Y).^2);
    psi = asin((Y-RX_pos(2,n))./R);

    %Weight function
    Wn = rectpuls((psi - psi_foc)./delta_psi_proc);
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

    plot(RX_pos(1,n),RX_pos(2,n),'dr')
   

    %% Path and pointing vector plot
%     x_vect = 0:dx:x_ax(end)-RX_pos(1,n);
%     pointing = RX_pos(2,n) + (x_vect)*tan(psi_path);
%     pointing = [ones(length(x_ax)-length(pointing),1)', pointing];
%     hold on, plot(x_ax,pointing,'y')
%     
%     pointing = RX_pos(2,n) + (x_vect)*tan(psi_point);
%     pointing = [ones(length(x_ax)-length(pointing),1)', pointing];
%     hold on, plot(x_ax,pointing,'r')
%     
   

%%
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