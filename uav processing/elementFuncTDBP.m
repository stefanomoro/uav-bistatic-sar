function [s0] = elementFuncTDBP(x0,y0,z0,TX_pos_x,TX_pos_y,TX_pos_z,RX_pos_x,RX_pos_y,RX_pos_z,lambda,Dk,RC,t,f0,k_rx_0)
%ELEMENTFUNCTDBP Summary of this function goes here
%   Detailed explanation goes here

R_tx = sqrt((TX_pos_x - x0).^2 + (TX_pos_y - y0).^2  + ...
    (TX_pos_z - z0).^2);                                           %Range distances from the tx antenna [m]

R_rx = sqrt((RX_pos_x-x0).^2 + (RX_pos_y - y0).^2  + ...
    (RX_pos_z - z0).^2);                                           %Range distances from the rx antenna [m]
distance = R_tx+R_rx;                                               %Total Tx-target-Rx distance [m]
delay = distance./physconst('LightSpeed');

%Compute target wave number
R = sqrt((RX_pos_x - x0).^2 + (RX_pos_y - y0).^2);
psi = asin((y0-RX_pos_y)./R);
k_rx = sin(psi).*(2*pi/lambda);

%Weight function
%         Wn = rectpuls((k_rx - k_rx_0)./psi_proc);
sigma = Dk/2;
Wn = gaussActivFunc(k_rx - k_rx_0,sigma);

cut = x0 < RX_pos_x ;                                       %Cut the back-lobe
Wn = Wn * cut;
% Backprojection of data from a single Radar position 
s0 = Wn.*interp1(t,RC(:),delay).*exp(+1i*2*pi*f0*delay);
end

