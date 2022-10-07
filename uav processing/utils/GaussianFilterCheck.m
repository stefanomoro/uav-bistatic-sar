clc, close all
%% Gaussian function
% x = 0:dt:length(tx_wave)*dt-dt;
dt = 1 / 56e6;
L = 101;
x_fil = ((0:L-1)-floor(L/2))*dt;
sigma = 2 * dt;
%sigma = 1/chirp_bw;

%Peak at the centre
% d = (x(end)-x(1)) / 2 + x(1);
% g = 1/(sigma*sqrt(2*pi)) * exp(-0.5*((x-d)./sigma).^2); 
g = 1/(sigma*sqrt(2*pi)) * exp(-0.5*(x_fil/sigma).^2); 

g = g/ max(g);
figure, plot(x_fil,abs(g))
title("g")

%% Gauss wavenumber window
lambda = physconst("LightSpeed")/1.65e9;
ang = deg2rad(-90:.1:90);
k_rx = sin(ang).*(2*pi/lambda);
pho_az = 1;
Dk = 2*pi / pho_az;
sigma = Dk/2;
%Peak at the centre
% d = (x(end)-x(1)) / 2 + x(1);
% g = 1/(sigma*sqrt(2*pi)) * exp(-0.5*((x-d)./sigma).^2); 
g = 1/(sigma*sqrt(2*pi)) * exp(-0.5*(k_rx/sigma).^2); 
max(g)
g = g/ max(g);
figure, plot(rad2deg(ang),abs(g))
title("g")

%% sinc raw
x =( (0:2^15 - 1) - 2^14) * dt;

f = sinc(x/ dt * .9);%conv(tx_wave,conj(tx_wave),'same');
f = f / max(abs(f));
figure, plot(x,abs(f))
title("f")

F = fftshift(fft(f));
figure,plot(abs(F))
title("F")


%% mult Freq domain

G = fftshift(fft(g,length(F)));
figure,plot(abs(G))
title("G")
Y = G.*F;
figure,plot(abs(Y))
title("Y")
y = ifftshift(ifft(Y));
figure,plot(abs(y))
title("y")
figure,plot(f),title("f")

%% Filter
figure
plot(x,abs(f))


% w = g;
% w = w/sum(w);
% w = g;
y = conv2(f,g,'same');

hold on, plot(x_fil,abs(g))
hold on, plot(x,abs(y))

legend('f','g','y')