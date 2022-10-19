%% Constants
c = physconst("LightSpeed");

Ptx_dbm = 10;
Ptx = 10^((Ptx_dbm-30)/10);             %watt
Gtx_db = 8;
Gtx = 10^(Gtx_db/10);

delta_elev = deg2rad(45);
delta_azim = deg2rad(60);
ant_elev_pointing_rx = deg2rad(-28);

fc = 2.4e9;
lambda = c/fc;
samp_rate = 56e6;
L_txwave = 2^15;
PRI = L_txwave/samp_rate;
Tchirp = PRI *.9;
B = samp_rate *.9;

A_eff = Gtx * lambda^2 /4 /pi;


pho_az = 1;
pho_R = c/B;
NF = 8;                                %dB
F_noise = 10^(NF/10);
%% Antenna Elevation Pattern
scaling = delta_elev /2.7689;           % scale for -3dB beam
directivity_elev = @(x) (sinc(x/scaling))^2; 
%% Antenna Azimuth Patterm
scaling = delta_azim /2.7689;           % scale for -3dB beam
directivity_azim  = @(x) (sinc(x/scaling))^2; 
%% Activation function
pho_az = .1;             % azimuth resolution
Dk = 2*pi / pho_az;
sigma = Dk / 2;
gaussActiv = @(x) exp(-0.5*((x)./sigma).^2);

%% 
x_ap = -20:pho_az/2:20;           % ground position of distributed target
x = x_ap;
s = complex(randn(size(x_ap)),randn(size(x_ap)) )/sqrt(2);      % distributed target
% s = zeros(size(x_ap)); s(floor(length(s)/2)) = 1;               % single target
y = 50;                     % ground distance of s

z = 30;

chi = -20:pho_az:20;            % drone sensor pos
d = zeros(size(chi));       % projected data
F = zeros(size(x));         % retroprojected data
illumination = zeros(size(x));
for i = 1:length(chi)
    for j = 1:length(x_ap)
        R = sqrt(  (x_ap(j) - chi(i))^2 + z^2 + y^2);
        psi = atan((x_ap(j)-chi(i)) /y);
        elev = atan(-z/y) - ant_elev_pointing_rx;
        f_ant = directivity_azim(psi) * directivity_elev(elev);
        d(i) = d(i) + s(j) * exp(-1i*4*pi / lambda * R) *f_ant^2;
        illumination(j) = illumination(j) +f_ant;
    end
end

% compute F(x)
for i = 1:length(x)
    for j = 1:length(chi)
        R = sqrt(  (chi(j) - x(i))^2 + z^2 + y^2);
        psi = atan(chi(j)-x(i)/y);
        k_rx = sin(psi).*(2*pi/lambda);
        F(i) = F(i) + d(j) * exp(1i*4*pi / lambda * R) *gaussActiv(k_rx);
    end
end
N_tau = length(x) * length(chi);

%% Compare output
disp(strcat("Mean s pow ", num2str(mean(abs(s).^2))))
disp(strcat("Mean F pow ", num2str(mean(abs(F).^2))))

foc_gain = mean(abs(F).^2);
ratio = max(abs(F))/max(abs(s));

disp(strcat("Ratio between max values ",num2str(ratio)))

figure,plot(x_ap,abs(s).^2),hold on, plot(x, abs(F).^2)
% figure,plot(x,illumination),title("Antenna Illumination")