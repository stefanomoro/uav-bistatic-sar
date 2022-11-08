clear variables
%% Constants
c = physconst("LightSpeed");

Ptx_dbm = 10;
Ptx = 10^((Ptx_dbm-30)/10);             %watt
Gtx_db = 8;
Gtx = 10^(Gtx_db/10);
speed = 4;                   %m/s

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

% SYSTEM Resolution
pho_az = .5;
real_samp = pho_az/4;
az_samp = lambda/4;
grid_res = pho_az/2;

% dx = speed * PRI;
pho_R = c/B;
NF = 8;                                %dB
F_noise = 10^(NF/10);

disp(strcat("Max speed is ",num2str(lambda/4/PRI),"m/s or ",num2str(lambda/4/PRI * 3.6),"km/h"))
%% Antenna Elevation Pattern
scaling = delta_elev /2.7689;           % scale for -3dB beam
directivity_elev = @(x) (sinc(x/scaling)).^2; 
%% Antenna Azimuth Patterm
scaling = delta_azim /2.7689;           % scale for -3dB beam
directivity_azim  = @(x) (sinc(x/scaling)).^2; 


%% Activation function
Dk = 2*pi / pho_az;

gaussActiv = @(x) exp(-0.5*((x)./Dk).^2);

%% Simulate

x_ap = -10:real_samp:10;           % ground position of distributed target
x = -10:grid_res:10;          % focalization grid
chi = -30:az_samp:30;            % drone sensor pos

N_sim = 200;
foc_gain = zeros(N_sim,1);
pow_orig = zeros(N_sim,1);
pow_foc = zeros(N_sim,1);
wbar = waitbar(0,"Simulation");

for sim_i = 1:N_sim

% compute s(x_ap)

s = complex(randn(size(x_ap)),randn(size(x_ap)) )/sqrt(2);      % distributed target
% s = zeros(size(x_ap)); s(floor(length(s)/2)) = 1;               % single target

y = 50;                     % ground distance of s
z = 30;                     % height of sensor

% compute d(chi)

x_ap_men_chi = x_ap(:) - chi(:).';
R = sqrt(  (x_ap_men_chi).^2 + z^2 + y^2);
psi = atan((x_ap_men_chi) /y);
elev = atan(-z/y) - ant_elev_pointing_rx;
f_ant = directivity_azim(psi).* directivity_elev(elev);
        
d = s(:).* exp(-1i*2*pi / lambda * 2 *R).*f_ant.^2;
d = sum(d,1);               %projected data

% compute F(x)

chi_men_x = chi(:) - x(:).';
R = sqrt(  (chi_men_x).^2 + z^2 + y^2);
psi = atan(chi_men_x/y);
elev = atan(z/y);
k_rx = cos(elev).* sin(psi).*(4*pi/lambda); % 4 pi for monostatic
F = d(:) .* exp(1i*2*pi / lambda * 2*R).*gaussActiv(k_rx);
F = sum(F,1);               % retroprojected data    

N_tau = length(x) * length(chi);

pow_orig(sim_i) = mean(abs(s)).^2;
pow_foc(sim_i) = mean(abs(F).^2);
foc_gain(sim_i) = mean(abs(F).^2)/mean(abs(s).^2);

waitbar(sim_i/N_sim,wbar);
    
end
close(wbar)
%% Compare output

disp(strcat("Mean s pow ", num2str(mean(pow_orig)))   )
disp(strcat("Mean F pow ", num2str(mean(pow_foc))))
disp(strcat("Mean Foc gain ", num2str(10*log10(mean(foc_gain))), " dB" ) )

figure,plot(10*log10(foc_gain)),hold on,yline(10*log10(mean(foc_gain)),'r--'), title("Focalization gain")


figure,
subplot(2,1,1)
stem(x_ap,abs(s).^2),title("Original target")
subplot(2,1,2)
plot(x, abs(F).^2),title("Focused image")
%% W function
tic
W = zeros(length(x_ap),length(x));
for chi_i = 1:length(chi)
        
    R1 = sqrt((chi(chi_i) -x_ap).^2 + y^2 + z^2);
    R2 = sqrt((chi(chi_i) -x).^2 + y^2 + z^2);
    DR = R1(:) - R2(:).';
    psi = atan((chi(chi_i)-x)/y);
    elev = atan(z/y);
    elev_ant = atan(-z/y) - ant_elev_pointing_rx;
    f_ant = directivity_azim(psi) * directivity_elev(elev_ant);
    
    k_rx = cos(elev) * sin(psi).*(4*pi/lambda); % 4 pi for monostatic
    W = W + exp(-1i*4*pi/lambda * DR).* f_ant.^2 .* gaussActiv(k_rx);
end

toc
figure,imagesc(x_ap,x,abs(W)),axis xy,xlabel("x'"),ylabel("x"),colormap jet

%%

G_foc = mean(sum(abs(W).^2,1));

disp(strcat("math Foc gain ", num2str(10*log10(G_foc)), " dB" ) )
disp(strcat("Mean Foc gain ", num2str(10*log10(mean(foc_gain))), " dB" ) )
