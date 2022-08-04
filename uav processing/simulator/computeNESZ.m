clear all
%% Constants
c = physconst("LightSpeed");

Ptx_dbm = 10;
Ptx = 10^((Ptx_dbm-30)/10);             %watt
Gtx_db = 8;
Gtx = 10^(Gtx_db/10);

delta_elev = deg2rad(45);
delta_azim = deg2rad(60);
ant_elev_pointing = deg2rad(-28);

fc = 1.65e9;
lambda = c/fc;
samp_rate = 56e6;
L_txwave = 2^15;
PRI = L_txwave/samp_rate;
Tchirp = PRI *.9;
B = samp_rate *.9;

A_eff = Gtx * lambda^2 /4 /pi;

drone_height = 15:15:90;                      %m

drone_speed = 5;                        %m/s
track_resolution = PRI * drone_speed;
flight_time = 20;                       %s
tau_ax = 0:PRI:flight_time;

pho_az = 1;
pho_R = c/2/B;
sigma_0_db = -6.4;
sigma_0 = 10^(sigma_0_db/10);
NF = 8;                                %dB
F_noise = 10^(NF/10);

for height_idx = 1:length(drone_height)

                     
max_pointing_2d_distance = - drone_height(height_idx) / tan(ant_elev_pointing);

RX_pos = [zeros(size(tau_ax));linspace(0,drone_speed * PRI * length(tau_ax),length(tau_ax)); drone_height(height_idx) * ones(size(tau_ax))];

%% antenna pattern elev
scaling = delta_elev /2.7689;           % scale for -3dB beam
directivity_elev = @(x) (sin(x/scaling)/(x/scaling))^2; 

ang = deg2rad(-90:.1:90);
for i = 1:length(ang)
    a(i) = directivity_elev(ang(i));
end
x_ax = -90:.1:90;
% figure,plot(x_ax,10*log10(a)),grid on,hold on, plot(x_ax,-3*ones(size(a))),title("Antenna Pattern Elevation")
%% antenna pattern azim
scaling = delta_azim /2.7689;           % scale for -3dB beam
directivity_azim  = @(x) (sin(x/scaling)/(x/scaling))^2; 

ang = deg2rad(-90:.1:90);
for i = 1:length(ang)
    a(i) = directivity_azim(ang(i));
end
x_ax = -90:.1:90;
% figure,plot(x_ax,10*log10(a)),grid on,hold on, plot(x_ax,-3*ones(size(a))), title("Antenna Pattern Azimuth")


%% CYCLE


tgt_2d_distance = 30:200;
SNR_db = zeros(size(tgt_2d_distance));
P_noise_db = zeros(size(tgt_2d_distance));
P_sig_db = zeros(size(tgt_2d_distance));
nesz_db = zeros(size(tgt_2d_distance));
synt_aperture_vec = zeros(size(tgt_2d_distance));

for dist_idx = 1:length(tgt_2d_distance)
tgt_pos = [tgt_2d_distance(dist_idx);50;0];


%% Power received

tgt_distance = sqrt(sum((RX_pos - tgt_pos).^2));
ant_azim = atan((RX_pos(2,:) - tgt_pos(2)) ./ (RX_pos(1,:) - tgt_pos(1)));
f_ant = zeros(size(ant_azim));
for i = 1:length(ant_azim)
    f_ant(i) = directivity_azim(ant_azim(i));
end
ant_elev = atan((RX_pos(3,1) - tgt_pos(3)) ./ (RX_pos(1,1) - tgt_pos(1))) ;
f_ant = f_ant * directivity_elev(ant_elev- ant_elev_pointing);

%monostatic case
A_ill = pho_az * pho_R / sin(pi/2 + ant_elev);
RCS = sigma_0 * A_ill;
% max_delay = 200/c;
% t = 0:1/samp_rate/4:max_delay;

% s_rc = zeros(length(t),length(tau_ax));
s_rc = zeros(length(tau_ax),1);
Prx = zeros(length(ant_azim));
for i= 1:length(tau_ax)
    Prx(i) = Ptx * Gtx * f_ant(i)^2 * RCS *A_eff / (4*pi*tgt_distance(i)^2)^2;
    delay = tgt_distance(i)/c;
%     s_rc(:,i) = sinc((t-delay).*B) .*sqrt(Prx(i)) .* Tchirp .*exp(-2*1i*pi*fc*delay);
    s_rc(i) = sqrt(Prx(i)) .* Tchirp;
end
%%
x_dist = tgt_pos(1) - RX_pos(1,1);
psi_proc = lambda / 2 / pho_az;
synt_apert = 2 * tan(psi_proc/2) * x_dist;
N_synt_apert = floor(synt_apert/track_resolution);
pos_0_idx = find(RX_pos(2,:) > tgt_pos(2),1);
selected_pos_idx =  pos_0_idx - floor(N_synt_apert/2): pos_0_idx+ floor(N_synt_apert/2);

N_tau = sum(sqrt(Prx(selected_pos_idx))) / sqrt(Prx(pos_0_idx));

s_foc_tgt = abs(s_rc(pos_0_idx)) * N_tau;

%% SNR
P_noise = physconst("Boltzman") * 290 * F_noise * N_tau * Tchirp;
P_sig = s_foc_tgt^2;
N0 = physconst("Boltzman")*290 * F_noise ;

nesz = N0 * (4*pi*tgt_distance(pos_0_idx)^2)^2 / (N_tau * Tchirp * Ptx * Gtx * f_ant(pos_0_idx)^2 * A_ill * A_eff) ;
nesz_db(dist_idx) = 10*log10(nesz);

SNR = P_sig / P_noise;
SNR_db(dist_idx) = 10*log10(SNR);
P_noise_db(dist_idx) = 10*log10(P_noise);
P_sig_db(dist_idx) = 10*log10(P_sig);
synt_aperture_vec(dist_idx) = synt_apert;
end
%%
% figure,plot(tgt_2d_distance,SNR_db,tgt_2d_distance,P_noise_db,tgt_2d_distance,P_sig_db)
% legend("SNR","P_{noise}","P_{signal}"),grid on

figure,plot(tgt_2d_distance,nesz_db),ylabel("NESZ [dB]")
grid on
yyaxis right
plot(tgt_2d_distance,synt_aperture_vec)
xline(max_pointing_2d_distance,'--r')
text([max_pointing_2d_distance + 1 ],[15],"Max antenna gain")
ylabel("Synthetic aperture [m]"),title(strcat("Drone heigth ", num2str(drone_height(height_idx))))
end
