clear all
%%
Ptx_dbm = 10;
Ptx = 10^((Ptx_dbm-30)/10);             %watt
Gtx_db = 8;
Gtx = 10^(Gtx_db/10);

delta_elev = deg2rad(45);
delta_azim = deg2rad(60);
fc = 1.65e9;
lambda = physconst("LightSpeed")/fc;
samp_rate = 56e6;
L_txwave = 2^15;
PRI = L_txwave/samp_rate;
Tchirp = PRI *.9;

Aeff = Gtx * lambda^2 /4 /pi;

%%
scaling = delta_elev /2.7689;           % scale for -3dB beam
directivity_elev = @(x) (sin(x/scaling)/(x/scaling))^2; 


ang = deg2rad(-90:.1:90);
for i = 1:length(ang)
    a(i) = directivity_elev(ang(i));
end
x_ax = -90:.1:90;
figure,plot(x_ax,10*log10(a)),grid on,hold on, plot(x_ax,-3*ones(size(a))),title("Antenna Pattern Elevation")
%%
scaling = delta_azim /2.7689;           % scale for -3dB beam
directivity_azim  = @(x) (sin(x/scaling)/(x/scaling))^2; 

ang = deg2rad(-90:.1:90);
for i = 1:length(ang)
    a(i) = directivity_azim(ang(i));
end
x_ax = -90:.1:90;
figure,plot(x_ax,10*log10(a)),grid on,hold on, plot(x_ax,-3*ones(size(a))), title("Antenna Pattern Azimuth")

