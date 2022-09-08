%% LIST FILES
file_list = dir("..\..\mat_files\uav_test\20220826\drone track\raw\tarot\IMU");
file_list = file_list(3:end);
i = 1;
ext = strsplit(file_list(i).name,".");
ext = ext{2};
while (not(strcmp(ext,"txt")))
    i = i+1;
    ext = strsplit(file_list(i).name,".");ext = ext{2};
end
fname = strcat(file_list(i).folder, "\",file_list(i).name);
imu_tarot = readtable(fname);
%%
file_list = dir("..\..\mat_files\uav_test\20220826\drone track\raw\tarot\gps");
file_list = file_list(3:end);
i = 1;
ext = strsplit(file_list(i).name,".");
ext = ext{2};
while (not(strcmp(ext,"txt")))
    i = i+1;
    ext = strsplit(file_list(i).name,".");ext = ext{2};
end
fname = strcat(file_list(i).folder, "\",file_list(i).name);
gps_tarot = readtable(fname);

%%
% read 2 files, IMU ->imu_tarot and GPS ->gps_tarot

gps_time = datetime(gps_tarot.datetime, ...
'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z','Format','yyyy-MM-dd HH:mm:ss.SSS') + hours(2);
gps_tarot.datetime = gps_time;

imu_time = datetime(imu_tarot.datetime, ...
'InputFormat','HH:mm:ss.SSS','Format','yyyy-MM-dd HH:mm:ss.SSS');
imu_time.Day = gps_time.Day(1);
imu_time.Month = gps_time.Month(1);
imu_time.Year = gps_time.Year(1);
imu_tarot.datetime = imu_time;

%% altitude

idxs = and(imu_tarot.datetime > datetime("2022-08-26 11:45:00"),imu_tarot.datetime < datetime("2022-08-26 11:50:00"));
press = movmean(imu_tarot.press,31);
P_0 = mean(press(idxs));

C = - 0.02896 * 9.8064 /(8.3143*288.15);
%P = P_0 * exp(C . h);
h = log(press/P_0)/C;
imu_tarot.altitude = h;
%% gps altitude
idxs = gps_tarot.alt > 150;
figure,plot(gps_time(idxs),gps_tarot.alt(idxs))

%% CUT only experiment
start_date = datetime("2022-08-26 11:53:00");
end_date = datetime("2022-08-26 11:59:00");
idxs = and(imu_tarot.datetime > start_date,imu_tarot.datetime < end_date);

figure,plot(imu_time(idxs),imu_tarot.altitude(idxs)),hold on,
tot_acc = imu_tarot.x_acc + imu_tarot.y_acc + imu_tarot.z_acc;
plot(imu_time(idxs),tot_acc(idxs))   

imu_out = imu_tarot(idxs,:);

idxs = and(gps_tarot.datetime > start_date,gps_tarot.datetime < end_date);
gps_out = gps_tarot(idxs,:);
%% FIX altitude
alt = imu_out.altitude;
mean_ground = mean([alt(1) alt(end)]);

alt = alt - max(alt(1),alt(end));
alt(alt<0) = 0;
alt = alt + mean_ground;
alt(alt== mean_ground) = 0;
imu_out.altitude = alt;
figure,plot(imu_out.datetime,alt)

%% SAVE output
dt = (gps_out.datetime(end) - gps_out.datetime(1))/length(gps_out.datetime);
dt = milliseconds(dt)*1e-3;
x = hour(gps_out.datetime(1))*3600 +  60*minute(gps_out.datetime(1)) +second(gps_out.datetime(1)) + [0:dt:dt*(length(gps_out.datetime)-1)  ];

dt = (imu_out.datetime(end) - imu_out.datetime(1))/length(imu_out.datetime);
dt = milliseconds(dt)*1e-3;
xq = hour(imu_out.datetime(1))*3600 + 60*minute(imu_out.datetime(1)) + second(imu_out.datetime(1)) + [0: dt:dt*(length(imu_out.datetime)-1)];


lat = interp1(x,gps_out.lat,xq,"linear","extrap");
lon = interp1(x,gps_out.lon,xq,"linear","extrap");
%%
tarot = imu_out;
tarot.lat = lat(:);
tarot.lon = lon(:);
%%
save("track_tarot","tarot")
% figure,plot(imu_out.datetime,lat),figure,plot(gps_out.datetime,gps_out.lat,'ko')


%% Static Transmitter
start_date = datetime("2022-08-26 12:11:00");
end_date = datetime("2022-08-26 12:28:00");

idxs = and(gps_tarot.datetime > start_date,gps_tarot.datetime < end_date);
stable_gps = gps_tarot(idxs,:);

idxs = and(imu_tarot.datetime > start_date,imu_tarot.datetime < end_date);
stable_imu = imu_tarot(idxs,:);
%% to x/y

R =  physconst("earth");%6.3781* 1e6;
lat  = deg2rad(stable_gps.lat);
lon = deg2rad(stable_gps.lon);

lat_0 = deg2rad(45.7873713);
lon_0 = deg2rad(9.3203187);
x = R * (lon - lon_0)* cos(lat_0);
y = R * (lat-lat_0);

figure,plot(x,y), grid on
%% To UTM

utmZ = utmzone(stable_gps.lat,stable_gps.lon)
[ellipsoid,estr] = utmgeoid(utmZ)

utmstruct = defaultm('utm');
utmstruct.zone = utmZ;
utmstruct.geoid = ellipsoid;
utmstruct = defaultm(utmstruct);

[stable_gps.utm_x,stable_gps.utm_y] = projfwd(utmstruct,stable_gps.lat,stable_gps.lon)
stable_gps.alt = 4*ones(length(stable_gps.alt),1) %fixed height of the pole 
save("stable_gps","stable_gps")



