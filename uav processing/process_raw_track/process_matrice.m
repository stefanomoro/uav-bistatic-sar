%%
file_list = dir("..\..\mat_files\uav_test\20220826\drone track\raw\matrice");
file_list = file_list(3:end);
i = 1; flight_n = 1;
while  flight_n <= length(file_list)

ext = strsplit(file_list(i).name,".");
ext = ext{2};
while (not(strcmp(ext,"csv")))
    i = i+1;
    ext = strsplit(file_list(i).name,".");ext = ext{2};
end
fname = strcat(file_list(i).folder, "\",file_list(i).name);

gps_raw = readtable(fname);
%%
gps_time = datetime(gps_raw.datetime_utc_(1), ...
'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z','Format','yyyy-MM-dd HH:mm:ss.SSS') + hours(2);


gps_time = gps_time + milliseconds(gps_raw.time_millisecond_) - milliseconds(gps_raw.time_millisecond_(1));

gps_matrice = table();
gps_matrice.datetime = gps_time;
gps_matrice.lat = gps_raw.latitude;
gps_matrice.lon = gps_raw.longitude;
gps_matrice.alt = gps_raw.height_above_takeoff_meters_;
gps_matrice.speed = gps_raw.speed_m_s_;
gps_matrice.yaw = gps_raw.compass_heading_degrees_;
gps_matrice.pitch = gps_raw.pitch_degrees_;
gps_matrice.roll = gps_raw.roll_degrees_;

%% to UTM

utmZ = utmzone(gps_matrice.lat,gps_matrice.lon)
[ellipsoid,estr] = utmgeoid(utmZ)

utmstruct = defaultm('utm');
utmstruct.zone = utmZ;
utmstruct.geoid = ellipsoid;
utmstruct = defaultm(utmstruct);

[gps_matrice.utm_x,gps_matrice.utm_y] = projfwd(utmstruct,gps_matrice.lat,gps_matrice.lon)


figure, plot3(gps_matrice.utm_x,gps_matrice.utm_y,gps_matrice.alt),grid on
xlabel("X"),ylabel("Y"), axis("xy")
%%
save(strcat("track_matrice_",num2str(flight_n)),"gps_matrice")
flight_n = flight_n + 1;
i = i +1;
end
return
%% targets
targets = [];
targets.lat(1) = 45.7875685; targets.lon(1)= 9.3201195;
targets.lat(2) = 45.7877989; targets.lon(2) = 9.3204629;
targets.lat(3) = 45.7877529; targets.lon(3) = 9.3202852;

utmZ = utmzone(targets.lat(1),targets.lon(1))
[ellipsoid,estr] = utmgeoid(utmZ)

utmstruct = defaultm('utm');
utmstruct.zone = utmZ;
utmstruct.geoid = ellipsoid;
utmstruct = defaultm(utmstruct);
for i = 1:length(targets.lat)
[targets.utm_x(i),targets.utm_y(i)] = projfwd(utmstruct,targets.lat(i),targets.lon(i))
end
figure,plot(targets.utm_x,targets.utm_y,'ro')
save("targets","targets")