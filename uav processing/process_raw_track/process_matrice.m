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


gps_time = gps_time + milliseconds(gps_raw.time_millisecond_);

gps_matrice = table();
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