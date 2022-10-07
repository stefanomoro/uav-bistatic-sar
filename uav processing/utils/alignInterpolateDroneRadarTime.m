function [RX,TX,radar] = alignInterpolateDroneRadarTime(const,drone,radar)
%ALIGNDRONERADARTIME align times from radar data and drone data
%   [RX,TX,radar] = alignInterpolateDroneRadarTime(const,drone,radar)
if(strcmp(const.setup_mode,"old"))
    RX.pos = [drone.x(:) drone.y(:)  drone.alt(:)].';
    RX.speed = drone.speed;
    RX.time = drone.time;
else
    RX.pos = [drone.rx.utm_x(:) drone.rx.utm_y(:) drone.rx.alt(:)].';
    RX.speed = drone.rx.speed;
    RX.time = drone.rx.datetime;
end

if(isfield(const,"radar_start_timestamp"))
    radar.time_raw.start = hour(const.radar_start_timestamp)*3600 + ...
        minute(const.radar_start_timestamp)*60 + second(const.radar_start_timestamp);
    radar.time_raw.end = radar.time_raw.start + radar.N_PRI_raw*const.PRI;
else
    % Find the time radar start/end
    radar.time_raw.end = hour(radar.last_mod_date)*3600 + ...
        minute(radar.last_mod_date)*60 + second(radar.last_mod_date);
    radar.time_raw.start = radar.time_raw.end - radar.N_PRI_raw*const.PRI;
end
% evaluate window for taking only processed samples
start_time = (radar.idx_wind_col(1)-1)*const.PRI +radar.time_raw.start;
end_time = (radar.idx_wind_col(end)-1)*const.PRI + radar.time_raw.start;

% find the related time of drone with radar and interpolate
[RX] = alignInterpolateTimes(radar,start_time,end_time,RX);

if(strcmp(const.setup_mode,"old"))
    TX.pos = [0;0;2].*ones(3,radar.N_PRI);
    TX.speed = zeros(radar.N_PRI,1);
else
    if(strcmp(const.setup_mode,"bistatic"))
        TX.pos = [drone.tx.utm_x(:) drone.tx.utm_y(:) drone.tx.alt(:)].';
        TX.speed = zeros(radar.N_PRI,1);
        TX.time = drone.tx.datetime;
        
        [TX] = alignInterpolateTimes(radar,start_time,end_time,TX);
        % average the position of Tarot, cause too noisy
        mean_x = mean(TX.pos(1,:)) .* ones(1,size(TX.pos,2));
        mean_y = mean(TX.pos(2,:)) .* ones(1,size(TX.pos,2));
        mean_z = mean(TX.pos(3,:)) .* ones(1,size(TX.pos,2));
        TX.pos = [mean_x;mean_y;mean_z];
    end
end

end

