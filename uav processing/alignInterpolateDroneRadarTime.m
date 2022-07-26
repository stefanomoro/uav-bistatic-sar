function [RX,TX,radar] = alignInterpolateDroneRadarTime(const,drone,radar)
%ALIGNDRONERADARTIME align times from radar data and drone data
%   [RX,TX,radar] = alignInterpolateDroneRadarTime(const,drone,radar)
RX.pos = [drone.x(:) drone.y(:)  drone.alt(:)].';
RX.speed = drone.speed;
RX.time = drone.time;
time_s = minute(drone.time)*60 + second(drone.time);


% Find the time radar start/end
radar.time_raw.end = minute(radar.last_mod_date)*60 + second(radar.last_mod_date);
radar.time_raw.start = radar.time_raw.end - radar.N_PRI_raw*const.PRI;

% evaluate window for taking only processed samples
start_time = (radar.idx_wind_col(1)-1)*const.PRI +radar.time_raw.start;
end_time = (radar.idx_wind_col(end)-1)*const.PRI + radar.time_raw.start;

% find the related time of drone with radar
idxs = and(time_s>start_time, time_s <end_time);

RX.time = drone.time(idxs);
RX.pos = RX.pos(:,idxs);
RX.speed = RX.speed(idxs);


% Interpolate RX axis with Radar N_PRI axis
old_time_ax = RX.time;
new_time_ax = linspace(RX.time(1),RX.time(end),radar.N_PRI);

RX.time = interp1(old_time_ax, RX.time, new_time_ax);

RX.speed = interp1(old_time_ax, RX.speed, new_time_ax);

RX_pos_interp = zeros(3,length(new_time_ax));
for n = 1:3
    % interpolate and average to reduce step-like behaviour
    RX_pos_interp(n,:) = movmean(interp1(old_time_ax, RX.pos(n,:), new_time_ax),1e3); 
end

RX.pos = RX_pos_interp;

TX.pos = [0;0;2].*ones(3,radar.N_PRI);
TX.speed = zeros(radar.N_PRI,1);

end

