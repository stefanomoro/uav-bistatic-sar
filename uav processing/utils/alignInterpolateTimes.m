function [OUT] = alignInterpolateTimes(radar,start_time,end_time,OUT)
%ALIGNINTERPOLATETIMES select correct time and interpolate
%   [OUT] = alignInterpolateTimes(radar,start_time,end_time,OUT)
time_s = hour(OUT.time)*3600 +  minute(OUT.time)*60 + second(OUT.time);
idxs = logical((time_s>=start_time) .* (time_s <=end_time));

OUT.time = OUT.time(idxs);
OUT.pos = OUT.pos(:,idxs);
OUT.speed = OUT.speed(idxs);


% Interpolate RX axis with Radar N_PRI axis
old_time_ax = OUT.time;
if(length(unique(old_time_ax)) < length(old_time_ax))
    old_time_ax = linspace(old_time_ax(1),old_time_ax(end),length(old_time_ax));
end

new_time_ax = linspace(OUT.time(1),OUT.time(end),radar.N_PRI);

OUT.time = interp1(old_time_ax, OUT.time, new_time_ax);

OUT.speed = interp1(old_time_ax, OUT.speed, new_time_ax);

pos_interp = zeros(3,length(new_time_ax));
for n = 1:3
    % interpolate and average to reduce step-like behaviour
    pos_interp(n,:) = movmean(interp1(old_time_ax, OUT.pos(n,:), new_time_ax),1e3); 
end
OUT.pos = pos_interp;
end

