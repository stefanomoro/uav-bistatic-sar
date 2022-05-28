
%% ========================================  Load data drone RX

time = drone_data.time;
speed = drone_data.speed;

RX_pos = [drone_data.x(:) drone_data.y(:)  drone_data.alt(:)].';
RX_speed = speed;
RX_timeDrone = minute(time)*60 + second(time);

%% Find the Swath samples before cut

% Find the time radar start/end
RX_timeRadar_total_end = minute(radar_last_mod_date)*60 + second(radar_last_mod_date);
RX_timeRadar_total_start = RX_timeRadar_total_end - N_PRI_tot*PRI;

% Find the idxs time drone start/end
start_idxs = find(RX_timeDrone>RX_timeRadar_total_start);
end_idxs = find(RX_timeDrone>RX_timeRadar_total_end);

RX_timeDrone_total_start_idx_ = start_idxs(1);
RX_timeDrone_total_end_idx = end_idxs(1)-1;

%% Find the Swath samples after cut

% Find the time radar end of the cutted radar data
RX_timeRadar_start = RX_timeRadar_total_start;
RX_timeRadar_end = RX_timeRadar_total_start + N_PRI*PRI;              

% Find the idxs time drone start/end
start_idxs = find(RX_timeDrone>RX_timeRadar_start);
end_idxs = find(RX_timeDrone>RX_timeRadar_end);

RX_timeDrone_start_idx = start_idxs(1);
RX_timeDrone_end_idx = end_idxs(1)-1;



%% ========================================== Load data TX
TX_pos = [0;0;2].*ones(3,N_PRI);
TX_speed = zeros(N_PRI,1);

%% ========================================== Load data Targets
cars = targets.cars;
humans = targets.human;

cars = [cars.x, cars.y, zeros(size(cars.x))]';
humans = [humans.x, humans.y, zeros(size(humans.x))]';
%% Plot
drone_start_pos = RX_pos(:,RX_timeDrone_start_idx);
drone_end_pos = RX_pos(:,RX_timeDrone_end_idx);

figure, plot3(RX_pos(1,:),RX_pos(2,:),RX_pos(3,:)), title('Flight path'), xlabel('x'), ylabel('y'), zlabel('z')

hold on,
plot3(drone_start_pos(1),drone_start_pos(2),drone_start_pos(3),'ro'), 
plot3(drone_end_pos(1),drone_end_pos(2),drone_end_pos(3),'go'),
plot3(TX_pos(1,1),TX_pos(2,1),TX_pos(3,1),'ko'), plot3(cars(1,:),cars(2,:),cars(3,:),'p'), plot3(humans(1,:),humans(2,:),humans(3,:),'h')
legend('RX drone path','start swath','end swath','TX','Cars','Humans');grid on


%% Cut the drone time axis 
RX_timeDrone = RX_timeDrone(RX_timeDrone_start_idx:RX_timeDrone_end_idx);
RX_pos = RX_pos(:,RX_timeDrone_start_idx:RX_timeDrone_end_idx);
RX_speed = speed(RX_timeDrone_start_idx:RX_timeDrone_end_idx);

% %%
% figure, plot3(RX_pos(1,:),RX_pos(2,:),RX_pos(3,:)), title('Flight path'), xlabel('x'), ylabel('y'), zlabel('z')
% %%
% figure, plot(RX_speed)

%% Interpolate RX axis with Radar N_PRI axis
old_time_ax = RX_timeDrone;
new_time_ax = linspace(RX_timeDrone(1),RX_timeDrone(end),N_PRI);

RX_timeDrone = interp1(old_time_ax, RX_timeDrone, new_time_ax);

RX_speed = interp1(old_time_ax, RX_speed, new_time_ax);

RX_pos_interp = [];
for n = 1:3
    % interpolate and average to reduce step-like behaviour
    RX_pos_interp(n,:) = movmean(interp1(old_time_ax, RX_pos(n,:), new_time_ax),1e3); 
end

RX_pos = RX_pos_interp;

clear old_time_ax new_time_ax RX_pos_interp
%% Plot
figure, plot3(RX_pos(1,:),RX_pos(2,:),RX_pos(3,:)), title('Flight path'), xlabel('x'), ylabel('y'), zlabel('z'),grid on


%% ========================================== Compute distance-TX-RX
distance_tx_rx = zeros(N_PRI,1);
for n = 1:N_PRI
distance_tx_rx(n) = sqrt( sum((TX_pos(:,n) - RX_pos(:,n)).^2)) ;

end
%% Plot
figure, plot(distance_tx_rx), xlabel('N_{PRI}'), ylabel('Distance [m]'), title('distance TX RX')

%% ========================================== Speed distance-TX-RX
speed_distance_tx_rx = [0;diff(distance_tx_rx(:))];

speed_distance_tx_rx = speed_distance_tx_rx ./ PRI;

%% Plot
figure, plot(speed_distance_tx_rx),hold on, plot(RX_speed) 
legend("delta Distance / PRI", "Speed from log"),xlabel('N_{PRI}'), ylabel('Velocity [m/s]'), title('speed distance TX RX')

%% ==========================================  distance Targets
%% cars
distance_cars = zeros(size(cars,2),N_PRI);
for n = 1:N_PRI
    for car_idx = 1:size(cars,2)
        distance_cars(car_idx,n) = sqrt(sum((TX_pos(:,n) - cars(:,car_idx)).^2)) + sqrt(sum( (RX_pos(:,n) - cars(:,car_idx)).^2));
    end
end

%% humans
distance_humans = zeros(size(humans,1),N_PRI);
for n = 1:N_PRI
    for hum_idx = 1:size(humans,2)
        distance_humans(hum_idx,n) = sqrt(sum((TX_pos(:,n) - humans(:,hum_idx)).^2)) + sqrt(sum((RX_pos(:,n) - humans(:,hum_idx)).^2));

    end
end
