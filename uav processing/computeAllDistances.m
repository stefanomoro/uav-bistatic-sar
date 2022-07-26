function [scenario] = computeAllDistances(const,radar,RX,TX,targets)
%computeAllDistances compute all distances for TX RX and targets
%   [scenario] = computeAllDistances(const,radar,RX,TX)

% TX-RX distance
distance_tx_rx = zeros(radar.N_PRI,1);
for n = 1:radar.N_PRI
    distance_tx_rx(n) = sqrt( sum((TX.pos(:,n) - RX.pos(:,n)).^2)) ;
end

% figure, plot(scneario.distance_tx_rx), xlabel('N_{PRI}'), ylabel('Distance [m]'), title('distance TX RX')

% Speed distance-TX-RX
speed_distance_tx_rx = [0;diff(distance_tx_rx(:))];

speed_distance_tx_rx = speed_distance_tx_rx ./ const.PRI;

% figure, plot(speed_distance_tx_rx),hold on, plot(RX.speed) 
% legend("delta Distance / PRI", "Speed from log"),xlabel('N_{PRI}'), ylabel('Velocity [m/s]'), title('speed distance TX RX')

%% Targets distances
%cars

distance_cars = zeros(size(targets.cars,2),radar.N_PRI);
for n = 1:radar.N_PRI
    for car_idx = 1:size(targets.cars,2)
        distance_cars(car_idx,n) = sqrt(sum((TX.pos(:,n) - targets.cars(:,car_idx)).^2)) + sqrt(sum( (RX.pos(:,n) - targets.cars(:,car_idx)).^2));
    end
end

%% humans
distance_humans = zeros(size(targets.humans,1),radar.N_PRI);
for n = 1:radar.N_PRI
    for hum_idx = 1:size(targets.humans,2)
        distance_humans(hum_idx,n) = sqrt(sum((TX.pos(:,n) - targets.humans(:,hum_idx)).^2)) + sqrt(sum((RX.pos(:,n) - targets.humans(:,hum_idx)).^2));

    end
end

scenario.distance.tx_rx = distance_tx_rx;
scenario.speed_distance_tx_rx = speed_distance_tx_rx;
scenario.distance.humans = distance_humans;
scenario.distance.cars = distance_cars;
end

