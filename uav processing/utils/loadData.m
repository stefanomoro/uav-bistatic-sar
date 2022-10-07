    function [radar,drone,targets] = loadData(const)
%LOADRAWRC load all the data from file and define useful parameters
%   [radar,drone,targets] = loadData(const)
radar.RC_raw = load_bin(const.RC_file);
radar.R_margin_raw = size(radar.RC_raw,1)*const.dR;
radar.R_ax_raw = -radar.R_margin_raw/2: const.dR: radar.R_margin_raw/2;
radar.tau_ax_raw = 0:const.PRI :const.PRI*size(radar.RC_raw,2);
radar.N_PRI_raw = length(radar.tau_ax_raw);

temp = load(const.drone_track_file);
if(isfield(temp,"drone"))       %backward compatibility
    drone = temp.drone;
end

if(strcmp(const.setup_mode,"bistatic"))
    drone.rx = temp.gps_matrice;
    if(strcmp(const.tx_mode,"drone"))
        drone.tx = load(const.drone_tx_track_file);
        drone.tx = drone.tx.tarot;
    elseif( strcmp(const.tx_mode,"fixed"))
        drone.tx = load(const.drone_tx_track_file);
        drone.tx = drone.tx.stable_gps;
    end
elseif (strcmp(const.setup_mode,"monostatic"))
    % TODO monostatic case
    drone = temp.gps_matrice;        
end

radar.last_mod_date = load(const.last_mod_date_file).last_mod_date;
targets = load(strcat(const.drone_track_folder,'targets.mat'));
if(isfield(targets,"cars"))
    cars = [targets.cars.x, targets.cars.y, zeros(size(targets.cars.x))]';
    humans = [targets.human.x, targets.human.y, zeros(size(targets.human.x))]';
    targets = [humans,cars];
else
    targets = [targets.targets.utm_x(:),targets.targets.utm_y(:),...
        zeros(length(targets.targets.utm_x),1)].';
end
end

