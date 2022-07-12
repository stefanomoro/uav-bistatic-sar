function [radar,drone,targets] = loadData(const)
%LOADRAWRC load all the data from file and define useful parameters
%   [radar,drone,targets] = loadData(const)
radar.RC_raw = load_bin(const.RC_file);
radar.R_margin_raw = size(radar.RC_raw,1)*const.dR;
radar.R_ax_raw = -radar.R_margin_raw/2: const.dR: radar.R_margin_raw/2;
radar.tau_ax_raw = 0:const.PRI :const.PRI*size(radar.RC_raw,2);
radar.N_PRI_raw = length(radar.tau_ax_raw);


drone = load(const.drone_track_file).drone;
radar.last_mod_date = load(const.last_mod_date_file).last_mod_date;
targets = load(strcat(const.drone_track_folder,'target.mat'));
cars = [targets.cars.x, targets.cars.y, zeros(size(targets.cars.x))]';
humans = [targets.human.x, targets.human.y, zeros(size(targets.human.x))]';
targets = [];
targets.humans = humans;
targets.cars = cars;
end

