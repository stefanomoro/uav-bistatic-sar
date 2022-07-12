function const = initializeConstants(param)
%INITIALIZECONSTANTS define all the needed constants for the xcript
%   const = initializeConstants(param)
%   param contains the hardcoded parameter of the test

const.chirp_bw = param.chirp_bw;                      % actual chirp bandwidth (B)
const.chirp_sr = param.chirp_sr;                          % SDR sample rate (fs)
const.norm_B = const.chirp_bw / const.chirp_sr;
const.f0 = param.f0;
const.OSF = 8;

const.experiment_name = param.experiment_name;
const.last_mod_date_file = param.last_mod_date_file;

const.radar_folder_name = param.radar_folder_name;
const.drone_track_folder = param.drone_track_folder;
const.drone_track_file = param.drone_track_file;

const.tx_wave = load(strcat('..\tx_waveform/tx_waveform_S56M.mat')).s_pad;
const.tx_wave = single(const.tx_wave);
const.samples_per_chirp = length(const.tx_wave);

const.lambda = physconst('LightSpeed')/const.f0;
const.PRI = const.samples_per_chirp / const.chirp_sr;
const.PRF = 1 / const.PRI;
const.dt = 1/const.chirp_sr;
const.dR = (physconst('LightSpeed') * const.dt);

const.RC_file = strcat(const.radar_folder_name,const.experiment_name);

end

