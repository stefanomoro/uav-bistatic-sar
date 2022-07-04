function const = initializeConstants()
%INITIALIZECONSTANTS define all the needed constants for the xcript
%   

const.chirp_bw = 56e6*0.9;                      % actual chirp bandwidth (B)
const.chirp_sr = 56e6;                          % SDR sample rate (fs)
const.norm_B = const.chirp_bw / const.chirp_sr;
const.f0 = 1.65e9;
const.OSF = 8;

const.experiment_name = "test1_cut1";

const.radar_folder_name = '..\mat_files\uav_test\20220502\radar\test1\';
const.drone_track_folder = '..\mat_files\uav_test\20220502\drone track\';

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

