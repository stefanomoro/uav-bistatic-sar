clear all
close all
clc
%% COSTANTS
chirp_sr = 56e6;                                % SDR sample rate
chirp_bw = .9*chirp_sr;                         % actual chirp bandwidth 

folder_name = '../mat_files/uav_test/20220826';

tx_wave = load(strcat('../tx_waveform/tx_waveform_S56M.mat')).s_pad;
tx_wave = single(tx_wave);
samples_per_chirp = length(tx_wave);            % 2^15 mew, 33002 for 30MSps(old)

PRI = samples_per_chirp / chirp_sr;
PRF = 1 / PRI;
dt = 1/chirp_sr;
dR = (physconst('LightSpeed') * dt)/2;
R_margin = 500;                                 % meter of margin in RC                 
samp_margin = round(R_margin / dR);

addpath("../lib");       % add path of lib
addpath("./utils")
%% LOAD
addpath("../lib");       % add path of lib
addpath("./utils");
file_paths = listOnlyBinFiles(strcat(folder_name,'/RC/complete'));

exp_num = 5;
rc_file = file_paths(exp_num).complete_path;
RC = load_bin(rc_file(1:end-3));

figure,imagesc(abs(RC(:,1:100)))
%% Cut tau ax with empirical values
peak_idx = 21878;
%  [~,peak_idx] = max(RC(:,1));

% check if the RC cut with margin is wrapped and need to be circshifted
shift = 0;
if peak_idx + samp_margin > size(RC,1)
    shift = -samp_margin;
elseif peak_idx - samp_margin < 1
    shift = samp_margin;
end
if shift ~=0
    RC= circshift(RC,shift, 1);
    peak_idx = peak_idx + shift;
end
idx_window = peak_idx - samp_margin: peak_idx + samp_margin;
RC_cut = RC(idx_window,:);

result_file = strcat(folder_name,...
    '/RC/cut/',file_paths(exp_num).exp_name,'_cut');
disp(['Saving result of exp ' num2str(exp_num)]),tic
save_bin(result_file,RC_cut);
disp(['Saved in  ' num2str(toc) ' s']);
%% PLOTTING
shift_dist= 0;              % actual phisical zero position
RC_centered = circshift(RC_cut,floor(shift_dist/dR),1);

t_ax = linspace(0,size(RC_centered,1) * dt,size(RC_centered,1));
R_ax = linspace(- R_margin,R_margin,size(RC_centered,1));

tau_ax = linspace(0,PRI * size(RC_centered,2),size(RC_centered,2));

figure,imagesc(tau_ax, R_ax ,abs(RC_centered) ),
title('RC_{total}' ),xlabel("Slow time [s]"),ylabel("Range [m]")
