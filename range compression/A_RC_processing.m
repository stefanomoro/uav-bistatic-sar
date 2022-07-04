clear all
close all
clc
%% COSTANTS
chirp_sr = 56e6;                                % SDR sample rate
chirp_bw = .9*chirp_sr;                         % actual chirp bandwidth 

experiment_name = 'test';
folder_name = 'mat_files/uav_test/20220502/';

tx_wave = load(strcat('tx_waveform/tx_waveform_S56M.mat')).s_pad;
tx_wave = single(tx_wave);
samples_per_chirp = length(tx_wave);            % 2^15 mew, 33002 for 30MSps(old)

PRI = samples_per_chirp / chirp_sr;
PRF = 1 / PRI;
dt = 1/chirp_sr;
dR = (physconst('LightSpeed') * dt)/2;
R_margin = 500;                                 % meter of margin in RC                 
samp_margin = round(R_margin / dR);
%% LOAD chunk
addpath(genpath([pwd, filesep, 'lib' ]));       % add path of lib

allProcTimer = tic;

for exp_num = 3:3
disp(' '),disp(['Loading raw data ' num2str(exp_num)]),tic
A =load_bin(strcat(folder_name,'raw/',experiment_name,num2str(exp_num))); 
disp(['Loaded in ' num2str(toc) ' s']);
%% RESHAPE 
% end_idx = samples_per_chirp * floor(length(A)/samples_per_chirp);
% A = A(1:end_idx);
% A = single(reshape(A,samples_per_chirp,[]));            % reshape as matrix of t X tau 
disp('Completed alignement of first chirp')
clear y Nchirp temp lags idx correlated start_idx end_idx

%% RANGE COMPRESSION
disp("Processing FFT-based range compression")
timer = tic;

N_fft = 2 ^ nextpow2(samples_per_chirp);
N_slice = 1e3;
RC_cut = zeros(N_fft,0,'like',A);
i = 0;
txWaveFFT = repmat(conj(fft(tx_wave,N_fft,1)),[1,size(N_slice,2)]);
for slice_idx = 1:N_slice:size(A,2)
    tic
    if slice_idx+ N_slice - 1 > size(A,2)
        matFFT = fft(A(:,slice_idx : end),N_fft,1);
        txWaveFFT = repmat(conj(fft(tx_wave,N_fft,1)),[1,size(matFFT,2)]);
    else 
        matFFT = fft(A(:,slice_idx : slice_idx+ N_slice - 1),N_fft,1);
    end
    A_rc_slice = ifft( matFFT .* txWaveFFT,N_fft,1); 
    RC_cut = cat(2,RC_cut,A_rc_slice);
    i = i+1;
    disp(['Slice ' num2str(i) ' done in ' num2str(toc) ' s'])
end
disp(['Finished Range Compression in ' num2str(toc(timer)) ' s'])

clear N_fft matFFT txWaveFFT firstHalf A_rc_slice i slice_idx ss
%% SAVING

disp(' '),disp('Start saving RC data')
tic
dest_file= strcat(folder_name,'RC/complete/',experiment_name,num2str(exp_num));

save_bin(dest_file,RC_cut)
disp(['Saved RC data of experiment ' experiment_name num2str(exp_num)])
disp(' ')
end
%% INSPECTION
% plot raw data
N_plot = 1e4;
for i = 1:floor(size(A,2)/N_plot)
    figure,imagesc(abs(A(:,(i-1)*N_plot +1:i*N_plot)))
end
%%
% plot RC data
N_plot = 1e4;
for i = 1:floor(size(RC_cut,2)/N_plot)
    figure,imagesc(abs(RC_cut(:,(i-1)*N_plot +1:i*N_plot)))
end

return
%% CUT contiguos partition
addpath(genpath([pwd, filesep, 'lib' ]));       % add path of lib

exp_num = 3;
RC_file= strcat(folder_name,'RC/complete/',experiment_name,num2str(exp_num));
RC = load_bin(RC_file);

%figure,imagesc(abs(RC))
start_idx = [1];
end_idx = [size(RC,2)];
%% Cut tau ax with empirical values
for cut_num = 1:length(end_idx)

[~,peak_idx] = max(abs(RC(:,end_idx(cut_num))));
RC_cut = RC(:,start_idx(cut_num):end_idx(cut_num));
% check if the RC cut with margin is wrapped and need to be circshifted
shift = 0;
if peak_idx + samp_margin > size(RC_cut,1)
    shift = -samp_margin;
elseif peak_idx - samp_margin < 1
    shift = samp_margin;
end
if shift ~=0
    RC_cut= circshift(RC_cut,shift, 1);
    peak_idx = peak_idx + shift;
end
idx_window = peak_idx - samp_margin: peak_idx + samp_margin;
RC_cut = RC_cut(idx_window,:);

result_file = strcat(folder_name,...
    '/RC/cut/',experiment_name,num2str(exp_num),'_cut',num2str(cut_num));
disp(['Saving result of exp ' num2str(exp_num) ' cut ' num2str(cut_num)]),tic
save_bin(result_file,RC_cut);
disp(['Saved in  ' num2str(toc) ' s']);
end
%% PLOTTING
shift_dist= 0;              % actual phisical zero position
RC_centered = circshift(RC_cut,floor(shift_dist/dR),1);

t_ax = linspace(0,size(RC_centered,1) * dt,size(RC_centered,1));
R_ax = linspace(- R_margin,R_margin,size(RC_centered,1));

tau_ax = linspace(0,PRI * size(RC_centered,2),size(RC_centered,2));

figure,imagesc(tau_ax, R_ax ,abs(RC_centered) ),
title('RC_{total}' ),xlabel("Slow time [s]"),ylabel("Range [m]")
