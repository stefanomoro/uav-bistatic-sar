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
%% LOAD RAW
addpath("../lib");       % add path of lib
addpath("./utils")
allProcTimer = tic;
%% List Files
file_paths = listOnlyBinFiles(strcat(folder_name,'/raw'));
%%
for exp_num =4 : 5%length(file_paths)
file_path = file_paths(exp_num).complete_path;
disp(' '),disp(['Loading raw data ' num2str(exp_num)]),tic
A =load_bin(file_path(1:end-3)); 
disp(['Loaded in ' num2str(toc) ' s']);

%% RANGE COMPRESSION
disp("Processing FFT-based range compression")
timer = tic;

N_fft = 2 ^ nextpow2(samples_per_chirp);
N_slice = 1e3;
RC_slice = zeros(N_fft,N_slice,'like',A);

experiment_name = file_paths(exp_num).exp_name;
dest_file= strcat(folder_name,'/RC/complete/',experiment_name);
if(isfile(strcat(dest_file,".bb")))
    error("ERROR Destination file already present. Delete before continue")
end
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
    RC_slice = ifft( matFFT .* txWaveFFT,N_fft,1); 
    i = i+1;
    disp(['Slice ' num2str(i) ' done in ' num2str(toc) ' s'])
    %save the slice, appending to the end
    save_bin(dest_file,RC_slice)

end
disp(['Finished Range Compression in ' num2str(toc(timer)) ' s'])

clear N_fft matFFT txWaveFFT firstHalf A_rc_slice i slice_idx ss

end
