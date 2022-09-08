clear all
clc
%close all
%%
addpath(genpath(['lib' ]));       % add path of lib

chirp_bw = 27e6;
samp_rate = 30e6;
file_n = input ("The file number ? ",'s');
file_name = strcat('radar_window_test/data',file_n');
% file_name = strcat('prova2/NOGPS_chirp10.bb');

if samp_rate == 30e6
    samples_per_chirp = 33002 ;
else
    samples_per_chirp = 11001 * samp_rate / 10e6;
end
A = load_bin(file_name);
L = samples_per_chirp * floor(length(A)/samples_per_chirp);
%%
tau_time = 5;%seconds
Nchirp = ceil(tau_time*samp_rate/samples_per_chirp);
chunk_len = Nchirp * samples_per_chirp;

A = reshape(A(1:chunk_len),samples_per_chirp,[]);
A = single(A);
max_V = max(max(real(A)));
disp(max_V)
%% PLOT
figure,imagesc(abs(A));

%% RANGE COMPRESSION
disp("Processing FFT-based range compression")

tx_wave = load("./tx_waveform/tx_waveform_B27M_S30M").s_pad;
N_fft = 2 ^ nextpow2(samples_per_chirp);
N_slice = 1e3;
A_rc = zeros(N_fft,0,'like',A);
i = 0;
txWaveFFT = repmat(conj(fft(tx_wave,N_fft,1)),[1,size(N_slice,2)]);
for slice_idx = 1:N_slice:size(A,2)-N_slice
    tic
    matFFT = fft(A(:,slice_idx : slice_idx+ N_slice - 1),N_fft,1);
    A_rc_slice = ifft( matFFT .* txWaveFFT,N_fft,1); 
    A_rc = cat(2,A_rc,A_rc_slice);
    i = i+1;
    disp(['Slice ' num2str(i) ' done in ' num2str(toc) ' s'])
end
disp(['Finished Range Compression'])
clear N_fft matFFT txWaveFFT firstHalf A_rc_slice i slice_idx ss

%% PLOT
if 0
dt = 1/samp_rate;
dR = (physconst('LightSpeed') * dt)/2;
R_margin = 5e2;                     % 1km of margin for RC
samp_margin = floor(R_margin/dR);

[~,peak_idx ]= max(A_rc(:,1e3));
plot_window = peak_idx-samp_margin:peak_idx+samp_margin;
if plot_window(1) <1
    A_rc = circshift(A_rc,-plot_window(1),1);
    plot_window = plot_window - plot_window(1);
elseif plot_window(end) > size(A_rc,1)
    A_rc = circshift(A_rc,size(A_rc,1) - plot_window(end),1);
    plot_window = plot_window - plot_window(end);
end
R_ax = -R_margin :dR:R_margin;

figure,imagesc([],R_ax,abs(A_rc(plot_window,:)))

end