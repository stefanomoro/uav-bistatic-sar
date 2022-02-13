clear all
clc
close all
%% COSTANTS
chirp_bw = 27e6;                                % actual chirp bandwidth 
chirp_sr = 30e6;                                % sdr sample rate
experiment_name = strcat('MULTISTATIC_GPS1_B',...
    num2str(chirp_bw/1e6),'M_S',num2str(chirp_sr/1e6),'M');

file_name = strcat('binary files/radar_window_test/',...
    experiment_name,'.bb');
tx_wave = load(strcat('tx_waveform/tx_waveform_B',...
    num2str(chirp_bw/1e6),'M_S',num2str(chirp_sr/1e6),'M.mat')).s_pad;

samples_per_chirp = length(tx_wave);            % 22002 for 20 MSps, 33002 for 30MSps

bbr = comm.BasebandFileReader('Filename',file_name,...
        'SamplesPerFrame',samples_per_chirp);
bbr_info = info(bbr);
total_size = bbr_info.NumSamplesInData;

%% Align CHIRP 
Nchirp = 1e3;                                   % take first Nchirp for computing cutting point
temp = zeros(Nchirp * samples_per_chirp,1);
for i = 1:Nchirp
    temp((i - 1)*samples_per_chirp + 1 : (i*samples_per_chirp)) = bbr();
end

%% Inspect signal
% Y = fftshift(abs(fft(temp)));
% f = linspace(-.5,.5,length(Y)) * bbr.SampleRate;
% figure,plot(f,Y);title(" FFT original")
% figure, plot(abs(temp));title("Rx signal")

%% Correlate with TX CHIRP
[correlated ,lags]= xcorr(temp,tx_wave);
% correlated = abs(correlated )/ max(abs(correlated));
% figure, plot(lags,abs(correlated)),title("Correlation")

[~, idx ] = max(correlated);
% we know the idle time is 10% of the active time, so we leave 5% of idle
% at start
cut_point = lags(idx) - ceil(5 /110 *samples_per_chirp);

%% Test Matrix assemble
% align_temp = temp( cut_point: end);
% figure,plot(abs(align_temp)),title("Aligned rx sign")
% tt = cat(1,align_temp,zeros(ceil(length(align_temp) / samples_per_chirp) * samples_per_chirp - length(align_temp),1));
% A = reshape(tt, samples_per_chirp,[ ]);
% figure,imagesc(abs(A)),title("ABS matrix");
% 
% figure,imagesc(fftshift(abs(fft(A,[],2)))),title("FFT on tau")


clear Nchirp temp f i Y align_temp lags idx correlated tt
disp('Completed alignement of first chirp')

%% CONSTANTS for file reading
chirps_per_chunk  = 3e4;                % predefined value to have circa 1e9 samples per chunk
samples_per_chunk = samples_per_chirp * chirps_per_chunk;
y = zeros(samples_per_chunk,1,'int8');

N_dropped_samples = ...                 % variable used to store unuseful samples
    ceil(cut_point / samples_per_chirp) * samples_per_chirp;
dropped_samples = complex(zeros(N_dropped_samples,1,'int8'),...
    zeros(N_dropped_samples,1,'int8')) ;

%% DROP first samples
release(bbr);                           % reset the file reader pointer
i= 0;last_idx = 1;

while i < cut_point / samples_per_chirp
    x = bbr();
    dropped_samples(last_idx:last_idx + length(x)-1) = x;
    i = i +1;
    last_idx = last_idx + length(x);
end
% drop the first samples, but recover part of the first nice chirp
initial_shift = N_dropped_samples - cut_point;


%% Read  all samples

for chunk_idx = 1: ceil(total_size / samples_per_chunk)

y(1:initial_shift) = ...
    dropped_samples(length(dropped_samples)-(initial_shift) + 1 : end);

last_idx = initial_shift + 1;
i = 0;
tic
while last_idx <= (samples_per_chunk - (samples_per_chirp - initial_shift)) && not(isDone(bbr))
    x = bbr();
    y(last_idx:last_idx + length(x)-1) = x ;
    last_idx = last_idx + length(x);
    i = i + 1;
    if mod(i,6e3)== 0
        disp(['Elaborated ' num2str(i * 100 /chirps_per_chunk )  ' % of chunk']);
    end
end
dropped_samples = bbr();            % get last exeeding chirp
y(last_idx : last_idx -1  + (samples_per_chirp - initial_shift)) = ...
    dropped_samples (1: samples_per_chirp - initial_shift);

y = single(y);
toc
%% SAVING CHUNK
disp(' ')
disp(['Starting Saving chunk n ' num2str(chunk_idx)])
tic
save(strcat('mat_files/radar_window_test/',experiment_name,'/chunk',num2str(chunk_idx)),'y','-v7.3');
toc
disp(['Saved chunk n ' num2str(chunk_idx)])
disp(' ')

y = zeros(samples_per_chunk,1,'int8');  % reset array


end
