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
tau_time = 60;%seconds
Nchirp = ceil(tau_time*samp_rate/samples_per_chirp);
chunk_len = Nchirp * samples_per_chirp;

A = reshape(A(chunk_len:end),samples_per_chirp,[]);
A = single(A);

%% PLOT

figure,imagesc(abs(A));


