% clear all
close all
clc
%% COSTANTS
chirp_bw = 27e6;                                % actual chirp bandwidth 
chirp_sr = 30e6;                                % SDR sample rate
experiment_name = strcat('MULTISTATIC_GPS1_B',...
    num2str(chirp_bw/1e6),'M_S',num2str(chirp_sr/1e6),'M');

mat_file = strcat('mat_files/radar_window_test/',...
    experiment_name,'/chunk');

tx_wave = load(strcat('tx_waveform/tx_waveform_B',...
    num2str(chirp_bw/1e6),'M_S',num2str(chirp_sr/1e6),'M.mat')).s_pad;
tx_wave = single(tx_wave);
samples_per_chirp = length(tx_wave);            % 22002 for 20 MSps, 33002 for 30MSps

PRI = samples_per_chirp / chirp_sr;
PRF = 1 / PRI;
dt = 1/chirp_sr;
dR = (physconst('LightSpeed') * dt)/2;

%% LOAD chunk
% chunk_idx = 4;
allProcTimer = tic;
for chunk_idx = 1:6
cicleTimer = tic;
    disp(' '),disp(['Loading chunk n ' num2str(chunk_idx)]),tic
load( [mat_file num2str(chunk_idx)] ); 
disp(['Loaded in ' num2str(toc) ' s']);
A = reshape(y,samples_per_chirp,[]);            % reshape as matrix of t X tau 
clear y
%% DF computing

% disp ("DF computing"),tic
% 
% Nchirps_small = 10e3;                           % small set of chirps only for computing the freq offset
% A_small = A(:,1 :Nchirps_small);
% A_fft = fftshift(fft(A_small,2^nextpow2(size(A_small,2)),2),2);
% [~, idx] = max(max(abs(A_fft)));
% slow_freq = linspace(-PRF/2,PRF/2,size(A_fft,2));
% DF(chunk_idx) = slow_freq(idx);
% disp( ['Found DF = ' num2str(slow_freq(idx)) ' , in ' num2str(toc) ' s']);
% clear A_small Nchirps_small idx

%% DF correction
% disp("Starting correction"),tic
% if chunk_idx == 1
%     T_start= 0;
% end
% t = single(linspace(T_start,T_start + dt * numel(A),numel(A)));
% T_start = t(end) + dt;
% t = reshape(t,size(A));
% fixer = single(exp(-1i*2*pi*DF(1)*t));           % use always the first DF
% A_fixed= A.* fixer;
% disp(['Finished correction in ' num2str(toc) ' s'])

% figure,imagesc(slow_freq,[],abs(A_fft)),title("FFT before correction")
% figure,imagesc(slow_freq,[],fftshift(abs(fft(A_fixed,length(slow_freq),2)) )) ,title("FFT after correction")
% clear  A_fft slow_freq fixer t
A_fixed = A;
%% RANGE COMPRESSION
disp("Processing FFT-based range compression")
timer = tic;

N_fft = 2 ^ nextpow2(samples_per_chirp);
N_slice = 1e3;
A_rc = zeros(N_fft,0,'like',A);
i = 0;
txWaveFFT = repmat(conj(fft(tx_wave,N_fft,1)),[1,size(N_slice,2)]);
for slice_idx = 1:N_slice:size(A,2)
    tic
    matFFT = fft(A_fixed(:,slice_idx : slice_idx+ N_slice - 1),N_fft,1);
    A_rc_slice = ifft( matFFT .* txWaveFFT,N_fft,1); 
    A_rc = cat(2,A_rc,A_rc_slice);
    i = i+1;
    disp(['Slice ' num2str(i) ' done in ' num2str(toc) ' s'])
end
disp(['Finished Range Compression in ' num2str(toc(timer)) ' s'])

clear N_fft matFFT txWaveFFT firstHalf A_rc_slice i slice_idx ss
%% SAVING
addpath(genpath([pwd, filesep, 'lib' ]));    % add path of lib
disp(' '),disp(['Start saving chunk n ' num2str(chunk_idx)])
tic
dest_file= strcat('mat_files/radar_window_test/',experiment_name,'/RC/chunk',num2str(chunk_idx));
% save(dest_file,'A_rc','-v7.3','-nocompression');
save_bin(dest_file,A_rc)
disp(['Saved chunk n ' num2str(chunk_idx) ' in ' num2str(toc) ' s'])
disp(['Completed chunk n ' num2str(chunk_idx) ' in ' num2str(toc(cicleTimer) / 60) ' minutes'])
disp(' ')
end
disp(' '),disp(['All processing took ' num2str(toc(allProcTimer)/ 60) ' minutes']);
% 
% dest_file= strcat('mat_files/radar_window_test/',experiment_name,'/RC/DF'); % save DF array
% save(dest_file,'DF','-v7.3');

%% CHECK Peak position for each chunk 
addpath(genpath([pwd, filesep, 'lib' ]));       % add path of lib

L_tau_chunk = 30e3;                             % tau ax length for each chunk
N_chunks = 6;
control_tau = 1:5e3:L_tau_chunk;                % tau points used for max checking
peak_t_idx = zeros(length(control_tau),N_chunks); 

for chunk_idx = 1:N_chunks
disp(['Loading chunk n ' num2str(chunk_idx)]),tic
A_chunk = load_bin( [mat_file num2str(chunk_idx)] ); 
disp(['Loaded in ' num2str(toc) ' s']);
for i = 1 : length(control_tau)
    [~,peak_t_idx(i,chunk_idx)] = max(A_chunk(:,control_tau(i)));
end
end
clear A_chunk i 


%% FIND stable RC 
diff_peak = [diff(peak_t_idx(:));0];
good_diff = abs(diff_peak) < 50;

counter = 0; good = 0;
st_idx = [];end_idx = [];
% cycle to find stable peak over at least 5 control points
for i = 1:length(good_diff)
    if good_diff(i)
        counter =counter + 1;
        if not(good)
            st_idx = [st_idx i];
        end
        good = 1;
    else 
        if counter > 4
            end_idx = [end_idx (i)];
        elseif good 
            st_idx = st_idx(1:end-1);
        end
        counter = 0; good = 0;
    end
end

if length(st_idx) > length(end_idx)
    st_idx = st_idx(1:end-1);
end
% mask used select the valuable control points of each chunk
mask = zeros(size(good_diff));
mask(st_idx:end_idx) = 1;
mask = reshape(mask,size(peak_t_idx)) == 1;

clear diff_peak good_diff counter good
%% PREPARE Cut of RC with Margin
good_peak_t = peak_t_idx.*mask;

% check if the RC cut with margin is wrapped and need to be circshifted
good_peak_idx = peak_t_idx(mask);good_peak_idx = good_peak_idx(1);
shift = 0;
if good_peak_idx + samp_margin > N_input_lags
    shift = -samp_margin;
elseif good_peak_idx - samp_margin < 1
    shift = samp_margin;
end
%% KEEP good chunk 
good_chunk = sum(good_peak_t,1) > 0 ;
RC_total = zeros(samp_margin*2+1,0,'single');
for chunk_idx = 1:N_chunks
    st_idx = 0; end_idx = 0; 
    if good_chunk(chunk_idx)
        A_chunk = load_bin( [mat_file num2str(chunk_idx)] ); 
        tau = control_tau'.* (good_peak_t(:,chunk_idx)>0 );
        if tau(1) == 0              
            % starting point
            ii = find(tau); ii = ii(1);
            st_idx = tau(ii);
            A_chunk = A_chunk(:,st_idx:end);
        elseif tau(end) == 0        
            % ending chunk
            ii = find(tau); ii = ii(end);
            end_idx = tau(ii);
            A_chunk = A_chunk(:,1:end_idx);
        end
        if shift ~=0
            A_chunk = circshift(A_chunk,shift, 1);
            good_peak_sh_idx = good_peak_idx + shift;
        end
        cut_rc = good_peak_sh_idx - samp_margin: good_peak_sh_idx + samp_margin;
        RC_total = cat(2,RC_total,A_chunk (cut_rc,:));

    end
end


result_file = strcat('mat_files/radar_window_test/',...
    experiment_name,'/RC/RC_total');
disp(['Saving result of tau L => ' num2str(size(RC_total,2)) ]),tic
save_bin(result_file,RC_total);
disp(['Saved in  ' num2str(toc) ' s']);
clear A_chunk
%% PLOTTING
shift_dist= 57;              % actual zero is at 60m
RC_centered = circshift(RC_total,floor(shift_dist/dR),1);

t_ax = linspace(0,size(RC_centered,1) * dt,size(RC_centered,1));
R_ax = linspace(- R_margin,R_margin,size(RC_centered,1));

tau_ax = linspace(0,PRI * size(RC_centered,2),size(RC_centered,2));

figure,imagesc(tau_ax, R_ax ,abs(RC_centered) ),
title('RC_{total}' ),xlabel("Slow time [s]"),ylabel("Range [m]")
