%% ======================================================================== INTERPOLATION
OSF = 8; 

t = 0:size(RC,1) -1;
tf = 0:1/OSF:(size(RC,1) -1);

tf_men_t =tf(:) - t;
W = sinc(norm_B*tf_men_t);

% interpolate whole RC signal
RC_OS = W*RC;
R_ax_OS = R_ax(1):dR/OSF:R_ax(end);

% figure
% imagesc(tau_ax,R_ax_OS,abs(RC_OS))
% title('RC interp - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

clear t tf tf_men_t W OSF
%%  =======================================================================  TIME SHIFT CORRECTION
%% Find correct idxs cross-talk
% define the idxs where the cross talk should be

cross_talk_idxs_corr = zeros(size(tau_ax));
for n = 1:length(tau_ax)
    idxs = find(R_ax_OS > distance_tx_rx(n));
    cross_talk_idxs_corr(n) = idxs(1);
end

cross_talk_idxs_corr = movmean(cross_talk_idxs_corr,1e3);
clear idxs
%% Cross-talk tracking
[~,cross_talk_idxs] = max(RC_OS);
% apply median filter and then moving average
cross_talk_idxs = movmean(medfilt1(cross_talk_idxs,1e3) ,1000); 

%% TIME SHIFT CORRECTION
cross_talk_shifts = cross_talk_idxs - cross_talk_idxs_corr;

Nf = 2^nextpow2(size(RC_OS,1));
X = fftshift(fft(RC_OS,Nf,1),1);

f = (-Nf/2:Nf/2 - 1)/Nf;
H = exp(1i*2*pi*f(:)*cross_talk_shifts);
H(1,:) = 0;                         

RC_Dt_fixed = ifft(ifftshift(X.* H,1),Nf,1);
RC_Dt_fixed = RC_Dt_fixed(1:size(RC_OS,1),:);

% figure
% imagesc(tau_ax,R_ax,abs(RC_Dt_fixed))
% title('RC Dt fixed - abs plot' ),xlabel("Slow time [s]"),ylabel("Range [m]")

clear Nf X f H
%%  ======================================================================= PHASE CORRECTION
%% find cross-talk value
cross_talk = max(RC_Dt_fixed);
cross_talk = movmean(cross_talk,1e3); 
%% Find correct cross talk phase
% Correct phase
c = physconst('Lightspeed');
phase_corr_old = -2*pi*f0*(distance_tx_rx/c); % - 2*pi*f0*(2*speed_distance_tx_rx/lambda);
phase_corr_old = phase_corr_old(:).';

% move first phase in 2pi interval
phase_0 = mod(phase_corr_old(1),2*pi) * 2*pi -pi;

phase_corr = phase_corr_old - phase_corr_old(1) + phase_0;
% 
% figure, plot(tau_ax, unwrap(angle(cross_talk))),
% hold on, plot(tau_ax, phase_corr),  xlabel('tau axis [s]'), ylabel('Phase axis [Hz]')
% legend('cross talk phase','phase corr')
% title('Cross-talk phase')
clear phase_0 phase_corr_old
%%   ================================================== PHASE CORRECTION BRUTAL
% %% Df FIX 
% computed_phase = movmean(unwrap(angle(cross_talk)),1e3);
% phase_shift = computed_phase(:) - phase_corr(:);
% phasor = exp(-1i*phase_shift);
% peak_fixed = cross_talk(:).* phasor(:);
% 
% %% FULL Df fix
% phasor_mat = repmat(phasor(:).',size(RC_Dt_fixed,1),1);
% RC_Df_fixed = RC_Dt_fixed .* phasor_mat;


%%  ======================================================================= PHASE CORRECTION Correct version

%% Cycle 
%number of shifts of the window over the peak
win_size = 2^9; %1024;
N_cycle = floor(length(cross_talk)/win_size) * 2 - 1;
Df = zeros(N_cycle,1);
N_fft = 8*win_size;
freq_ax = (-N_fft/2:N_fft/2-1)/N_fft * PRF;
win_idx = 1 : win_size;

check_win_size = zeros(N_cycle,1);
for i = 1:N_cycle
    aa = cross_talk(win_idx);
    ff = fftshift(fft(aa,N_fft));
    [~,ii] = max(abs(ff));
    Df(i) = freq_ax(ii);
    win_idx = win_idx + win_size/2;
    
    check_win_size(i) = max(abs(ff))/sum(abs(aa));
end

figure,plot(Df)
title('Df from fft')
xlabel('Slow time [samples]')


% figure,plot(check_win_size)
% title('check window size')
% xlabel('Slow time [samples]')

clear N_cycle N_fft freq_ax win_idx check_win_size aa
%% INTEGRATE Df to get phase
computed_phase = [];
pp = angle(cross_talk(1));
t = 0:PRI:(win_size/2-1)*PRI;
for i = 1:length(Df)
    pp = pp(end) + 2*pi*Df(i)*t(:);
    computed_phase = [computed_phase; pp];
end
if length(computed_phase)<length(cross_talk)
    temp = repmat(computed_phase(end),length(cross_talk)-length(computed_phase),1);
    computed_phase = [computed_phase(:); temp]; 
end
% computed_phase = movmean(computed_phase,1e3);

% figure,plot(unwrap(angle(cross_talk))),hold on, plot(computed_phase),  
% legend('unwrap angle peak','computed phase'), 
% title('Check if phase value is correct'), xlabel('Slow time [samples]')
clear temp t pp
%% Df FIX 
phase_shift = computed_phase(:) - phase_corr(:);
phasor = exp(-1i*phase_shift);
peak_fixed = cross_talk(:).* phasor(:);

% figure,plot(computed_phase),hold on,plot(phase_corr),plot(phase_shift),legend("Computed phase","Phase correct","Phase shift")
%% FULL Df fix
phasor_mat = repmat(phasor(:).',size(RC_Dt_fixed,1),1);
RC_Df_fixed = RC_Dt_fixed .* phasor_mat;

clear phasor
%%  ======================================================================= FINAL CORRECTION RESULTS
%% Amplitude
figure,subplot(2,1,1)
imagesc(tau_ax,R_ax_OS,abs(RC_OS));
title('RC original' ),xlabel("Slow time [s]"),ylabel("Range [m]")
subplot(2,1,2)
imagesc(tau_ax,R_ax_OS,abs(RC_Df_fixed));
title('RC fixed' ),xlabel("Slow time [s]"),ylabel("Range [m]")
%% Phase 
figure,subplot(2,1,1)
imagesc(tau_ax,R_ax_OS,angle(RC_OS));
title('RC Dt fixed' ),xlabel("Slow time [s]"),ylabel("Range [m]")
subplot(2,1,2)
imagesc(tau_ax,R_ax_OS,angle(RC_Df_fixed));
title('RC freq fixed' ),xlabel("Slow time [s]"),ylabel("Range [m]")


%% Peak Phase
figure,plot(angle(cross_talk)),hold on,plot(angle(peak_fixed)); 
legend("Original","Df fixed"),xlabel("Slow time [samples]"), title('Crosstalk RC')

figure,plot(unwrap(angle(cross_talk))),hold on,plot(unwrap(angle(peak_fixed))),plot(phase_corr) 
legend("Original unwrap","Df fixed unwrap","Distance based phase"),xlabel("Slow time [samples]"), title('Crosstalk RC')



%% ======================================================================== OUTPUT
RC = RC_Df_fixed;
R_ax = R_ax_OS;

