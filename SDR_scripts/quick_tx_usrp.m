%% COSTANTS
SampleRate = 30e6;
fc=1.65e9;
disp("Begin setup!")


% Chirp Generation (Prof)
fs = SampleRate; % f campionamento
dt = 1/fs;
B = fs * 0.9; % banda totale
samples_per_chirp = 2^15;
T_tot = samples_per_chirp * dt; % durata del chirp
T = T_tot * .9;
K = B/T;
t = (-T/2:dt:T/2);

s = exp(1i*pi*K*t.^2);
s_pad=zeros(samples_per_chirp,1);
s_pad(1:length(s))=s;

% save(strcat("tx_waveform_B",num2str(B/1e6),"M_S",num2str(fs/1e6),"M"),"s_pad")
% to plot spectrum
S = fftshift(fft(s));
plot(abs(S))

%% Normalize (max in USRP is 1)

txNorm=0.99*s_pad/max(abs(s_pad));

%%%%%%%%%%
TxPow = 3;              % power in dBm
%%%%%%%%%%
MaxPow = 10;            % max power USRP ca transmit is 10dbm
MaxGain = 89.75;        % gain associated to max tx power
TxGain = MaxGain - (MaxPow - TxPow);
%% RADIO SETUP

tx = comm.SDRuTransmitter('Platform','B210', ...
          'SerialNum','3218D18', ...
          'Gain',TxGain,...  %min 0 max 89.75 (step 0.25)
          'CenterFrequency',fc, ...
          'MasterClockRate',fs, ...
          ...'SamplesPerFrame',frame_len,...
          ...'EnableBurstMode',1,...
          ...'NumFramesInBurst',burst_n,...
          'PPSSource','GPSDO',...
          'ClockSource','GPSDO',...
          'TransportDataType','int8',...
          'InterpolationFactor',1 ...
);
          

%% TRANSMISSION
disp("Begin Transmission!")     
while 1
    res = tx(txNorm);
    if res 
        disp("underrun detected")
    end
end
