SampleRate = 30e6;
fc=2.385e9;
snGTX='sn:1044735411960009090038004f371e8d27';
disp("Begin setup!")
% Old Chirp Generation
% t1=3e5;
% t=1:1:t1;
% f0=1;
% f1=4e6;
% s1= chirp(t,f0,t1,f1,'complex').';
% s=s1;
% 
% % Chirp Generation (Prof)
% fs = SampleRate; % f campionamento
% dt = 1/fs;
% B = fs/1.5; % banda totale
% T = 1e-3; % durata del chirp
% K = B/T;
% t = (-T/2:dt:T/2);
% Tc = T * 0.8;
% s = exp(1i*pi*K*t.^2);% .* rectpuls(t/Tc);
% s_pad=zeros(length(s)+ceil(length(s)*0.1),1);
% s_pad(1:length(s))=s;

% % to plot spectrum
% S = fftshift(fft(s));
% plot(abs(S))

% Normalize (check max in USRP?)
txWave = load("../tx_waveform/tx_waveform_B27M_S30M").s_pad;
% txNorm=transpose(0.5*s_pad/max(abs(s_pad)));
% S = fftshift(fft(s));
% plot(abs(S))



tx = comm.SDRuTransmitter('Platform','B210', ...
              'SerialNum','3218D18', ...
              'Gain',60,...  %min 0 max 89.75 (step 0.25)
              'CenterFrequency',fc, ...
              'MasterClockRate',fs, ...
              ...'SamplesPerFrame',frame_len,...
              ...%'EnableBurstMode',1,...
              ...%'NumFramesInBurst',burst_n,...
              'TransportDataType','int8',...
              'InterpolationFactor',1 ...
              );



disp("Begin Transmission!")     
while 1
    res = tx(txNorm);
    if res 
        disp("underrun detected")
    end
end
%transmitRepeat(tx,txNorm); %Non va sulle USRP
