SampleRate = 30e6;
fc=2e9;
snTX='sn:1044739659930001100027008b3e294aaa';
disp("Begin setup!")

% Chirp Generation (Prof)
fs = SampleRate; % f campionamento
dt = 1/fs;
B = fs*.9; % banda totale
T = 1e-3; % durata del chirp
K = B/T;
t = (-T/2:dt:T/2);
Tc = T * 0.8;
s = exp(1i*pi*K*t.^2);% .* rectpuls(t/Tc);
s_pad=zeros(length(s)+ceil(length(s)*0.1),1);
s_pad(1:length(s))=s;

% % to plot spectrum
% S = fftshift(fft(s));
% plot(abs(S))


% Normalize (check max in USRP?)
txNorm=transpose(s_pad/max(abs(s_pad))).';

% Pluto Txer
GtxPluto = sdrtx('Pluto',...
       'RadioID',snTX,...
       'CenterFrequency',fc,...
       'Gain',-10,... %must be between -89 and 0
       ...'SamplesPerFrame',SPF,...
       'BasebandSampleRate',SampleRate);
% release(GtxPluto);
transmitRepeat(GtxPluto,txNorm);

