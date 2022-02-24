clear all
close all
clc
%%
fc = 2.2649e9;
fs=16e6;
dec=2;
frame_len=1e5;     %Samples per burst
% test_l=(0.1*fs)/dec;
% burst_n=ceil(test_l/frame_len);

%% RADIO SETUP
disp("SDR setup")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chiama findsdru per trovare la USRP  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx = comm.SDRuReceiver(...
              'Platform','B210', ...
              'SerialNum','31A3D0E', ...
              'Gain',10,...
              'CenterFrequency',fc, ...
              'MasterClockRate',fs, ...
              'SamplesPerFrame',frame_len,...
              ...%'EnableBurstMode',1,...
              ...%'NumFramesInBurst',burst_n,...
              'TransportDataType','int16',...000
              ...'OutputDataType','Single',... %in alternativa
              'DecimationFactor',dec...
              ...'PPSSource','GPSDO',...
              ...'ClockSource','GPSDO'...
              );




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_name = 'cubeSAT/test_21_11_16';
test_idx = 1;
name_string=strcat(folder_name,'/data',num2str(test_idx),'.bb');
while isfile(name_string)
    test_idx = test_idx + 1;
    name_string=strcat(folder_name,'/data',num2str(test_idx),'.bb'); %This is the name of the file that is saved!!
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rxWriter = comm.BasebandFileWriter(name_string, ...
         fs/dec, rx.CenterFrequency);
     
capture_time = 15; %Set the capture length
maxL = ceil(capture_time * (fs/dec) / frame_len); %Number of bursts to acquire (frame_len=1e5; <-- campioni per burst)


%% Option for test capture 
isTest = 0;

if isTest
    maxL = 0.05*maxL; %reduce the size of array to do only testing in ram
    prova_rx=zeros( maxL * frame_len,1,'int16');    %<<<< This saves to RAM, careful size!
else
    data=zeros(frame_len,1,'int16');
end

len=zeros(maxL,1,'single');
%% Acquisition

disp("Begin acquisition")
tic;
test_overrun=zeros(maxL,1);

i=1;
while(i<=maxL)      
   [data, acq_len, overrun] = rx();
   if overrun
       disp("overrun detected")
   end
   
    test_overrun(i)=overrun;
   if acq_len
    idxsav=((i-1)*frame_len+1):1:i*frame_len;
	if isTest
        prova_rx(idxsav) = data;             %<<<< This saves to RAM (careful size!)
    else
        rxWriter(data);                   %<<<< This saves to file
    end
    len(i)=acq_len;
    %disp(len(i))
    i=i+1;
   toc
   end
   
end
release(rx)

disp("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
%% END
% check overrun 
if find(test_overrun)
    disp('OVERRUN HAPPEND')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attenzione che crasha MATLAB VVVVVVVVVVVVVVVVVV %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isTest 
%     fs_fft = fs/dec;
%     Y = fftshift(abs(fft(prova_rx)));
%     figure
%     freq_ax = linspace(-fs_fft/2,fs_fft/2,length(prova_rx));
%     %plot(freq_ax,Y);
%     title("FFT")

    figure,plot(abs(single(prova_rx(1:1e6)))),title("time abs")
end

%%
% plot(1:1e5,real(data),1:1e5,imag(data))
% span = 1e6;
% figure
% plot(abs(single(prova_rx(1*span:2*span))))
