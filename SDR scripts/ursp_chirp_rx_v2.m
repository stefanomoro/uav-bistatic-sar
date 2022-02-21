clear all
close all
clc
%% CONSTANTS
addpath(genpath('lib'));       % add path of lib

fc = 2e9;
fs=60e6;
dec=2;
frame_len=1e4;     %Samples per burst


%% RADIO SETUP
disp("SDR setup")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chiama findsdru per trovare la USRP  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx = comm.SDRuReceiver(...
              'Platform','B210', ...
              'SerialNum','31A3D0E', ...
              'Gain',40,...   % max gain 76
              'CenterFrequency',fc, ...
              'MasterClockRate',fs, ...
              'SamplesPerFrame',frame_len,...
              ...%'EnableBurstMode',1,...
              ...%'NumFramesInBurst',burst_n,...
              'TransportDataType','int8',...000
              ...'OutputDataType','Single',... %in alternativa
              'DecimationFactor',dec,...
              'PPSSource','GPSDO',...
              'ClockSource','GPSDO'...
              );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_name = 'radar_window_test';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_idx = 1;
name_string=strcat(folder_name,'/data',num2str(test_idx),'.bb');
while isfile(name_string)
    test_idx = test_idx + 1;
    name_string=strcat(folder_name,'/data',num2str(test_idx),'.bb'); %This is the name of the file that is saved!!
end

% rxWriter = comm.BasebandFileWriter(name_string, ...
%          fs/dec, rx.CenterFrequency);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
capture_time = 120; %Set the capture length in s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxL = ceil(capture_time * (fs/dec) / frame_len); %Number of bursts to acquire (frame_len=1e5; <-- campioni per burst)


%% RAM capture
RAMsave=complex(zeros(maxL,frame_len,'int8'),zeros(maxL,frame_len,'int8'));
data = complex(zeros(frame_len,1,'int8'),zeros(frame_len,1,'int8'));

len=zeros(maxL * 2,1);
test_overrun=zeros(maxL*2,1);
times_store=zeros(maxL,1);
%% Acquisition

disp("Begin acquisition")
start_t = datetime("now");
end_t = datetime("now") + seconds(capture_time);
disp("Acquisition end at ")
disp(end_t)
disp("")

i = 1;
j = 1;
tic;
while i <= maxL
   [data, len(j), test_overrun(j)] = rx();
   if len(j) 
        RAMsave(i,:)=data;  %<<<< This saves to RAM
        times_store(i)=toc;
        tic;
        i = i + 1;
   end
   j = j+1;
end
end_t = datetime("now");
release(rx)
toc   
disp("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
disp("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
disp("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
%% END
% check overrun 
if find(test_overrun)
    disp('OVERRUN HAPPEND')
    figure,stem(test_overrun),hold on, 
    stem(len/(frame_len * 2)),
    stem(times_store*5)
    legend("OVERRUN","CaptureLenght","times");
end

%% Save to file
disp('Start saving')
RAMsave = reshape(RAMsave.',[],1);
save_bin(name_string(1:end-3),RAMsave);
disp('Save complete')
timestamp.start = start_t;
timestamp.end = end_t;
save(strcat(name_string(1:end-3),'_timestamp'),"timestamp");


% disp('Start saving')
% for ii = 1:maxL
%     rxWriter(RAMsave(ii,:).');                   %<<<< This saves to file
%     if mod(ii,floor(maxL/10)) == 0
%         disp(['Saved ' num2str(ii) ' of ' num2str(maxL) ' frames'])
%     end
% end
% disp('Save complete')