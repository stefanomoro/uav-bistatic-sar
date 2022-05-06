file_path = '../../mat_files/uav/20220502/raw/test1.dat';
last_mod_date = datetime(dir(file_path).date,'Locale','it_IT','InputFormat','dd-MMM-yyyy HH:mm:ss','Format','yyyy-MM-dd HH:mm:ss.SSS');
f = fopen(file_path, 'rb');

read_L = 1e9;                       % length to be read each iteration
values = zeros(read_L,1,"int16");
idx = 1;
y = zeros(0,1,'int16');
LL = dir(file_path).bytes /4;          % num complex samples
while ~isempty(values) 

values= fread(f,read_L,"int16=>int16");
C = complex(values(1:2:end-1),values(2:2:end));
y = [y;C];
disp(['Read ' num2str(length(y)/LL*100) '%'])
end
fclose(f);

disp('Start matrix reshape...')
samples_per_chirp = 2^15;
L = floor(length(y) / samples_per_chirp) * samples_per_chirp;
A = reshape(y(1:L),samples_per_chirp,[]);
clear y
disp('Done')
addpath(genpath([pwd, filesep, '..\..\lib' ]));       % add path of lib
A = single(A);
disp('Start saving...')
save_bin(file_path(1:end-4),A);
save([file_path(1:end-4) '_last_mod_date'],'last_mod_date')
disp('Done')
