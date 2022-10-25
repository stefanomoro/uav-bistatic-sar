folder = '../../mat_files/giuriati_test/20221020/raw/gnuradio/';
dir_out = dir(folder);
if(length(dir_out) == 2)
    disp("No file presents")
    return ;
end
dir_out = dir_out(3:end);

for i =1 : length(dir_out)
    if(dir_out(i).isdir)
        continue
    end
%%
file_path = strcat(folder,dir_out(i).name);
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
%%
disp('Start matrix reshape...')
samples_per_chirp = 2^15;
L = floor(length(y) / samples_per_chirp) * samples_per_chirp;
A = reshape(y(1:L),samples_per_chirp,[]);
clear y
disp('Done')
%%
addpath(genpath([pwd, filesep, '..\..\lib' ]));       % add path of lib
A = single(A);
disp('Start saving...')
save_bin(file_path(1:end-4),A);
save([file_path(1:end-4) '_last_mod_date'],'last_mod_date')
disp(strcat('Done ',num2str(i), "/",num2str(length(dir_out))))

end