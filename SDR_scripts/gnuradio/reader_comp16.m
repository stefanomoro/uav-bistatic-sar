f = fopen('../../mat_files/giuriati_test/22_03_18/test1.dat', 'rb');
read_L = 1e9;
values = zeros(read_L,1,"int16");
idx = 1;
y = zeros(0,1,'int16');
LL = 90 * 56e6;
while ~isempty(values) 

values= fread(f,read_L,"int16=>int16");
C = complex(values(1:2:end-1),values(2:2:end));
y = [y;C];
disp(['Read ' num2str(length(y)/LL*100) '%'])
end
fclose(f);

samples_per_chirp = 2^15;
L = floor(length(y) / samples_per_chirp) * samples_per_chirp;
A = reshape(y(1:L),samples_per_chirp,[]);
clear y
addpath(genpath([pwd, filesep, '..\..\lib' ]));       % add path of lib
A = single(A);
save_bin('../../mat_files/giuriati_test/22_03_18/test1',A);
