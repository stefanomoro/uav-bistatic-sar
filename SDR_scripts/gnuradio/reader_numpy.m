f = fopen('../../../giuriati_test/giuriati_gnuradio.dat', 'rb');
read_L = 1e9;
values = zeros(read_L,1,"single");
idx = 1;
while ~isempty(values) 
values= fread(f,read_L,"float32=>float32");

C = complex(values(1:2:end-1),values(2:2:end));


per = 122334;
L = floor(length(C)/per) * per;
A = C(1:L);
fs = 56e6;
PRI = per/fs;
A = reshape(A,per,[]);
tau_ax = 0:PRI:(size(A,2)-1)*PRI;

figure,imagesc(tau_ax,[],abs(A))
title(["plot " num2str(idx)])
idx = idx + 1;
end
fclose(f);
