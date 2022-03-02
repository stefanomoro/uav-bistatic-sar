f = fopen('../../../gnuradio_data/gnuradio_56.dat', 'rb');
read_L = 1e9;
idx = 1;
values= fread(f,read_L,"float32=>float32");
C = complex(values(1:2:end-1),values(2:2:end));
fclose(f);


%%
tx_wave = load("../../tx_waveform/tx_waveform_B27M_S30M.mat").s_pad;
tx_wave = single(tx_wave);

f = fopen("./tx_waveform_gnuradio_v2.dat","w");
out = zeros(length(tx_wave) * 2,1,"single");
out(1:2:end-1) = real(tx_wave);
out(2:2:end) = imag(tx_wave);
fwrite(f,out,"float32");
fclose(f);