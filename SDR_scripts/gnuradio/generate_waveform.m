%%
tx_wave = load("../../tx_waveform/tx_waveform_S56M.mat").s_pad;
tx_wave = single(tx_wave);

f = fopen("./tx_waveform_gnuradio_S56M.dat","w");
out = zeros(length(tx_wave) * 2,1,"single");
out(1:2:end-1) = real(tx_wave);
out(2:2:end) = imag(tx_wave);
fwrite(f,out,"float32");
fclose(f);