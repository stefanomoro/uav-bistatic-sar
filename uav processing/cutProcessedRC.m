function radar = cutProcessedRC(radar)
%CUTPROCESSEDRC cut raw data to process only a subset
%   radar = cutProcessedRC(radar)
radar.idx_wind_row = 150:300;%180:230;
%test2_1 = 1:4e4;test1_1 = 1:2.4e4
radar.idx_wind_col = 1:2.4e4;

radar.RC = radar.RC_raw(radar.idx_wind_row,radar.idx_wind_col);
radar.R_ax = radar.R_ax_raw(radar.idx_wind_row); 
radar.tau_ax = radar.tau_ax_raw(radar.idx_wind_col); 

radar.N_PRI = length(radar.tau_ax);

end

