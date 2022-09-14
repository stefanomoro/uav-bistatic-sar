function radar = cutProcessedRC(radar,param)
%CUTPROCESSEDRC cut raw data to process only a subset
%   radar = cutProcessedRC(radar,param)

if(param.radar.idx_wind_row == -1) %we load all the data
    radar.idx_wind_row = 1:size(radar.RC_raw,1);
    radar.idx_wind_col = 1:size(radar.RC_raw,2);
else
    if length(param.radar.idx_wind_row) == 2
        radar.idx_wind_row = param.radar.idx_wind_row(1):param.radar.idx_wind_row(2);
        radar.idx_wind_col = param.radar.idx_wind_col(1):param.radar.idx_wind_col(2);
    else
        radar.idx_wind_row = param.radar.idx_wind_row;
        radar.idx_wind_col = param.radar.idx_wind_col;
    end
    
end

radar.RC = radar.RC_raw(radar.idx_wind_row,radar.idx_wind_col);
radar.R_ax = radar.R_ax_raw(radar.idx_wind_row); 
radar.tau_ax = radar.tau_ax_raw(radar.idx_wind_col); 

radar.N_PRI = length(radar.tau_ax);

end