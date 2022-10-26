function [] = plotAllFocusedSquints(focus,scenario,RX,TX,targets,mode,col_axis)
%PLOTALLFOCUSEDSQUINTS plot all squints in 2 figures. mode 1:not EQ
% mode 2:EQ
%   [] = plotAllFocusedSquints(focus,scenario,RX,TX,targets,mode,col_axis)
switch mode
    case 1
        %plot not equalized focused image
        F_vec = focus.Focused_vec;
    case 2
        %plot equalized image
        F_vec = focus.Focus_eq;
    otherwise
        return
end

plot_per_wind = 6;
N_wind = ceil(size(F_vec,3)/plot_per_wind);

for wind_i = 1:N_wind
    figure("WindowState","maximized");
    for idx = 1:plot_per_wind

        subplot(2,3,idx)
        real_idx = idx + (wind_i-1) * plot_per_wind;
        if(real_idx >size(F_vec,3))
            return
        end
        F = filterHammingFocus(F_vec(:,:,real_idx),3);

        plotFocusedWithTargets(scenario,RX,TX,targets,20*log10(F),...
        strcat("Squint ",num2str(focus.angle_vec(real_idx)) ,"°, Hamming, dB"));  
        caxis(col_axis)
    end
end

end

