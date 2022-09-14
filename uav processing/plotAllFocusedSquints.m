function [] = plotAllFocusedSquints(focus,mode,col_axis)
%PLOTALLFOCUSEDSQUINTS plot all squints in 2 figures. mode 1:not EQ
% mode 2:EQ
%   [] = plotAllFocusedSquints(focus,mode,col_axis)
switch mode
    case 1
        %plot not equalized focused image
        F = focus.Focused_vec;
    case 2
        %plot equalized image
        F = focus.Focus_eq;
    otherwise
        return
end
figure("WindowState","maximized");
for idx = 1:8
subplot(2,4,idx)
F = filterHammingFocus(F(:,:,idx),3);
plotFocusedWithTargets(scenario,RX,TX,targets,20*log10(F),...
        strcat("Squint ",num2str(focus.angle_vec(idx)) ,"°, dB"));  
    caxis(col_axis)
end
figure("WindowState","maximized");
for idx = 8:15
subplot(2,4,idx-7)
F = filterHammingFocus(F(:,:,idx),3);
plotFocusedWithTargets(scenario,RX,TX,targets,20*log10(F),...
        strcat("Squint ",num2str(focus.angle_vec(idx)) ,"°, dB"));  
    caxis(col_axis)
end

end

