function [] = makeGIF(const,scenario,RX,TX,targets,focus,mode,c_axis)
%MAKEGIF Make gif with multiple focused image, Mode => 1 not eq, 2 eq
%   [] = makeGIF(const,scenario,RX,TX,targets,focus,mode,caxis)


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
fig = figure("WindowState","maximized");
F_vec = abs(F_vec);
idxs = [ 1:length(focus.angle_vec) length(focus.angle_vec)-1 :-1:2];
images = cell(size(idxs));

for i = 1:length(idxs)
    F =20*log10 (filterHammingFocus(F_vec(:,:,idxs(i)),3) );
    plotFocusedWithTargets(scenario,RX,TX,targets,F,...
        strcat("Squint angle ",num2str(focus.angle_vec(idxs(i))),"°" ));
    caxis(c_axis)
    frame = getframe(fig);
    images{i} = frame2im(frame);
end
% MAKE GIF
    filename = strcat(const.experiment_name,'.gif'); % Specify the output file name
    for idx = 1:length(images)
        [A,map] = rgb2ind(images{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.5);
        end
    end
end

