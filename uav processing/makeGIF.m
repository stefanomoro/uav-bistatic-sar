function [] = makeGIF(const,scenario,RX,TX,targets,focus)
%MAKEGIF Summary of this function goes here
%   [] = makeGIF(const,scenario,RX,TX,targets,focus)
Focused_vec = focus.Focused_vec;
fig = figure("WindowState","maximized");
    images = cell(size(Focused_vec));
    for i = 1:length(focus.angle_vec)
        plotFocusedWithTargets(scenario,RX,TX,targets,Focused_vec{i},...
            strcat("Focused image with angle ",num2str(focus.angle_vec(i)),"°" ));
        frame = getframe(fig);
        images{i} = frame2im(frame);
    end
% MAKE GIF
    filename = strcat(const.experiment_name,'.gif'); % Specify the output file name
    for idx = 1:length(Focused_vec)
        [A,map] = rgb2ind(images{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.5);
        end
    end
end

