function [] = makeGIF(const,scenario,RX,TX,targets,focus,mode)
%MAKEGIF Make gif with multiple focused image, Mode => 1 not eq, 2 eq
%   [] = makeGIF(const,scenario,RX,TX,targets,focus,mode)


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
images = cell(size(F_vec));

        
for i = 1:length(focus.angle_vec)
    plotFocusedWithTargets(scenario,RX,TX,targets,F_vec{i},...
        strcat("Focused image with angle ",num2str(focus.angle_vec(i)),"°" ));
    frame = getframe(fig);
    images{i} = frame2im(frame);
end
% MAKE GIF
    filename = strcat(const.experiment_name,'.gif'); % Specify the output file name
    for idx = 1:length(F_vec)
        [A,map] = rgb2ind(images{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.5);
        end
    end
end

