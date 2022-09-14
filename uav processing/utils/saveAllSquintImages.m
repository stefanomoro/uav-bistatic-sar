function [] = saveAllSquintImages(const,scenario,RX,TX,targets,focus,mode)
%SAVEALLSQUINTIMAGES save all single focused image, Mode => 1 not eq, 2 eq
%   [] = saveAllSquintImages(const,scenario,RX,TX,targets,focus,mode)


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
%     F = 20*log10(filterHammingFocus( F_vec{i},3));
    F = 20*log10(F_vec{i});
%     F = F_vec{i};
    plotFocusedWithTargets(scenario,RX,TX,targets,F,...
        strcat("Squint angle ",num2str(focus.angle_vec(i)),"°" ));
    caxis([160 220])
    frame = getframe(fig);
    images{i} = frame2im(frame);
end
% save images

for idx = 1:length(F_vec)
    [A,map] = rgb2ind(images{idx},256);
    filename = strcat('img/',const.experiment_name,'_',num2str(idx),'.jpg');
    imwrite(A,map,filename);
end
end