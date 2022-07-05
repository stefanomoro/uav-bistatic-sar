function [file_name] = getResultsFileName(focus)
%SAVEFOCUSED get first available filename for saving results
%   [file_name] = getSaveFileName(focus)

i = 0;
last_ang = focus.angle_vec(end);
fname = strcat("focused_",num2str(last_ang),"_",num2str(i));
while isfile(strcat("focused_result/",fname))
    i = i +1;
end
file_name = strcat("focused_result/",fname);
end

