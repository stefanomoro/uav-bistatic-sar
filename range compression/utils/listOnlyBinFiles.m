function [file_paths] = listOnlyBinFiles(folder_path)
%LISTONLYBINFILES list all .bb files in the folder path 
%   [file_paths] = listOnlyBinFiles(folder_path)
file_list = dir(folder_path);
if(length(file_list) == 2)
    disp("No file presents")
    return ;
end
file_list = file_list(3:end);
file_paths = [];
for i = 1 :length(file_list)
    %load only .bb files
    file_path = strcat(file_list(i).folder,'\',file_list(i).name);
    [~,~,ext] = fileparts(file_path) ;
    if or(file_list(i).isdir, not(strcmp(ext,'.bb')) )
        continue
    end

    file_paths(end+1).complete_path = file_path;
    file_paths(end).exp_name = file_list(i).name(1:end-3);
end
end

