function [C] = load_bin(filename)
%LOAD_BIN load binary file with info found in meta mat file
%   It is the simmetric of save_bin. Read complex data where samples are
%   saved first with all real and then all image part

meta =load(strcat(filename,'_meta')).meta;
fid = fopen(strcat(filename,'.bb'),'r');
L = prod(meta.size);
C = complex(fread(fid,L,strcat(meta.class,'=>',meta.class)),...
    fread(fid,L,strcat(meta.class,'=>',meta.class)));
C = reshape(C,meta.size);
fclose(fid);

end

