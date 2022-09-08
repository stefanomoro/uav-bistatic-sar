function [C] = load_bin(filename)
%LOAD_BIN load binary file with info found in meta mat file
%   It is the simmetric of save_bin. Read complex data where samples are
%   saved : v1) first with all real and then all image part
%           v2) interleaved 1 real and 1 image

meta =load(strcat(filename,'_meta')).meta;
fid = fopen(strcat(filename,'.bb'),'r');
L = prod(meta.size);
if (not(isfield(meta,"version")))
    I = fread(fid,L,strcat(meta.class,'=>',meta.class));
    Q = fread(fid,L,strcat(meta.class,'=>',meta.class));
    C = complex(I,Q);
elseif(meta.version == 2)
    A = fread(fid,2*L,strcat(meta.class,'=>',meta.class));
    C = complex(A(1:2:end-1),A(2:2:end));    
end

C = reshape(C,meta.size);
fclose(fid);

end

