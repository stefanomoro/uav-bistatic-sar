function save_bin(filename,var)
%SAVE_BIN(filename,variable) save complex data in binary file. Filename
%without extension
%   v1 save all real and then all image
%   v2 save 1 real and 1 image, interleaved (default)

v_class = class(var);
fid = fopen(strcat(filename,'.bb'),'a');

version = 2;
if(version == 1)
    fwrite(fid,real(var),v_class);
    fwrite(fid,imag(var),v_class);
elseif(version ==2)
    A = zeros(2*numel(var),1,v_class);
    A(1:2:end-1) = real(var);
    A(2:2:end) = imag(var);
    count = fwrite(fid,A,v_class);
    if(count ~= numel(A))
        disp("ERROR Not enough space on disk")
    end
end
fclose(fid);

if(isfile(strcat(filename,"_meta.mat")))
    meta = load(strcat(filename,"_meta.mat")).meta;
    v_size = size(var);
    meta.size(2) = meta.size(2) + v_size(2);
    if(not(strcmp(meta.class,v_class)))
        error("save_bin ERROR Class of variable to append is wrong")
    end
else
    meta.size = size(var);
    meta.class = v_class;
    
end
meta.version = 2;
save(strcat(filename,"_meta"),'meta');

end

