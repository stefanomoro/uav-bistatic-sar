function save_bin(filename,var)
%SAVE_BIN(filename,variable) save complex data in binary file. Filename
%without extension
%   Detailed explanation goes here
v_class = class(var);
fid = fopen(strcat(filename,'.bb'),'w');
fwrite(fid,real(var),v_class);
fwrite(fid,imag(var),v_class);
fclose(fid);

meta.size = size(var);
meta.class = v_class;
save(strcat(filename,"_meta"),'meta');
end

