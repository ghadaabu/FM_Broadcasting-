function []=WR_bin_file(name,data)
fid = fopen(name, 'w');
r=real(data);
im=imag(data);
z=zeros(1,length(data)*2);
z(1:length(data))=r;
z(length(data)+1:length(data)*2)=im;
fwrite(fid, z, 'double');
fclose(fid);
end

