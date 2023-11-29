function [cplx_data]=RD_bin_file(name,data_len)
fid = fopen(name, 'r');
data = fread(fid, [1 2*data_len],'double') ;
cplx_data=data(1:data_len)+1j*data(data_len+1:data_len*2);
fclose(fid);
end
