% Convert data from netCDF to single files for CHRS CONNECT code
% Phu Nguyen 03112016
clear all

ncid=netcdf.open('./Data/IVT_201212.nc','nowrite');
    
    % CDF Data information
    [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
    
    for i = 0:numvars-1
        [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
        flag = 0;
        for j = 0:numatts - 1
            attname1 = netcdf.inqAttName(ncid,i,j);
            attname2 = netcdf.getAtt(ncid,i,attname1);
            disp([attname1 ':  ' num2str(attname2)])
            if strmatch('add_offset',attname1)
                offset = attname2;
            end
            if strmatch('scale_factor',attname1)
                scale = attname2;
                flag = 1;
            end
        end
        
        if flag
            eval([varname '= double(double(netcdf.getVar(ncid,i))*scale + offset);'])
        else
            eval([varname '= double(netcdf.getVar(ncid,i));'])
        end
    end
    netcdf.close(ncid);
    
    n=size(IVT,3);
    for i=1:n
        fn_bin=['./Data/IVT' sprintf('%03d',i) '.bin'];
        data=IVT(:,:,i)';
    
        write_bn(fn_bin,data,'float32');
        system(['gzip ' fn_bin]);
    end