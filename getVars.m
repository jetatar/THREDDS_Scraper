function getVars( fl, varargin )

    % Directory to store the newly created .nc files with variables of interest.
    % Directory name must end with '/'
    storeIn = '/home/jtatar/Work/CHRS/CHRS_Connect/TestData/merra/Vars/';
    logfl   = 'processedDBFiles.txt';

    if isempty( varargin )

        error( 'Function requires at least 2 arguments.' )

    end

    newfl = strsplit( fl, '/' ); 
    newfl = char( newfl(end) );
    flname= newfl;
                   
    dpos = find( fl == '/', 1, 'last');
    addr = fl(1:dpos);

    if ~isempty( storeIn )
        newfl = strcat( storeIn, newfl );
    end

    newfl = newfl( 1:end-4 );

    for i = 1:numel( varargin )

        vinfo       = ncinfo( fl, varargin{i} );

        dimNames    = { vinfo.Dimensions.Name };
        ndims       = length( dimNames );
        dimLen      = { vinfo.Dimensions.Length };
        dType       = { vinfo.Datatype };
        dType       = ncVarType( dType );

        if ndims == 2

            if strmatch( 'lon', dimNames{1} ) && strmatch( 'lat', dimNames{2} )
                % Do nothing, write variable to disk.

                vdata = ncread( fl, varargin{i} );

                infstr = ['_', varargin{i}, '_VAR.nc'];

                savefl = strcat( newfl, infstr );

                ncid        = netcdf.create( savefl, 'NC_WRITE' );
                dimidrow    = netcdf.defDim( ncid, 'lon', dimLen{1} );
                dimidcol    = netcdf.defDim( ncid, 'lat', dimLen{2} );
                varid       = netcdf.defVar( ncid, varargin{i}, ...
                                dType, [dimidrow, dimidcol] );
                netcdf.endDef( ncid );
                netcdf.putVar( ncid, varid, vdata );
                netcdf.close( ncid );
                disp( 'Saved file: ' );
                disp( savefl );

            else
                error( 'Unknown dimension names.' );
            end

        elseif ndims == 3

            if strmatch( 'lon', dimNames{1} ) && strmatch( 'lat', dimNames{2} ) ...
                                            && strmatch( 'lev', dimNames{3} )

                % Get one of the levels and copy data to disk.
                vdata = ncread( fl, varargin{i} );

                level = 1;

                tfdata = vdata( :, :, level );

                infstr = ['_', varargin{i}, num2str(level, '_lev%d'), ...
                                            '.nc'];

                savefl = strcat( newfl, infstr );


                ncid        = netcdf.create( savefl, 'NC_WRITE' );
                dimidrow    = netcdf.defDim( ncid, 'lon', dimLen{1} );
                dimidcol    = netcdf.defDim( ncid, 'lat', dimLen{2} );
                varid       = netcdf.defVar( ncid, varargin{i}, ...
                                dType, [dimidrow, dimidcol] );
                netcdf.endDef( ncid );
                netcdf.putVar( ncid, varid, tfdata );
                netcdf.close( ncid );
                disp( 'Saved file: ' );
                disp( savefl );

            else
                error( 'Unknown dimension names.' );
            end

        elseif ndims == 4

            if strmatch( 'lon', dimNames{1} ) && strmatch( 'lat', dimNames{2} ) ...
                && strmatch( 'lev', dimNames{3} ) && strmatch( 'time', dimNames{4} )

                if dimLen{4} == 8
                    vdata = ncread( fl, varargin{i} );

                    level = 1;


                    for t = 1:8
                        tfdata = vdata( :, :, level, t );

                        hr = (t - 1) * 3;

                        infstr = ['_', varargin{i}, num2str(level, '_lev%d_'), ...
                                                    num2str(hr, '%02d'),'HR.nc'];

                        savefl = strcat( newfl, infstr );


                        ncid        = netcdf.create( savefl, 'NC_WRITE' );
                        dimidrow    = netcdf.defDim( ncid, 'lon', dimLen{1} );
                        dimidcol    = netcdf.defDim( ncid, 'lat', dimLen{2} );
                        varid       = netcdf.defVar( ncid, varargin{i}, ...
                                        dType, [dimidrow, dimidcol] );
                        netcdf.endDef( ncid );
                        netcdf.putVar( ncid, varid, tfdata );
                        netcdf.close( ncid );
                        disp( 'Saved file: ' );
                        disp( savefl );
                    end

                else
                    error( 'Time dimension does not match expected length.' );
                end


            else
                error( 'Unknown dimension names.' );
            end

        else
            error( 'Cannot deal with %d dimensions.', ndims );

        end

    end

    fid_proc = fopen( logfl, 'a' );
    fprintf( fid_proc, '%s%s\n', addr, flname );
    fclose( fid_proc );

    fprintf ( 'File %s.nc4 has been split into time frames.\n', newfl );

end
