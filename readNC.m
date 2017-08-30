function vdata = readNC( fl, vname )

    vinfo = ncinfo( fl, vname );

    nDims = numel( vinfo.Size );

    if nDims > 4 || nDims < 2 
        
        error( 'Variable %s has %d dimensions.  Code can handle between 1 and 4 dimensions', vname, nDims );

    end

    if nDims == 2

        vdata = ncread( fl, vname, [1 1], [Inf Inf], [1 1] ); 

        if size( vdata )

            vdata = vdata( 1:576, 1:360 )';

            vdata( vdata > 10000000 ) = 0;

            vdata = [ vdata(:, 289:576) vdata(:, 1:288) ];

            vdata = flipud( vdata );

            return;

        else

            error( 'Variable %s is empty', vname );

        end

    else
        msg = sprintf( 'Variable %s can not be handled with %d dimensions.', ...
                                                                vname, nDims );
        disp( msg );
    end

end
