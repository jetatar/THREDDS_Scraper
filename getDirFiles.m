function fls = getDirFiles( ndir, varname, bFull )

    allfls = { };

    flstr = [ '*_', varname, '_*.nc' ];
    flstr = strcat( ndir, flstr );

    fs  = dir( flstr );

    if int8(bFull)
        fs  = strcat( ndir, {fs.name} );
    else
        fs  = {fs.name};
    end

    fls = fs;

end
