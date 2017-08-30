%function res = loadNC( fl, var1, var2, var3, var4, var5, var6, var7 )
function res = loadNC( fl, varargin )

    if isempty( varargin )

        error( 'Function requires at least 2 arguments.' )

    end

    ncid = netcdf.open( fl, 'NC_NOWRITE' );

    [ ndims, nvars, nglobalatts, unlimdimID ] = netcdf.inq( ncid );

    %for i = 1:nargin-1
    for i = 1:length(varargin)

        %vname   = eval( ['var', num2str(i)] );
        vname   = varargin{ i };

        vid     = netcdf.inqVarID( ncid, vname );

        xvarname{ i } = vname;
        
        try 
            scale_factor = netcdf.getAtt( ncid, vid, 'scale_factor' );
        catch me
            scale_factor = 1.0;
        end

        try
            add_offset = netcdf.getAtt( ncid, vid, 'add_offset' );
        catch me
            add_offset = 0.0;
        end

        var         = netcdf.getVar( ncid, vid );

        disp( var );
%{
        xvar{ i }   = double(double(var) * scale_factor + add_offset);

        %assignin( 'caller', xvarname{i}, xvar{i} );
        assignin( 'caller', xvarname{i}, xvar{i} );
%}
        clear var; 
    end

    res = 0;
    %res = xvarname;

    %netcdf.close( ncid );

end
