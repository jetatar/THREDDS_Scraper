function nctype = ncVarType( name )

    if isequal( name, 'double' )

        nctype = 'NC_DOUBLE';

    elseif isequal( name, 'single' )

        nctype = 'NC_FLOAT';

    elseif isequal( name, 'int64' )

        nctype = 'NC_INT64';

    elseif isequal( name, 'uint64' )

        nctype = 'NC_UINT64';

    elseif isequal( name, 'int32' )

        nctype = 'NC_INT';

    elseif isequal( name, 'uint32' )

        nctype = 'NC_UINT';

    elseif isequal( name, 'int16' )

        nctype = 'NC_SHORT';

    elseif isequal( name, 'uint16' )

        nctype = 'NC_USHORT';

    elseif isequal( name, 'int8' )

        nctype = 'NC_BYTE';

    elseif isequal( name, 'uint8' )

        nctype = 'NC_UBYTE';

    elseif isequal( name, 'char' )

        nctype = 'NC_CHAR';

    else

        error( 'Unknown data type!.' ); 

end
