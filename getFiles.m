%
%   Function should return list of urls to all files found by parsing the THREDDS 
%   database.
%

function urls = getFiles( )

    sourceList  = formURL( datetime(2011, 1, 1, 'Format', 'yyyy/MM'), ...
                           datetime(2011, 1, 1, 'Format', 'yyyy/MM') );
    %sourceList  = cellstr( sourceList );

    urlPath     = {};

    for s = 1 : length(sourceList)

        try 
            xmlDoc      = xmlread( sourceList{s} );
            msg         = sprintf( '-> Parsing %s.', sourceList{s} );
            disp( msg );

            allListitems= xmlDoc.getElementsByTagName( 'dataset' );
            len         = allListitems.getLength() - 1;

            % Loop over attribute list
            for i = 0 : len

                uPath   = allListitems.item(i).getAttribute( 'urlPath' );
                uPath   = char( uPath );    % Java string to Matlab char array

                if( ~isempty(uPath) )
                    uPath = ['http://67.58.48.50:8080/thredds/dodsC/' uPath];
                    uPath = cellstr(uPath);

                    urlPath(end + 1) = uPath;
                end

            end

            %if( ~isempty(urlPath) )
            %    msg = sprintf( '-> Found %d files', length(urlPath) );
            %    disp( msg );
            %end

        catch
            msg     = sprintf( '!> Source %s not found, skipping.', sourceList{s} );
            disp( msg );
            % TODO: Skip and go to next file.
        end

    end % for source list loop

    fl = fopen( 'datafiles.txt', 'w' );

    for row = 1:length(urlPath)
        fprintf( fl, '%s\n', urlPath{row} );
    end
    
    fclose( fl );

    urls = urlPath;

    msg  = sprintf( '-> Total number of individual files found: %d', ...
                                                            length(urlPath) );
    disp( msg );

end


function sourceList = formURL( s_date, e_date )

    i = { 5 };
    j = { 'M2I3NPASM.5.12.4/' };  % NOTE: add slash at end

    dt    = months( datestr(s_date), datestr(e_date) );

    sourceList = {};

    for a = i
        for b = j
            for c = 0 : dt  % +1 month iter

                date    = datetime( s_date ) + calmonths( c );
                i_addr  = [ 'http://67.58.48.50:8080/thredds/catalog/testAll/' ...
                            'goldsmr' num2str(a{:}) ...
                            '.gesdisc.eosdis.nasa.gov/data/MERRA2/' b{:} ... 
                            datestr(date, 'yyyy/mm') '/catalog.xml' ];
                i_addr  = cellstr( i_addr );

                sourceList(end + 1) = i_addr;

            end
        end
    end
end

%{
function getFilesFromURL( url, s_date, e_date )

    % Form search string.
    addr1 = 'http://67.58.54.105:8080/thredds/catalog';
    addr2 = '/testAll/goldsmr{4}.gesdisc.eosdis.nasa.gov/';
    addr3 = '/data/MERRA2/M2T1NXCHM.5.12.4/';
    addr4 = '$yyyy/$mm/catalog.html'

    serv    = {4, 5};
    chans   = {'M2T1NXAER.5.12.4', 'M2T1NXCSP.5.12.4'};

    if isnan( e_date )
        e_date = Date();
    else
        e_date = edate.getYYYY() edate.getMM()
    end

    % Iterate over choices
    for i = 1 : length( serv )
        for j = 1 : length( chans )
            iterrate over date           
                read xml file
                extract data files
    addr = a

    % Parse xml file and generate urls.

    %   Get list of files to process.  Write it.  Read it, compared it
    %   against a 'already processed' files list and only continue where left off.

end
%}
