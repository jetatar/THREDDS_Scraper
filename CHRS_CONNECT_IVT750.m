% Code for segmenting extreme events from satellite precipiation data
% ***** by Phu Nguyen ***** 11/22/2011
% updated 09/05/2015
% MERRA IVT data for Scott Sellars
clear all
tic

%% INPUT DATA
% Threshold for segmenting objects 
thres_precip = 750; % 
% Threshold for choosing extreme events
thres_maxprecip = 1; % mm/d
% Number of pixels (area) for extreme events
n_pixels_thres = 1;
% Number of connectivities (6, 26)
conn=26;
% Threshold duration for segmenting extreme events (hours, days)
thres_duration = 1; % hours
% Number of images in 1 subdivision (split the dataset into small chunks for limited RAM)
n_image = 500;
% Number of buffering images (when merging subdivisions we need to process the buffering areas)
n_buffer =100; % assuming 100days overlapping
% Number of rows in the image
n_row=360;
% Number of columns in the image
n_col=576;
% Size of image
n_size=n_row*n_col;
% Number of extra columns in longitude buffering areas
n_colx=180; %90 degree
% Size of image with extra columns
n_sizex = n_row*(n_col+n_colx);
% Data resolution (degree)
resx=0.625;
resy=0.5;
% Data directory
Data_path='/home/jtatar/Work/CHRS/CHRS_Connect/TestData/';
% Directory where files containing individual timeframes will be stored.
TFFilesDir = '/home/jtatar/Work/CHRS/CHRS_Connect/TestData/merra/Vars/';
% Result path:
Results_path='/home/jtatar/CHRS/data/OrignalCodeOutput/';

%% PROCESSING DATA
% Creating directories for storing results
if (exist([Results_path 'Object_Table'],'dir')<1)
    mkdir([Results_path 'Object_Table'])
end

if (exist([Results_path 'Boundary'],'dir')<1)
    mkdir([Results_path 'Boundary'])
end

if (exist([Results_path 'Track','dir'])<1)
    mkdir([Results_path 'Track'])
end

if (exist([Results_path 'Track_DWR'],'dir')<1)
    mkdir([Results_path 'Track_DWR'])
end


ivtVar = containers.Map( ...
        {'name', 'lev', 't_start', 't_stop', 'valLow', 'valHi'}, ...
        {'EPV',      1,       NaN,      NaN,      750,  NaN} );

TVar = containers.Map( ...
        {'name', 'lev', 't_start', 't_stop', 'valLow', 'valHi'}, ...
        {'T',      1,       NaN,      NaN,      750,  NaN} );

VAR  = { ivtVar, TVar };

runPeriod = 60; % how often to run the code.


%
% The main infinite loop.
%
while true

    % Data files
    allDBFiles = getFiles( );

    varnames = { };

    for v = 1:length( VAR )
        el = VAR{v};
        varnames = [ varnames; {el('name')} ];
    end

    % Name of file with processed file names.
    processedDBFiles = 'processedDBFiles.txt';
    procDBFiles      = { };

    s = dir( processedDBFiles );

    if exist( processedDBFiles, 'file' ) == 2

        if s.bytes ~= 0

            %fl          = fopen( processedDBFiles, 'r' );
            %procDBFiles = fscanf( fl, '%s' );
            %fclose( fl );
            procDBFiles = importdata( processedDBFiles );

        end

    end

    % A list of addresses of THREDDS files that have not been broken down to
    % individual time frames.

    flsToBreak = setdiff( allDBFiles, procDBFiles );

    if ~isempty( flsToBreak )

        for ff = 1:length(flsToBreak)

            fprintf( 'Breaking up %s\n', flsToBreak{ff} );

            getVars( flsToBreak{ff}, varnames{:} );

        end

    else

        disp( 'No new DB files to reduce to individual time frames.' );

    end

    % Get all available single time frame files and run CONNECT algorithm only on
    % the ones that have not been processed yet.
    %allTFFls = { };

    for vs = 1:length(varnames)

        allTFFls    = getDirFiles( TFFilesDir, varnames{vs}, 0 );% 0 = file names only
        allTFProcFls= getDirFiles( Results_path, varnames{vs}, 0 );

        files = setdiff( allTFFls, allTFProcFls );
        files = strcat( TFFilesDir, files );

        % 'files' is a cell array with rows corresponding to the variable number
        % while the columns are the stored filenames for that variable
        % If any file needs to be processed, 'files' shouldn't be empty.
        if ~isempty( files )

            % number of files in the dataset
            n = length( files );

%n = length(files)

% Number of subdivisions
n_sub = ceil(n/n_image);

% Initializing variables
id=1;
id1=1;
id2=1;
maxxx=(n-1)*n_sizex;
vox1=[];
vox_dat1=[];
tx1=[];

% Main loop
t0 = datenum(1980,01,01,00,00,00); %(Y,M,D,H,MN,SS)

delta_t=3; % 3hourly

if n==0
    disp('No data found');
    
    return;
end


%% Main Loop
for t=1:n_sub
    t
    dat0=[];
    % Reading and preliminarily processing data
    str = sprintf('Processing sub %d of %d', t, n_sub);
    disp(str);
    str = sprintf('Reading data');
tic    
    disp(str);
    
    if n_sub==1
        dat0=zeros(n_row,n_col,n);
        for i=1:n
            %filename = files(i,1).name;
            filename = files{i};

            outstr = sprintf( 'Processing file %s', filename );
            disp( outstr );

            %dat0(:,:,i) = readNC([Data_path filename]);
            dat0(:,:,i) = readNC(filename, varnames{vs});

            e = ['cp -rs ' filename ' ' Results_path ];
            system( e );

            %disp( dat0 );
            %disp(size(dat0));
            %disp(class(dat0));
        end
        
    elseif t<n_sub
        dat0=zeros(n_row,n_col,n_image);
        for i=1:n_image
            filename = files{id};
            dat0(:,:,i) = readNC( filename, varnames{vs} );
            e = ['cp -rs ' filename ' ' Results_path ];
            system( e );
            id=id+1;
        end
        
        if t==n_sub-1
            for i=1:min(n_buffer,n-(n_sub-1)*n_image)
                filename = files{i+id-1};
                dat0(:,:,n_image+i) = readNC(filename, varnames{vs});
                e = ['cp -rs ' filename ' ' Results_path ];
                system( e );
            end
        else
            for i=1:n_buffer
                filename = files{i+id-1};
                dat0(:,:,n_image+i) = readNC(filename, varnames{vs});
                e = ['cp -rs ' filename ' ' Results_path ];
                system( e );
            end
            
        end
    elseif (t==n_sub) && (mod(n,n_image)>=1)
        dat0=zeros(n_row,n_col,mod(n,n_image));
        for i=1:mod(n,n_image)
            filename = files{id};
            dat0(:,:,i) = readNC(filename, varnames{vs});
            e = ['cp -rs ' filename ' ' Results_path ];
            system( e );
            id=id+1;
        end
        
    elseif (t==n_sub) && (mod(n,n_image)==0)
        dat0=zeros(n_row,n_col,n_image);
        for i=1:n_image
            filename = files{id};
            dat0(:,:,i) = readNC(filename, varnames{vs});
            e = ['cp -rs ' filename ' ' Results_path ];
            system( e );
            id=id+1;
        end
    end
toc    
    % Using BWCONNCOMP function -> list of objects in structure format - tmp (vox)
    str = sprintf('Using BWCONNCOMP function to extract objects');
    disp(str);
tic    
    %disp( dat0 );
    dat=[dat0 dat0(:,1:n_colx,:)];
    %disp( dat );
    %vox=vox_fn(dat,conn,thres_precip);
    vox=vox_fn(dat,conn,VAR{vs});
    
    disp( vox );
toc    


    %% Matching objects in time buffering areas
    disp('Matching objects in buffering areas');
tic
    check=1;
    maxx0=(n_image+1)*n_sizex;
    nx=n_buffer;
    idx=id-1+n_buffer;
    
    if t<n_sub
        while check==1
            
            check=0;
            maxx1=(n_image+nx-1)*n_sizex;
            minn=cellfun(@min,vox);
            maxx=cellfun(@max,vox);
            
            test=minn(maxx>maxx1);
            if any(test<=maxx0)
                check=1;
            end
            
            if check==1
                if (idx+n_buffer)>=n
                    datx=zeros(n_row,n_col,n-idx);
                    for j=1:(n-idx)
                        filename = files{j+idx};
                        datx(:,:,j) = readNC(filename, varnames{vs});
                        e = ['cp -rs ' filename ' ' Results_path ];
                        system( e );
                    end
                    check=0;
                    
                else
                    datx=zeros(n_row,n_col,n_buffer);
                    for j=1:n_buffer
                        filename = files{j+idx};
                        datx(:,:,j) = readNC(filename, varnames{vs});
                        e = ['cp -rs ' filename ' ' Results_path ];
                        system( e );
                    end
                    nx=nx+n_buffer;
                    idx=idx+n_buffer;
                end
                dat0=cat(3,dat0,datx);
                dat=[dat0 dat0(:,1:n_colx,:)];
                vox=vox_fn(dat,conn,thres_precip);
            end
            
        end
        
    end
    %clear dat0
toc    
    % Removing objects having pixels on beginning section
    if t>1
        minn=cellfun(@min,vox);
        vox=vox(minn>n_sizex);
    end
    
    % Removing objects having  pixels on ending section
    if t<n_sub
        minn=cellfun(@min,vox);
        vox(minn>maxx0)=[];
    end
    if isempty(vox)
        return;
    end
    
    
%% Removing redundant objects in longitude buffering areas
    
    disp('Removing redundant objects in longitude buffering areas');
    
    temp=cellfun(@(x) floor(mod(x,n_sizex)/n_row)+1,vox,'uni',false);
    
    vox_lon=cellfun(@(x) -resx+x*360/n_col,temp,'uni',false);
    
    minn=cellfun(@min,vox_lon);
    maxx=cellfun(@max,vox_lon);
    temp=union(intersect(find(minn>0),find(maxx<=360-resx)),intersect(find(minn<=360),find(maxx>=360)));
    vox=(vox(temp));
    
    if isempty(vox)
        return;
    end
    
    % Removing small objects using threshold
    disp('Removing small objects based on duration threshold');
    
    minn=cellfun(@min,vox);
    maxx=cellfun(@max,vox);
    
    duration=ceil((maxx - minn)/n_sizex*delta_t);
    vox=(vox(duration>=thres_duration));

    disp( vox );    

    if isempty(vox)
        return;
    end
    
    %% Calculating incomplete objects
    disp('Calculating incomplete objects');
    
    maxx=cellfun(@(x) max(x)+(t-1)*n_image*n_sizex,vox);
    temp=find(maxx>maxxx);
    
    if ~isempty(temp)
        
        vox0=[{vox{temp}}];
        tx0=ones(length(vox0),1)*t;
        
        vox(temp)=[];
        vox_dat0=cellfun(@(x) dat(x),vox0,'uni',false);
        vox1=[vox1 vox0];
        tx1=[tx1;tx0];
        vox_dat1=[vox_dat1 vox_dat0];
    end
    
    %% Calculating Object Tables
    disp('Calculating object tables');
    if ~isempty(vox)
        
    Z=sum(cellfun(@numel,vox));
    Table=zeros(Z,4);
    Table_time=zeros(Z,6);
    z=sum(cellfun(@(x) (ceil(max(x)/n_sizex)-ceil(min(x)/n_sizex)),vox))+length(vox);
    Track=zeros(z,9);
    
    % object ID
    obj_id=cellfun(@(x) ones(length(x),1),vox,'uni',false);
    for i=1:length(obj_id)
        obj_id{i}=(i+id1-1)*obj_id{i};
    end
    id1=id1+length(vox);
    
    % latitude
    temp=cellfun(@(x) mod(mod(x,n_sizex),n_row),vox,'uni',false);
    temp=cellfun(@(x) cellreplace_func(x,0,n_row),temp,'uni',false);
    lat=cellfun(@(x) 90 - x/n_row*180,temp,'uni',false);
    
    % longitude
    temp=cellfun(@(x) ceil( mod(x,n_sizex)/n_row),vox,'uni',false);
    temp=cellfun(@(x) cellreplace_func(x,0,n_col+n_colx),temp,'uni',false);
    lon=cellfun(@(x) -resx + x*360/n_col,temp,'uni',false);
    
    % intensity
    intensity=cellfun(@(x) dat(x),vox,'uni',false);
    
    % date and time
    time= cellfun(@(x) t0 + ( ceil(x/n_sizex) -1 + (t-1)*n_image )*delta_t/24,vox,'uni',false);
    time= cellfun(@datevec,time,'uni',false);
    
    % Calculating Tracking table
    i1=1;
    i2=1;
    for i=1:length(vox)
        L=length(vox{i});
        temp=ceil((vox{i})/n_sizex);
        k1=min(temp);
        k2=max(temp);
        sub_track=zeros(k2-k1+1,9);
        timei=time{i};
        y=timei(:,1);
        m=timei(:,2);
        d=timei(:,3);
        H=timei(:,4);
        M=timei(:,5);
        S=timei(:,6);
        ID=obj_id{i};
        Lat=lat{i};
        Lon=lon{i};
        
        for j=k1:k2
            
            sub_track(j-k1+1,1)=mean( ID(temp==j));
            sub_track(j-k1+1,2)=mean( Lat(temp==j));
            sub_track(j-k1+1,3)=mean( Lon(temp==j));
            sub_track(j-k1+1,4)=mean( y(temp==j));
            sub_track(j-k1+1,5)=mean( m(temp==j));
            sub_track(j-k1+1,6)=mean( d(temp==j));
            sub_track(j-k1+1,7)=mean( H(temp==j));
            sub_track(j-k1+1,8)=mean( M(temp==j));
            sub_track(j-k1+1,9)=mean( S(temp==j));
            
        end
        
        Track(i2:i2+k2-k1,:)=sub_track;
        i1=i1+L;
        i2=i2+k2-k1+1;
        
    end
    
    
    %% Exporting Voxel Tables
    obj_id=[cat(1,obj_id{:})];
    lat=[cat(1,lat{:})];
    lon=[cat(1,lon{:})];
    time=[cat(1,time{:})];
    intensity=[cat(1,intensity{:})];
    
    disp('Exporting Tables');
    fn=[Results_path sprintf('Object_Table/Object_Table%d.asc',t)];
    fid=fopen(fn,'w');
    fprintf(fid,'%i,%8.4f,%9.4f, %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d, %8.4f\n',[ obj_id'; lat';lon';time(:,1)';time(:,2)';time(:,3)';time(:,4)';time(:,5)';time(:,6)';intensity']);
    fclose(fid);
    
    %% Exporting Tracking Tables
    
    formatSpec1='(%i,ST_GeomFromText(''LINESTRING(';
    formatSpec2='%7.3f %7.3f,';
    formatSpec3='%7.3f %7.3f)'',4326)),\n';
    formatSpec4='%7.3f %7.3f)'',4326))';
    
    hmin=min(Track(:,1));
    hmax=max(Track(:,1));
    h=hmax-hmin+1;
    
    fn=[Results_path sprintf('Track/Track%d.sql',t)];
    fid_track=fopen(fn,'w');
    
    for i=1:h
        idx2=Track(:,1)==(i+hmin-1);
        Tracki=Track(idx2,:);
        
        fprintf(fid_track,formatSpec1,Tracki(1,1));
        [l k]=size(Tracki);
        if i<h
            if l>1
                for j=1:l-1
                    fprintf(fid_track,formatSpec2,Tracki(j,3),Tracki(j,2));
                end
                fprintf(fid_track,formatSpec3,Tracki(j+1,3),Tracki(j+1,2));
            else
                fprintf(fid_track,formatSpec2,Tracki(1,3),Tracki(1,2));
                fprintf(fid_track,formatSpec3,Tracki(1,3),Tracki(1,2));
            end
        else
            if l>1
                for j=1:l-1
                    fprintf(fid_track,formatSpec2,Tracki(j,3),Tracki(j,2));
                end
                fprintf(fid_track,formatSpec4,Tracki(j+1,3),Tracki(j+1,2));
            else
                fprintf(fid_track,formatSpec2,Tracki(1,3),Tracki(1,2));
                fprintf(fid_track,formatSpec4,Tracki(1,3),Tracki(1,2));
            end
            
        end
    end
    fclose(fid_track);
    
    %% Writing tracking for DWR
    disp('Exporting tracking for DWR');
    fn=[Results_path sprintf('Track_DWR/Track_DWR%d.asc',t)];
    fid_track_DWR=fopen(fn,'w');
    
    fprintf(fid_track_DWR,'%i,%8.4f,%9.4f, %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d\n',[Track(:,1)'; Track(:,2)';Track(:,3)';Track(:,4)';Track(:,5)';Track(:,6)';Track(:,7)';Track(:,8)';Track(:,9)']);      
    fclose(fid_track_DWR);
    %% Calculating Object Boundaries
    disp('Calculating object boundary');
    
    Boundary={};
    
    for i=1:length(vox)
        
        dat_vox = zeros(n_row,(n_col + n_colx));
        
        temp=ceil( mod(vox{i},n_sizex)/n_row);
        temp(temp==0)=n_col+n_colx;
        vox_x= temp;
        
        vox_y= mod( mod(vox{i},n_sizex),n_row);
        vox_y(vox_y==0)=n_row;
        
        for j=1:length(vox_x)
            dat_vox(vox_y(j),vox_x(j))=1;
        end
        
        [B,L] = bwboundaries(dat_vox,'noholes');
        m0=1;n0=1;
        for m=1:length(B)
            if length(B{m})>n0
                n0=length(B{m});
                m0=m;
            end
        end
        
        K = B{m0};
        n_points = length(B{m0});
        temp = zeros(n_points,3);
        
        temp(:,1)=id2;
        temp(:,2) = 90 - K(:,1)/n_row*180;                      % Lat
        temp(:,3) = -resx + K(:,2)*360/n_col;                         % Lon
        
        if i==1
            Boundary=temp;
        else
            Boundary=[Boundary; temp];
        end
        
        id2=id2+1;
    end
    
    
    %% Exporting object boundary tables
    disp('Exporting object boundary tables');
    
    formatSpec1='(%i,ST_MakePolygon(ST_GeomFromText(''LINESTRING(';
    formatSpec2='%7.3f %7.3f,';
    formatSpec3='%7.3f %7.3f';
    formatSpec4=')'',4326)),''%4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d'',';
    formatSpec5='''%4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d''),\n';
    formatSpec6='''%4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d'')';
    
    hmin=min(Track(:,1));
    hmax=max(Track(:,1));
    h=hmax-hmin+1;
    
    fn=[Results_path sprintf('Boundary/Boundary%d.sql',t)];
    
    fid_boundary=fopen(fn,'w');
    
    
    for i=1:h
        idx2=Track(:,1)==(i+hmin-1);
        Tracki=Track(idx2,:);
        
        idx3=Boundary(:,1)==(i+hmin-1);
        Boundaryi=Boundary(idx3,:);
        
        fprintf(fid_boundary,formatSpec1,Boundaryi(1,1));
        [l k]=size(Boundaryi);
        [g k]=size(Tracki);
        if i<h
            if l>1
                for j=1:l-1
                    fprintf(fid_boundary,formatSpec2,Boundaryi(j,3),Boundaryi(j,2));
                end
                fprintf(fid_boundary,formatSpec3,Boundaryi(l,3),Boundaryi(l,2));
                fprintf(fid_boundary,formatSpec4,Tracki(1,4),Tracki(1,5),Tracki(1,6),Tracki(1,7),Tracki(1,8),Tracki(1,9));
                fprintf(fid_boundary,formatSpec5,Tracki(g,4),Tracki(g,5),Tracki(g,6),Tracki(g,7),Tracki(g,8),Tracki(g,9));
            else
                fprintf(fid_boundary,formatSpec4,Tracki(1,4),Tracki(1,5),Tracki(1,6),Tracki(1,7),Tracki(1,8),Tracki(1,9));
                fprintf(fid_boundary,formatSpec5,Tracki(g,4),Tracki(g,5),Tracki(g,6),Tracki(g,7),Tracki(g,8),Tracki(g,9));
            end
        else
            if l>1
                for j=1:l-1
                    fprintf(fid_boundary,formatSpec2,Boundaryi(j,3),Boundaryi(j,2));
                end
                fprintf(fid_boundary,formatSpec3,Boundaryi(l,3),Boundaryi(l,2));
                fprintf(fid_boundary,formatSpec4,Tracki(1,4),Tracki(1,5),Tracki(1,6),Tracki(1,7),Tracki(1,8),Tracki(1,9));
                fprintf(fid_boundary,formatSpec6,Tracki(g,4),Tracki(g,5),Tracki(g,6),Tracki(g,7),Tracki(g,8),Tracki(g,9));
            else
                fprintf(fid_boundary,formatSpec4,Tracki(1,4),Tracki(1,5),Tracki(1,6),Tracki(1,7),Tracki(1,8),Tracki(1,9));
                fprintf(fid_boundary,formatSpec6,Tracki(g,4),Tracki(g,5),Tracki(g,6),Tracki(g,7),Tracki(g,8),Tracki(g,9));
            end
            
        end
    end

    fclose(fid_boundary);
    end
end



%% Calculating Incomplete Objects
disp('Calculating incomplete object tables');
if ~isempty(vox1)
    
    Z=sum(cellfun(@numel,vox1));
    Table=zeros(Z,4);
    Table_time=zeros(Z,6);
    z=sum(cellfun(@(x) (ceil(max(x)/n_sizex)-ceil(min(x)/n_sizex)),vox1))+length(vox1);
    Track=zeros(z,9);
    
    % object ID
    obj_id=cellfun(@(x) ones(length(x),1),vox1,'uni',false);
    for i=1:length(obj_id)
        obj_id{i}=(i+id1-1)*obj_id{i};
    end
    id1=id1+length(vox1);
    
    % latitude
    temp=cellfun(@(x) mod(mod(x,n_sizex),n_row),vox1,'uni',false);
    temp=cellfun(@(x) cellreplace_func(x,0,n_row),temp,'uni',false);
    lat=cellfun(@(x) 90 - x/n_row*180,temp,'uni',false);
    
    % longitude
    temp=cellfun(@(x) ceil( mod(x,n_sizex)/n_row),vox1,'uni',false);
    temp=cellfun(@(x) cellreplace_func(x,0,n_col+n_colx),temp,'uni',false);
    lon=cellfun(@(x) -resx + x*360/n_col,temp,'uni',false);
    
    % intensity
    intensity=vox_dat1;
    
    % date and time
    time= cellfun(@(x) t0 + ( ceil(x/n_sizex) -1 + (t-1)*n_image )*delta_t/24,vox1,'uni',false);
    time= cellfun(@datevec,time,'uni',false);
    
    % Calculating Tracking table
    i1=1;
    i2=1;
    for i=1:length(vox1)
        L=length(vox1{i});
        temp=ceil((vox1{i})/n_sizex);
        k1=min(temp);
        k2=max(temp);
        sub_track=zeros(k2-k1+1,9);
        timei=time{i};
        y=timei(:,1);
        m=timei(:,2);
        d=timei(:,3);
        H=timei(:,4);
        M=timei(:,5);
        S=timei(:,6);
        ID=obj_id{i};
        Lat=lat{i};
        Lon=lon{i};
        
        for j=k1:k2

            %sub_track(j-k1+1,4:9)=mean( timei(temp==j,:));
            
            sub_track(j-k1+1,1)=mean( ID(temp==j));
            sub_track(j-k1+1,2)=mean( Lat(temp==j));
            sub_track(j-k1+1,3)=mean( Lon(temp==j));
            sub_track(j-k1+1,4)=mean( y(temp==j));
            sub_track(j-k1+1,5)=mean( m(temp==j));
            sub_track(j-k1+1,6)=mean( d(temp==j));
            sub_track(j-k1+1,7)=mean( H(temp==j));
            sub_track(j-k1+1,8)=mean( M(temp==j));
            sub_track(j-k1+1,9)=mean( S(temp==j));
            
        end
        
        Track(i2:i2+k2-k1,:)=sub_track;
        i1=i1+L;
        i2=i2+k2-k1+1;
        
    end
    
    
    %% Exporting voxel Tables
    disp('Exporting Tables');
    obj_id=[cat(1,obj_id{:})];
    lat=[cat(1,lat{:})];
    lon=[cat(1,lon{:})];
    time=[cat(1,time{:})];
    intensity=[cat(1,intensity{:})];
    t1=t+1;
    
    disp('Exporting Tables');
    fn=[Results_path sprintf('Object_Table/Object_Table%d.asc',t+1)];
    fid=fopen(fn,'w');
    fprintf(fid,'%i,%8.4f,%9.4f, %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d, %8.4f\n',[ obj_id'; lat';lon';time(:,1)';time(:,2)';time(:,3)';time(:,4)';time(:,5)';time(:,6)';intensity']);
    fclose(fid);
    
    %% Exporting Tracking Tables
    
    formatSpec1='(%i,ST_GeomFromText(''LINESTRING(';
    formatSpec2='%7.3f %7.3f,';
    formatSpec3='%7.3f %7.3f)'',4326)),\n';
    formatSpec4='%7.3f %7.3f)'',4326))';
    
    hmin=min(Track(:,1));
    hmax=max(Track(:,1));
    h=hmax-hmin+1;
    
    fn=[Results_path sprintf('Track/Track%d.sql',t+1)];
    fid_track=fopen(fn,'w');
    for i=1:h
        idx2=Track(:,1)==(i+hmin-1);
        Tracki=Track(idx2,:);
        
        fprintf(fid_track,formatSpec1,Tracki(1,1));
        [l k]=size(Tracki);
        if i<h
            if l>1
                for j=1:l-1
                    fprintf(fid_track,formatSpec2,Tracki(j,3),Tracki(j,2));
                end
                fprintf(fid_track,formatSpec3,Tracki(j+1,3),Tracki(j+1,2));
            else
                fprintf(fid_track,formatSpec2,Tracki(1,3),Tracki(1,2));
                fprintf(fid_track,formatSpec3,Tracki(1,3),Tracki(1,2));
            end
        else
            if l>1
                for j=1:l-1
                    fprintf(fid_track,formatSpec2,Tracki(j,3),Tracki(j,2));
                end
                fprintf(fid_track,formatSpec4,Tracki(j+1,3),Tracki(j+1,2));
            else
                fprintf(fid_track,formatSpec2,Tracki(1,3),Tracki(1,2));
                fprintf(fid_track,formatSpec4,Tracki(1,3),Tracki(1,2));
            end
            
        end
    end
    fclose(fid_track);
    
    %% Writing tracking for DWR
    disp('Exporting tracking for DWR');
    fn=[Results_path sprintf('Track_DWR/Track_DWR%d.asc',t+1)];
    fid_track_DWR=fopen(fn,'w');
    fprintf(fid_track_DWR,'%i,%8.4f,%9.4f, %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d\n',[Track(:,1)'; Track(:,2)';Track(:,3)';Track(:,4)';Track(:,5)';Track(:,6)';Track(:,7)';Track(:,8)';Track(:,9)']);
    fclose(fid_track_DWR);
    
    %% Calculating Object Boundaries
    disp('Calculating object boundary');
    
    Boundary={};
    
    for i=1:length(vox1)
        
        dat_vox1 = zeros(n_row,(n_col + n_colx));
        
        temp=ceil( mod(vox1{i},n_sizex)/n_row);
        temp(temp==0)=n_col+n_colx;
        vox_x= temp;
        
        vox_y= mod( mod(vox1{i},n_sizex),n_row);
        vox_y(vox_y==0)=n_row;
        
        for j=1:length(vox_x)
            dat_vox1(vox_y(j),vox_x(j))=1;
        end
        
        [B,L] = bwboundaries(dat_vox1,'noholes');
        m0=1;n0=1;
        for m=1:length(B)
            if length(B{m})>n0
                n0=length(B{m});
                m0=m;
            end
        end
        
        K = B{m0};
        n_points = length(B{m0});
        temp = zeros(n_points,3);
        
        temp(:,1)=id2;
        temp(:,2) = 90 - K(:,1)/n_row*180;                      % Lat
        temp(:,3) = -resx + K(:,2)*360/n_col;                      % Lon
        
        if i==1
            Boundary=temp;
        else
            Boundary=[Boundary; temp];
        end
        
        id2=id2+1;
    end
    
    
    %% Exporting object boundary tables
    disp('Exporting object boundary tables');
    
    formatSpec1='(%i,ST_MakePolygon(ST_GeomFromText(''LINESTRING(';
    formatSpec2='%7.3f %7.3f,';
    formatSpec3='%7.3f %7.3f';
    formatSpec4=')'',4326)),''%4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d'',';
    formatSpec5='''%4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d''),\n';
    formatSpec6='''%4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d'')';
    
    hmin=min(Track(:,1));
    hmax=max(Track(:,1));
    h=hmax-hmin+1;
    
    fn=[Results_path sprintf('Boundary/Boundary%d.sql',t+1)];
    fid_boundary=fopen(fn,'w');
    
    for i=1:h
        idx2=find(Track(:,1)==(i+hmin-1));
        
        Tracki=Track(idx2,:);
        
        idx3=find(Boundary(:,1)==(i+hmin-1));
        Boundaryi=Boundary(idx3,:);
        
        fprintf(fid_boundary,formatSpec1,Boundaryi(1,1));
        
        [l k]=size(Boundaryi);
        [g k]=size(Tracki);
        if i<h
            if l>1
                for j=1:l-1
                    fprintf(fid_boundary,formatSpec2,Boundaryi(j,3),Boundaryi(j,2));
                end
                fprintf(fid_boundary,formatSpec3,Boundaryi(l,3),Boundaryi(l,2));
                fprintf(fid_boundary,formatSpec4,Tracki(1,4),Tracki(1,5),Tracki(1,6),Tracki(1,7),Tracki(1,8),Tracki(1,9));
                fprintf(fid_boundary,formatSpec5,Tracki(g,4),Tracki(g,5),Tracki(g,6),Tracki(g,7),Tracki(g,8),Tracki(g,9));
            else
                fprintf(fid_boundary,formatSpec4,Tracki(1,4),Tracki(1,5),Tracki(1,6),Tracki(1,7),Tracki(1,8),Tracki(1,9));
                fprintf(fid_boundary,formatSpec5,Tracki(g,4),Tracki(g,5),Tracki(g,6),Tracki(g,7),Tracki(g,8),Tracki(g,9));
            end
        else
            if l>1
                for j=1:l-1
                    fprintf(fid_boundary,formatSpec2,Boundaryi(j,3),Boundaryi(j,2));
                end
                fprintf(fid_boundary,formatSpec3,Boundaryi(l,3),Boundaryi(l,2));
                fprintf(fid_boundary,formatSpec4,Tracki(1,4),Tracki(1,5),Tracki(1,6),Tracki(1,7),Tracki(1,8),Tracki(1,9));
                fprintf(fid_boundary,formatSpec6,Tracki(g,4),Tracki(g,5),Tracki(g,6),Tracki(g,7),Tracki(g,8),Tracki(g,9));
            else
                fprintf(fid_boundary,formatSpec4,Tracki(1,4),Tracki(1,5),Tracki(1,6),Tracki(1,7),Tracki(1,8),Tracki(1,9));
                fprintf(fid_boundary,formatSpec6,Tracki(g,4),Tracki(g,5),Tracki(g,6),Tracki(g,7),Tracki(g,8),Tracki(g,9));
            end
            
        end
    end
    fclose(fid_boundary);
    
end
else % if ~isempty( files )
        disp( 'All individual timeframe files have been processed.' );

        msg = sprintf( 'Sleeping for %d seconds.', runPeriod );
        disp( msg );

        pause( runPeriod );
end % if ~isempty( files )
end % 'for vs = 1:length(varnames)' : iteration over files with single variable
end % 'while true' loop

toc
