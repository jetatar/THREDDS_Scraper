function voxx = vox_fn( data, conn, var )
% This code is modifided from bwconncomp.m function in Matlab
% by Phu Nguyen 11/22/2011
% dat: data in 3D matrix format
% conn: connectivity defined to segment objects
% conn= 6            three-dimensional six-connected neighborhood
%      18            three-dimensional 18-connected neighborhood
%      26            three-dimensional 26-connected neighborhood

    fL = isnan( var('valLow') );
    fH = isnan( var('valHi') );

    if not(fL) & not(fH)    % if both limits are set...

        data = data >= var( 'valLow' );
        data = data <= var( 'valHi' );

    elseif not(fL) & fH     % if only lower limit set

        data = data >= var( 'valLow' );

    elseif fL & not(fH)     % if only higher limit set

        data = data <= var( 'valHi' );

    end

    temp = bwconncomp(data,conn);
    voxx = temp.PixelIdxList;
