function write_bn(fn,dat,sz)
% savebfn -- saves a binary file fn of dim
%       ([nrows ncols]) of sz 'float','int', etc.
%
%       by DKB on 2/5/96
%

if(nargin < 3)
        sz = 'float32';
end

f = fopen(fn,'w','l');
fwrite(f,dat',sz);
fclose(f);
