function image = ml_readimage( filename,tmpdir)

% IMAGE = ML_READIMAGE( FILENAME,TMPDIR)
%
% Reads a 2D image in any format, including TCL files.
% If FILENAME is an emtpy string, then an empty matrix
% is returned.
%

% Copyright (C) 2006  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

% ml_readimage written by Meel Velliste

if(nargin<2)
    tmpdir = tempdir;
end

if( tmpdir(end) ~= filesep) 
    tmpdir(end+1) = filesep; 
end;

if( isempty( filename))
  image = [];
  return;
end

L = length( filename);
extension = filename(L-2:L);


switch( extension)
 case '.gz'
    [f,tempfilename] = ml_fopentemp(tmpdir);
    tempfullpath = [tmpdir tempfilename];
    fclose(f);
    unix(['gunzip -c ' filename ' > ' tempfullpath]);
    image = imread( tempfullpath);
    unix(['rm ' tempfullpath]);
 case 'bz2'
    [f,tempfilename] = ml_fopentemp(tmpdir);
    tempfullpath = [tmpdir tempfilename];
    fclose(f);
    unix(['bunzip2 -c ' filename ' > ' tempfullpath]);
    image = imread( tempfullpath);
    unix(['rm ' tempfullpath]);
 case 'dat'
    image = ml_tclread( filename);
 otherwise
    image = imread( filename);
end

