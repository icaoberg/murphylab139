function status = ml_createfile_excl( filefullpath) 

% STATUS = ML_CREATEFILE_EXCL( FILEFULLPATH)
%
% Create the file FILEFULLPATH only if it does not already exist.
% Return value is 0 if file was successfully created, or 1 if it
% already exists or 2 if it could not be created for some other
% reason (such as non-existent directory or not enough permission).

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

% Meel Velliste 5/20/02
% Modified 3/03 EJSR
%     binpath
% Modified 7/05 T. Zhao
%     bincmd

bincmd=which('ml_createfile_excl');
slashpos=find(bincmd=='/');
bincmd=[bincmd(1:slashpos(end-1)) 'bin' bincmd(slashpos(end):end-2)];

[status,output] = unix([bincmd ' ' filefullpath]);

if status~=0
    warning('Failed to create the file');
end