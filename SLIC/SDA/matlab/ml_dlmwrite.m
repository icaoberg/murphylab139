function ml_dlmwrite(filename, m, dlm, r, c)
% ML_DLMWRITE Write ASCII delimited file.
%
% ML_DLMWRITE( FILENAME, M, DLM, R, C)
%    
%   This is like matlab's standard dlmwrite, except
%   ELEMENTS WHOSE VALUE IS 0 WILL NOT BE OMITTED, and
%   up to 7 significant digits are used (instead of 5).
%   The function returns 1 if the file was successfully written, or
%   0 if the file cannot be written for any reason (e.g. not enough
%   permissions, or disk full or whatever).
%
%   See also DLMWRITE, DLMREAD, CSVREAD, WK1READ, WK1WRITE.

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

%   Brian M. Bourgault 10/22/93
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2006/06/27 13:33:47 $
%
    
%
% test for proper filename
%
if ~isstr(filename),
    error('FILENAME must be a string.');
end;

if nargin < 2, error('Requires at least 2 input arguments.'); end

NEWLINE = sprintf('\n');

% delimiter defaults to Comma for CSV
if nargin < 3, dlm = ','; end
dlm = sprintf(dlm); % Handles special characters.


% open the file
if strncmp(computer,'MAC',3)
  fid = fopen(filename ,'wt');
else
  fid = fopen(filename ,'wb');
end

if fid == (-1), error(['Could not open file ' filename]); end

% check for row,col offsets
if nargin < 4, r = 0; end
if nargin < 5, c = 0; end

% dimensions size of matrix
[br,bc] = size(m);

% start with offsetting row of matrix
for i = 1:r
    for j = 1:bc+c-1
        fwrite(fid, dlm, 'uchar');    
    end
    fwrite(fid, NEWLINE, 'char');
end

% start dumping the array, for now number format float
for i = 1:br

    % start with offsetting col of matrix
    for j = 1:c
        fwrite(fid, dlm, 'uchar');    
    end

    for j = 1:bc
        %if(m(i,j) ~= 0)
            str = num2str(m(i,j),'%0.7g');
            fwrite(fid, str, 'uchar');    
        %end
        if(j < bc)
            fwrite(fid, dlm, 'uchar');    
        end
    end
    fwrite(fid, NEWLINE, 'char'); % this may \r\n for DOS 
end

% close files
fclose(fid);
