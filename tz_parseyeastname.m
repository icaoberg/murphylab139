function [chr,fileprefix,channel] = tz_parseyeastname(yeastname)
%TZ_PARSEYEASTNAME parse the combined yeast name.
%   CHR = TZ_PARSEYEASTNAME(YEASTNAME) parse the string returned from
%   TZ_YEATIMG2FEATNAME. It returns the chormosome name.
%   
%   [CHR,FILEPREFIX,CHANNEL] = TZ_PARSEYEASTNAME(...) also returns file 
%   prefix and channel of the [yeast image]. 
%   
%   See also

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

%   16-Feb-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

dotPosition = find(yeastname=='.');
if ~isempty(dotPosition)
    nameWithNoExtension = yeastname(1:dotPosition(1)-1);
else
    nameWithNoExtension = yeastname;
end
separator = '_';
tokens = tz_strtok(nameWithNoExtension,separator);

chr = tokens{1};
fileprefix = tz_cell2str(tokens(2:end-1),separator);
channel = tokens{end};