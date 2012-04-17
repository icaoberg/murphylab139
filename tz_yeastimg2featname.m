function featname = tz_yeastimg2featname(imgpath)
%TZ_YEASTIMG2FEATNAME Image file path to feature file name
%   FEATNAME = TZ_YEASTIMG2FEATNAME(IMGPATH) converts the path of a [yeast
%   image] into a unique file name with the information of chromosome and
%   image file name.
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

%   12-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

imagePathTokens = tz_strtok(imgpath,filesep);
chromosomeName = imagePathTokens{end-1};
filename = strtok(imagePathTokens{end},'.'); %without extension

featname = [chromosomeName '_' filename];
