function regions = tz_maskimgs(img,masks)
%TZ_MASKIMGS Extract masked regions from an image
%   REGIONS = TZ_MASKIMGS(IMG,MASKS) returns a cell array of images that
%   are regions extracted from IMG specified by MASKS. If MASKS is also a 
%   cell array, then each element of MASKS should be a matrix that 
%   is as large as IMG. A sparse matrix is also allowable. Any pixels no
%   greater than 0 will be taken as background. If MASKS is a matrix, the
%   regions will be defined as a group of pixels that have the same gray
%   level, which must be a postitive integer. The integers must be
%   continuous and start from 1.
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

%   10-Nov-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

if iscell(masks) %cell array mask format
    nregion = length(masks);
else %matrix mask format
    if any(size(masks)~=size(img))
        error('Size unmatched.');
    end 
    nregion = double(max(masks(:)));
end

regions = {};

for k=1:nregion
    img2 = img;
    if iscell(masks)
        mask = masks{k};
    else
        mask = (masks==k);
    end
    img2(mask<=0) = 0;
    regions{end+1} = tz_imcropbg2(img2,mask);
end
