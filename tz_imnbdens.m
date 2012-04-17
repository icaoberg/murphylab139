function fs = tz_imnbdens(img,wndlength)
%TZ_IMNBDENS Total intensities between every pixel pair in an image.
%   FS = TZ_IMNBDENS(IMG,WNDLENGTH) returns a 4-D matrix with dimension row X col X W X W. 
%   to access the value at (50,20), use squeeze(edgepot(50,20,:,:))
%   edge potential is calculated as sum of intensities along the line between
%   the center nodes and the neighbor node
%   FS(i,j,u,v) is the total intensity between pixel IMG(i,j) and 
%   IMG(i+u-wndlength,j+v-wndlength). The neighbor window for each pixel
%   is (2*wndlength+1)x(2*wndlength+1)
%   edgepot is a 4-D matrix with dimension row X col X W X W
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

%   25-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

img = double(img);

wndSize = [2*wndlength+1,2*wndlength+1];
[filters,offsets]=tz_makelinefilters(wndlength);

fs = zeros(size(img,1),size(img,2),wndSize(1),wndSize(2));

for k=1:length(filters)
    fs(:,:,offsets(k,1)+wndlength+1,offsets(k,2)+wndlength+1) = ...
        imfilter(img,filters{k},'same');
end