function [filters,offsets] = tz_makelinefilters(wndsize)
%TZ_MAKELINEILTERS Make line filters
%   FILTERS = TZ_MAKELINEFILTERS(WNDSIZE) returns a set of filters which
%   reprepent a set of lines in an image. WNDSIZE determines the size of
%   largest filter, which is (WNDSIZE*2+1)x(WNDSIZE*2+1)
%  
%   [FILTERS,OFFSETS] = TZ_MAKELINEFILTERS(...) also retuns end point of
%   linnes corresponding to filters.
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

%   18-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

wndWidth = wndsize*2+1;
wndHeight = wndWidth;
wndCenter = [wndsize,wndsize];
m=1;
offsets = [];
for j=0:wndWidth-1
    for k=0:wndHeight-1
        currentPoint = [j,k];
        if any(currentPoint~=wndCenter)
            offset = currentPoint-wndCenter;
            offsets =[offsets;offset];
            filtersize = abs(offset)*2+1;
            filter = zeros(filtersize);
            newCenter = abs(offset)+1;
            cornerPoint = newCenter+offset;
            filters{m} = ...
                tz_setimglnpixel(filter,newCenter,cornerPoint);
            m=m+1;
        end
    end
end
