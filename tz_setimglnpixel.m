function img2=tz_setimglnpixel(img,s,t,val)
%TZ_SETIMGLNPIXEL Draw a line in an image.
%   IMG2 = TZ_SETIMGLNPIXEL(IMG,S,T) set all pixels on a line to 1 in the
%   image IMG. The line is from S to T. See ML_GETLINEPTS for S and T.
%   
%   IMG2 = TZ_SETIMGLNPIXEL(IMG,S,T,VAL) set the value VAL instead of 1.
%
%   See also TZ_SETIMGLNPIXEL2 ML_GETLINEPTS

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

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin<4
    val=1;
end

pts=ml_getlinepts(s,t);
img2=ml_setimgptspixel(img,[pts,zeros(size(pts,1),1)+val]);