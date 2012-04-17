function RGB = mv_superimpose( img, mask, my_colormap)

% MV_SUPERIMPOSE( IMG, MASK, COLORMAP)
% SUPERIMPOSE MASK on IMG and display it

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

if( nargin < 3)
  my_colormap = gray;
end
  
% Turn img into and RGB image with colormap hot
Max = double(max(img(:)));
img = double(img);
img = round(63*img/Max + 1);
colormap_r = my_colormap(:,1);
colormap_g = my_colormap(:,2);
colormap_b = my_colormap(:,3);
R = colormap_r(img);
G = colormap_g(img);
B = colormap_b(img);

% superimpose the mask image in a distinct color
mask_rgn = find(mask);
R(mask_rgn) = 0;
G(mask_rgn) = 0.5;
B(mask_rgn) = 1;
RGB = cat(3,R,G,B);

% Display the rsultant overlaid image
imagesc(RGB);


function f2

Max = max(img(:));
R = img;
G = img;
B = img;
mask_rgn = find(mask);
R(mask_rgn) = 3*double(Max);
G(mask_rgn) = 0;
B(mask_rgn) = 0;
RGB = cat(3,R,G,B);
imagesc(RGB);

