function edgeimg = tz_imstdedge(img)
%TZ_IMCORREDGE
%   EDGEIMG = TZ_IMCORREDGE(IMG) returns the edge image 
%   based on identifying the location which has high
%   variance of pixel intensity
%
%   See also sc_cutcells

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

if nargin < 1
    error('Exactly 1 argument is required')
end
img = double(img);
maskimg = fspecial('gaussian',5,2);
% meandicimg=filter2(maskimg>0,img,'same')./sum(maskimg(:)>0);
meandicimg=imfilter(img,maskimg,'symmetric','same');
shiftdicimg = img-meandicimg;
vardicimg=double(shiftdicimg).*double(shiftdicimg);
vardicimg=imfilter(vardicimg,maskimg,'symmetric','same');
edgeimg=sqrt(vardicimg)./meandicimg;

edgeimg = uint8(round(mat2gray(edgeimg)*255));
thresh = ml_rcthreshold(edgeimg);
edgeimg(edgeimg<thresh)=0;
