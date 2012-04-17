function img2 = tz_normreg(img)
%TZ_NORMREG Convert a [region image] to a [normalized region image].
%   IMG2 = TZ_NORMREG(IMG) returns a [normalized region image] which
%   has the same regions as those in the [region image] IMG.
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

%   23-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

bwimg = img>0;
newLabeledImage = bwlabel(bwimg);

newLabeledImage = immultiply(newLabeledImage, ...
    zeros( size(newLabeledImage) )+double( max( img(:) ) ) );

img2 = imadd(double(img),newLabeledImage);

eval(['img2 = ' class(img) '(img2);']);

img2(img2>0) = tz_squeezenum(img2(img2>0));
