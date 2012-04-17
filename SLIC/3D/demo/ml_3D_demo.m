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

addpath('/home/ejr/ml/3D/matlab/mex');
addpath('/home/ejr/ml/3D/matlab');
imName = '/home/ejr/ml/3D/Tub/Tub20010319ss09s08q1i050.tif';
featName = [imName '.feat.mat'];
imgs = ml_tiffread(imName, [2 3]);
[masks] = ml_3dsegment_cells(imgs{2}, imgs{1});
compressed_masks{1} = ml_rle(masks{1});
save([imName '.crop.mat'], 'compressed_masks');
ml_3dfeatset(imName, featName, 'Prot', [], [])
