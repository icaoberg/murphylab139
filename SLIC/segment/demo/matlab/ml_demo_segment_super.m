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


dir1 = '/home/ejr/atto/02-12-12-ChoCells-16Bit-Confocal/';
dir2 = '/home/ejr/atto/02-12-19_ChoCells/confocal/';
dnaIm = 'E1 confocal Hoechst 20x.tif';
mitIm = 'E1 confocal MitoTracker Red 20x.tif';
alxIm = 'E1 confocal AlexaPhalloidin 20x.tif';

addpath(genpath('/home/ejr/ml/input/'));

addpath('/home/ejr/ml/segment/matlab/');
addpath('/home/ejr/ml/segment/matlab/mex/');

dna = ml_readimage([dir1 'multi/' dnaIm]);
im = ml_readimage([dir1 'multi/' alxIm]);

[lines1, masks] = ml_segment_cells_voronoi_super(im, dna);


figure, imshow(lines1, []);

