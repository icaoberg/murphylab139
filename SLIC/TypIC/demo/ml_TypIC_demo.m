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

 
addpath(genpath('/home/ejr/ml/input'))
addpath(genpath('/home/ejr/ml/featcalc'))
prot = '/imaging/mb_hela_images/giant/prot';
crop = '/imaging/mb_hela_images/giant/crop';
disp('TypIC demo extracts features and ranks images');
[s, i, d] = ml_TypIC2(-1, prot, crop, 'nodna_nohar', 'none', 0, 0, ...
		      'false', prot, crop)
