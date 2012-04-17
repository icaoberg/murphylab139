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

addpath ('/home/ejr/ml/SImEC/matlab');
addpath ('/home/ejr/ml/SImEC/matlab/mex');

im_dir1 = '/imaging/ejr/images/giant/prot';
im_dir2 = '/imaging/ejr/images/gpp130/prot';
crop_dir1 = '/imaging/ejr/images/giant/crop';
crop_dir2 = '/imaging/ejr/images/gpp130/crop';
dna_dir1 = '/imaging/ejr/images/giant/dna';
dna_dir2 = '/imaging/ejr/images/gpp130/dna';
featsetname = 'SLF7dna';
conf = 0.05;
saver = 1;
im_name1 = 'giantin';
im_name2 = 'gpp130';

fprintf('\nSImEC Demo.\n');
fprintf('   A demonstration of SImEC comparing giantin and \n');
fprintf('   gpp130 with their corresponding DNA and crop \n'); 
fprintf('   images using nodna_nohar feature set.\n');


[res, F, crit, both, features1, features2, slf_names, names, start_now] = ...
    ml_simec(im_dir1, im_dir2, crop_dir1, crop_dir2, dna_dir1, dna_dir2, ...
	     featsetname, conf, saver, im_name1, im_name2)
