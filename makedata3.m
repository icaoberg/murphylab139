% This script is calculating field level features

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

fprintf(1,'Calculating field level features...\n\n');

addpath(genpath([pwd filesep 'SLIC']));
method.procfilesmatname = [pwd filesep 'data' filesep 'procfiles.mat'];
method.celllevel = 0;
method.har_intbins = 128;   
method.downsamplerate = 1;
sc_yeastfeatcalc(method);

fprintf(1,'Field level feature calculation is done!\n\n');
