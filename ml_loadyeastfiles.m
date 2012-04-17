%This script loads all the paths of image files and save them into a cell array stored in procfiles.mat

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

%Get directories of chromosomes
cloneDirs = tz_ls(yeastRoot,'dir');

procfiles = {{},{},{}};

for i=1:length(cloneDirs)
    %DAPI(DNA) channel
    procfiles{1} =  ml_cellcat({procfiles{1}' ...
        tz_ls([yeastRoot filesep cloneDirs{i} ...
        filesep '*_DAPI.png'],'all -r')'})';

    %DIC(cell) channel
    procfiles{2} = ml_cellcat({procfiles{2}' ...
        tz_ls([yeastRoot filesep cloneDirs{i} ...
        filesep '*_DIC.png'],'all -r')'});

    %GFP(protein) channel
    procfiles{3} = ml_cellcat({procfiles{3}' ...
        tz_ls([yeastRoot filesep cloneDirs{i} ...
        filesep '*_GFP.png'],'all -r')'});
end
tz_save(['.' filesep 'data' filesep 'procfiles.mat'],{procfiles},{'procfiles'},'ml_loadyeastfiles', ...
        'yeast image file paths');
