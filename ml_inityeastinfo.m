%This script initializes location information of the yeast images.
%Notice: here gene name is also an ORF name
%Created a structure variable yeastdata with the following fields:
%    procfiles - the list of image file names, including 3 channels.
%            procfiles{1} is a [string array] of all DAPI image paths
%            procfiles{2} is a [string array] of all DIC image paths
%            procfiles{3} is a [string array] of all GFP image paths
%    geneNames - the list of gene(ORF) names of images.
%    imageIndices - the indices of images.

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

disp('Load yeast image file paths ...');
load(['data' filesep 'procfiles.mat']);

if exist('yeastdata','var')
    rmfield(yeastdata,'geneNames');
end

%Get ORF names from the files
for k=1:length(procfiles{1})    
    yeastdata.geneNames{k} = tz_getgenename(procfiles{1}{k});
end

yeastdata.imageIndices = 1:length(procfiles{1});
yeastdata.procfiles = procfiles;

disp('Creating yeastdata.mat...');



