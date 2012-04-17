% This script generates several intermeidate files about the yeast image information and UCSF labeling
% All these information will be stored into yeastdata.mat that will be used later
% yeastdata.mat has two variables, yeastdata and yeastdataunique, 
% which describe all the information about yeast data with/with unique labels

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

addpath(genpath([pwd filesep 'SLIC']));

% make data/procfiles.mat to store yeast images information
make_procfile

%Initialize location information of the yeast images and save them in a structure variable yeastdata
ml_inityeastinfo

%Load UCSF labeling information into yeastdata
ml_initucsfinfo

%Load CYGD labeling information into yeastdata
ml_initcygdinfo

%Create a structure variable yeastdataunigene for unique gene names:
yeastdataunigene = ml_uniyeastinfo(yeastdata);

%Store the results into yeastdata.mat
save(['data' filesep 'yeastdata.mat'],'yeastdata','yeastdataunigene');

disp(['Initialization done; all yeast information is stored in ' 'data' filesep 'yeastdata.mat']);
