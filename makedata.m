% This script generates all the results

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

% generates several intermeidate files about the yeast image information and UCSF labeling
makedata0

% segment all the yeast images into cell regions
makedata1

% calculate cell level features
makedata2

% calculate field level features
makedata3

% prepare feature matrix for cell level classification 
makedata4

% run cell level classification on 20 unique location classes with purality voting scheme
makedata5

% prepare feature matrix for field level classification 
makedata6

% run field level classification on 4 and 21 unique location classes with 6 fold cross-validation
makedata7

% generate figures and tables
makedata8
