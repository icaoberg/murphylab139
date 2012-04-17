% This script is creating tables and figures in the ISMB paper

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

% make data/procfiles.mat to store yeast images information
make_procfile

% creates the field level confusion matrix for Table 2 (data/Table2.html)
createfieldconfusionmat

% generate figures in Figure 3 (data/Fig3.tif), Figure 4 (data/Fig4.tif), and Figure 5 (data/Fig5.tif)
fig345

% generate the confusion matrix for cell level classification (data/Table3.html)
% draw precision-recall plot (data/Fig2.tif)
% create Supplementary material for list of 501 proteins whose label from visual assignment differs from that by automated classification
supplement1

% generate supplementary material for a list proteins in the ambiguous category (data/supplement2.txt)
% and a list proteins in the punctate_composite category (data/supplement3.txt)
supplement23

% summarize the classification result
ClassifySummary
