function yeastdata = ml_uniyeastinfo(yeastdata)
%ML_UNIYEASTINFO Remove duplicated ORF from yeast information. 
%   YEASTDATA = ML_UNIYEASTINFO(YEASTDATA) returns a structure that contain
%   unique ORFs in the structure YEASTDATA.
%   
%   See also

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

%   14-Feb-2007 Initial write T. Zhao
%   Copyright (c) 2007 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 1
    error('Exactly 1 argument is required');
end

[names,indices] = unique(yeastdata.geneNames);

yeastdata.geneNames = yeastdata.geneNames(indices);
yeastdata.imageIndices = 1:length(indices);
for i=1:length(yeastdata.procfiles)
    yeastdata.procfiles{i} = yeastdata.procfiles{i}(indices);
end

if isfield(yeastdata,'ucsf')
    yeastdata.ucsf.labelMatrix = yeastdata.ucsf.labelMatrix(indices,:);
    yeastdata.ucsf.locNums = yeastdata.ucsf.locNums(indices);
    yeastdata.ucsf.classLabels = yeastdata.ucsf.classLabels(indices);
end

if isfield(yeastdata,'cygd')
     yeastdata.cygd.locNums = yeastdata.cygd.locNums(indices);
     yeastdata.cygd.classLabels = yeastdata.cygd.classLabels(indices);
end

