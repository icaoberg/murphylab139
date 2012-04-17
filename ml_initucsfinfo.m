%This script loads UCSF labeling information into yeastdata

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

disp('Loading UCSF labels from allOrfData.txt...');
ml_makeucsfyeaststruct

allFields = fieldnames(yeaststruct);

yeastdata.ucsf.locationNames = allFields(10:end);

disp('Creating label matrix ...');
yeastdata.ucsf.labelMatrix = [];
for i=1:length(yeastdata.ucsf.locationNames)
    curLocLabels = getfield(yeaststruct,yeastdata.ucsf.locationNames{i});
    yeastdata.ucsf.labelMatrix(:,i) = (cell2mat(curLocLabels)=='T')';
end

idx = ml_strmatch(yeastdata.geneNames,yeaststruct.yORF);

yeastdata.ucsf.labelMatrix = yeastdata.ucsf.labelMatrix(idx,:);
% ucsfLabelNums = sum(ucsfLabelMatrix,2);
yeastdata.ucsf.locNums = sum(yeastdata.ucsf.labelMatrix,2);

[ucsfUniqueMixLabelMatrix,tmplabels,yeastdata.ucsf.classLabels] = ...
    unique(yeastdata.ucsf.labelMatrix,'rows');
yeastdata.ucsf.classLocNums = sum(ucsfUniqueMixLabelMatrix,2);

ucsfMcfImageIndices = ...
    ml_combfeats2mcf([1:length(procfiles{1})]',yeastdata.ucsf.classLabels);

disp('Making classes ...');
for i=1:size(ucsfUniqueMixLabelMatrix,1)
    labelCode = ucsfUniqueMixLabelMatrix(i,:);
    if all(labelCode==0) %For unvisualized classes
        yeastdata.ucsf.strClassNames{i} = 'unvisualized';
    else
        yeastdata.ucsf.strClassNames{i} = ...
            tz_cell2str(yeastdata.ucsf.locationNames( ...
            find(labelCode>0)),',');
    end
end

disp('All information has been stored in the structure yeastdata');
