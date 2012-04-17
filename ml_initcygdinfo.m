%This script is to initialize location information of the [yeast image]s.
%The information is from CYGD.
%Notice: here gene name means an ORF name
%Created variables:
%   cygdinfo (locs previously) - a structure of the list of gene names and 
%       their corresponding location labels:
%           geneNames - 
%       
%   yeastdata - This is a structure containing information about yeast data.
%       It is created in tz_run_inityeastinfo. But in this scripct a field
%       called cygd will be added. yeastdata.cygd is also a structure and
%       has the following fields:
%           'locLabels' - a cell array of cell array. Each element is the 
%               assigned [location label]s of the corresponding gene.
%           'locNums' - the number of possible locations for each gene.
%           'classLabels' - a [label vector] based on locLabels. Each  
%               element is the class label of the corresponding [yeast image].             
%           'mainLabels' - a cell array of cell array. Each element is the 
%               assigned main [location label]s of the corresponding gene.
%           'strLocLabels' - a [string array]. Each element is a string of
%              the combination of [location label]s.
%           'allLocLabels' - all unique combination of [location labels]s.
%           'classNames' - a cell array of cell array. Each element is the 
%               cell array of [location name]s of the corresponding class.
%           'strClassNames' - a [string array]. Each element is a string of
%               the combination of [location name]s of the corresponding
%               class.
%           'classLocNums' - number of main locations of each class.
%           'classGeneNums' - Number of genes in each class.
%   matchedIndices - the indices of mathching cygdinfo.geneNames and 
%       geneNames. The ith element of geneNames is the same as       
%       cygdinfo.geneNames{matchedIndices(i)} except that matchedIndices(i)
%       is 0, which means that there is no such match.
%   invalidIndices - indices of unmatched gene names

tz_run_initcygdinfo_status = 0;

%load location labels
if ~exist('cygdinfo','var')
    cygdinfo = tz_loadlocdata('data/subcellulardata.txt');
end


reply = 'n';
if ~exist('isinteractive','var')
    isinteractive = 0;
end
if isinteractive==1
    reply = input( ['Do you want to exclude images'  ...
        'without significant signal? y/n [n]: '],'s');
end


%Match ORF names with location labels
matchedIndices = strmatch_multi(yeastdata.geneNames,cygdinfo.geneNames);

%Find ORFs that do not have assigned location labels
invalidIndices = find(matchedIndices==0);

yeastdata.cygd = struct('islocunique',[]);

%Find proteins that have only one location label
for k=1:length(matchedIndices)
    if matchedIndices(k)>0
        [islocunique(k),yeastdata.cygd.locLabels{k}, ...
            yeastdata.cygd.mainLabels{k}] = ...
            tz_islocunique(cygdinfo.locationLabels{matchedIndices(k)});       
        yeastdata.cygd.locNums(k) = length(yeastdata.cygd.locLabels{k});
    else
        islocunique(k) = 0;
        yeastdata.cygd.locLabels{k} = {'NaN'}; %Not Available
        yeastdata.cygd.mainLabels{k} = {'NaN'};
        yeastdata.cygd.locNums(k) = 0;
    end
    yeastdata.cygd.strLocLabels{k} = ...
        tz_cell2str(yeastdata.cygd.locLabels{k});
end

[yeastdata.cygd.allLocLabels,tmplabels,yeastdata.cygd.classLabels] = ...
    unique(yeastdata.cygd.strLocLabels);
yeastdata.cygd.classLabels = yeastdata.cygd.classLabels';

%Load location schemes
cygdinfo.locSchemes = tz_loadlocschm('data/schemarr.txt');

sortedLabels = ml_combfeats2mcf(yeastdata.cygd.classLabels, ...
    yeastdata.cygd.classLabels);

mcfImageIndices = ml_combfeats2mcf([1:length(procfiles{1})]', ...
    yeastdata.cygd.classLabels);

for i=1:length(yeastdata.cygd.allLocLabels)
    yeastdata.cygd.classNames{i} = ...
        tz_getlocnames(cygdinfo.locSchemes, ...
        tz_strtok(yeastdata.cygd.allLocLabels{i},','));
    yeastdata.cygd.strClassNames{i}  = ...
        tz_cell2str(yeastdata.cygd.classNames{i});
    if strcmp(yeastdata.cygd.strClassNames{i},'unknown')
        yeastdata.cygd.classLocNums(i) = 0;
    else
        yeastdata.cygd.classLocNums(i) = ...
            length(yeastdata.cygd.classNames{i});
    end
    yeastdata.cygd.classGeneNums(i) = size(sortedLabels{i},1);
end

tz_run_initcygdinfo_status = 1;
