function locs = tz_loadlocdata(filepath)
%TZ_LOADLOCDATA Load location data from a file.
%   LOCS = TZ_LOADLOCDATA(FILEPATH) returns a structure with two fields:
%       geneName: a [string array] of ORF names
%       locationLabels: a [string array] of location labels
%   To get location label for a given ORF name GENE, first search the name 
%   in LOCS.geneName and get the index I, i.e. LOCS.genename(I) is the same
%   as GENE, then LOCS.locationLabels(I) are the labels.
%   
%   Notice:
%       One orf could have multiple labels.
%   
%   See also

%   04-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

str = textread(filepath,'%s','delimiter','\n');
inc = 1;
for k=1:length(str)
    linepos = find(str{k}=='|');
    if length(linepos)>0
        tokens = tz_strtok(str{k},'|');
        if length(tokens)>1
            if str{k}(linepos(1)+1)~='|'
                geneNames{inc} = upper(tokens{1});
                locationLabels{inc} = tokens{2};
                inc = inc+1;
            end
        end
    end
end

[locs.geneNames,indices] = unique(geneNames);
[indices rankIndices] = sort(indices);
locs.geneNames = upper(locs.geneNames(rankIndices));
indices = [0 indices];

for k=1:length(indices)-1
    locs.locationLabels{k} = locationLabels(indices(k)+1:indices(k+1));
end