function [links, outarr] = ml_readconsensus(inputfile, feat, method, weight)
% FUNCTION [LINKS, OUTARR] = ML_READCNSENSUS(INPUTFILE, FEAT, METHOD, WEIGHT)
% Read in the majority consensus tree output from the Java program and 
% calculate the distance for clones.  The output wuld be used by 
% ml_drawconsensus function.
% inputfile: the output of the Majority consensus Java program
% feat: the features used in creating the set of dendrograms
% method: 0 for 'euclidean' and 1 for 'mahalanobis'
% weight: an optional array for weights, only apply to Euclidean distance
% Xiang Chen Jan 4, 2005

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


     %Read in the text from input file
     f = fopen(inputfile);
     alltext = fread(f);
     alltext = char(alltext');

     %Read in taxa labels
     taxpos = findstr(alltext, 'taxlabels');
     endpos = findstr(alltext, 'end;');
     labeltext = alltext(taxpos + 10 : endpos(1) - 3);
     spacepos = [1 findstr(labeltext, ' ') length(labeltext)];
     spacepos(end) = spacepos(end) + 1;
     for m = 1 : length(spacepos) - 1
         text = labeltext(spacepos(m) : spacepos(m + 1) - 1);
         text(findstr(text, ' ')) = [];
         text(findstr(text, char(13))) = [];
         text(findstr(text, char(9))) = [];
         labels{m} = text;
     end

     %Read in tree structures
     treetext = alltext(findstr(alltext, '[&R]') + 6 : endpos(2) - 3);
     returnpos = findstr(treetext, char(13));
     treetext(returnpos) = [];
     tabpos = findstr(treetext, char(9));
     treetext(tabpos) = [];

     ll = length(labels);
     for m = 1 : ll
         textpos = findstr(treetext, [labels{m}, ',']);
         if (isempty(textpos))
             textpos =  findstr(treetext, [labels{m} ')']);
         end
    	 treetext(textpos:textpos + length(num2str(m)) - 1) = num2str(m);
	     treetext(textpos + length(num2str(m)) : textpos + length(labels{m}) - 1) = [];
     end
     treetext(findstr(treetext, '>')) = [];
     
  
     allfeat = [];
     for m = 1 : ll
         arr{m} = m;
         allfeat = [allfeat
             feat{m}];
         meanfeat(m,:) = mean(feat{m}, 1);
     end
     meanf = meanfeat;
     zs = zscore(allfeat);
     
     mf = mean(allfeat, 1);
     sf = std(allfeat);
     sf(find(sf==0)) = 1;
     [text, id, outarr] = xc_consensus2array(arr, treetext);
     for m = (ll + 1) : (length(outarr))
         meanfeat(m, :) = mean(meanfeat(outarr{m}, :), 1);
         arr{m} = arr{outarr{m}(1)};
         for n = 2 : length(outarr{m})
             arr{m} = [arr{m} arr{outarr{m}(n)}];
         end
         meanf(m, :) = mean(meanfeat(arr{m}, :), 1);
     end
     meanfeat = (meanfeat - repmat(mf, [size(meanfeat, 1) 1])) ./ repmat(sf, [size(meanfeat, 1) 1]);
     meanf = (meanf - repmat(mf, [size(meanfeat, 1) 1])) ./ repmat(sf, [size(meanfeat, 1) 1]);
     if (exist('method', 'var') & (~isempty(method)) & method)
         if (method == 1)
             %dists = squareform(pdist(meanfeat, 'mahal')); 
             dists = squareform(ml_pdist(meanfeat, 'mah', (cov(zs))));
             %dist2 = squareform(pdist(meanf, 'mahal')); 
             dist2 = squareform(ml_pdist(meanf, 'mah', (cov(zs))));
         elseif (method == 2)
             dists = squareform(ml_pdist(meanfeat, 'city'));
             dist2 = squareform(ml_pdist(meanf, 'city'));
         elseif (method == 3)
             dists = squareform(ml_pdist(meanfeat, 'ca'));
             dist2 = squareform(ml_pdist(meanf, 'ca'));
         end
     else
         if (exist('weight', 'var'))
             meanfeat = meanfeat .* repmat(weight, [size(meanfeat, 1) 1]);
             meanf = meanf .* repmat(weight, [size(meanf, 1) 1]);
         end
         dists = squareform(pdist(meanfeat));
         dist2 = squareform(pdist(meanf));
     end
     
     for m = (ll + 1) : (length(outarr))
         d = dist2(outarr{m}, outarr{m});
         as = sum(d(:)) / (length(outarr{m}) * (length(outarr{m}) - 1)) / 2;
         s = [];
         smax = [];
         if (length(outarr{m}) <= 2)
             s(1:length(outarr{m})) = as; % / length(outarr{m});;
             smax = s;
             for n = 1 : length(outarr{m})
                 if (outarr{m}(n) > ll)
                     s(n) = s(n) - links(outarr{m}(n) - ll).Average;
                     if (s(n) < 0)
                         s(n) = 0.01;
                     end
                     %if ((s(n) + links(outarr{m}(n) -ll).Max) > smax(n))
                     smax(n) = (s(n) + links(outarr{m}(n) -ll).Max);
                     %end
                 end
             end
         else
             for n = 1 : length(outarr{m})
                 s(n) = (sum(d(n, :)) - sum(d(:)) / (2 * length(outarr{m}) - 2)) / (length(outarr{m}) - 2);
                 if (outarr{m}(n) > ll)
                     s(n) = s(n) - links(outarr{m}(n) - ll).Average;
                     if (s(n) < 0)
                         s(n) = 0.01;
                     end
                     smax(n) = (s(n) + links(outarr{m}(n) -ll).Max);
                 else
                     smax(n) = s(n);
                 end

             end
         end
         d1 = dists(arr{m}, arr{m});
         s1 = sum(d1(:)) / (length(arr{m}) * (length(arr{m}) - 1));
         d2 = dist2(outarr{m}, outarr{m});
         s2 = sum(d2(:)) / (length(outarr{m}) * (length(outarr{m}) - 1));
         links(m - ll) = struct('Clusters', outarr{m}, 'Distance', s, 'Distance_on_clone', s1, 'Weighted_distance', s2, 'Average', as, 'Max', max(smax));
     end
     
function tree = xc_link2tree(links, names)
% FUNCTION TREE = XC_LINK2TREE(LINKS, NAMES)
% Convert a linkage matrix to a NEWICK format tree
% links: the linkage matrix from linkage function for dendrogram or from
%        xc_neighborjoining for additive tree
% names: a cell array containing the taxa name for each leave
% tree: NEWICK formatted output tree

taxadimension = size(links, 1) + 1;

if (~exist('names', 'var'))
     for m = 1 : taxadimension
         names{m} = num2str(m);
     end
end

if (length(names) < taxadimension)
     for m = length(names) + 1 : taxadimension
         names{m} = num2str(m);
     end
end

for m = 1 : taxadimension
     tmp{m} = names{m};
end

for m = 1 : (taxadimension - 1)
     tmp{taxadimension + m} = ['(' tmp{links(m, 1)} ', ' tmp{links(m, 2)} 
')'];
end

tree = tmp{2 * taxadimension - 1};


function [set1, set2] = xc_randomhalf(feat)

% FUNCTION [SET1 SET2] = XC_RANDOMHALF(FEAT)
% Take in a cell array of features (each element of the cell array is a class)
% and randomly devide it into half (for each element) and return as set1 and
% set2

rand('state',sum(100*clock));
for m = 1: length(feat)
     featsize = size(feat{m}, 1);
     r = randperm(featsize);
     set1{m} = feat{m}(r(1 : round(featsize / 2)) , :);
     set2{m} = feat{m}(r(round(featsize / 2) + 1 : featsize), :);
end


function [text, id, outarray] = xc_consensus2array(array, inputtree)
     treetext = inputtree;
     first = 1;
     tmp = [];
     o1 = array;
     while ((~isempty(treetext)) & (first | (treetext(1) == ',')))
         if (first)
             first = 0;
         else
             treetext(1) = [];
         end
         if (treetext(1) == '(')
             [treetext, id1, o1] = xc_consensus2array(o1, treetext(2:end));
         else
             comma = findstr(treetext, ',');
             comma2 = findstr(treetext, ')');
             if (isempty(comma) & isempty(comma2))
                 t = treetext;
                 i = 0;
             else if (isempty(comma))
                     t = treetext(1:comma2(1) -1 );
                     i = 2;
                 else if (isempty(comma2))
                         t = treetext(1:comma(1) - 1);
                         i = 1;
                     else if (comma(1) < comma2(1))
                             t = treetext(1:comma(1) - 1);
                             i = 1;
                         else
                             t = treetext(1:comma2(1) - 1);
                             i = 2;
                         end
                     end
                 end
             end
             
             id1 = str2num(t);
             %o1 = array;
             if (i == 0)
                 treetext = [];
             end
             if (i == 1)
                 treetext = treetext(comma(1):end);
             end
             if(i == 2)
                 treetext = treetext(comma2(1):end);
             end
         end
         tmp = [tmp id1];
     end
     id = length(o1) + 1;
     o1{id} = tmp;
     outarray = o1;
     if (length(treetext > 1))
         text = treetext(2:end);
     else
         text = [];
     end


