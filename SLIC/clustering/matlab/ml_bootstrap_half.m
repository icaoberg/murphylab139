function ml_bootstrap_half(filename, treename, feat, name, cycle, meth, ...
    covs, lmeth, comment)
% FUNCTION ML_BOOTSTRAP_HALF(FILENAME, TREENAME, FEAT, NAME, CYCLE, METH, COVS, LMETH, COMMENT)
%   Create a bootstrap analysis using randomly halfing of the samples. Output 
%   to a file for consensus analysis program
% FILENAME: output file
% TREENAME: name of the tree
% FEAT: a cell with n (num of taxa) elements, each containing m features
% NAME: a cell array containing the names of the taxa
% CYCLE: number of cycles in bootstrap analysis, 100 by default
% METH: 0 for euclidean (default) and 1 for mahn
% COVS: if meth == 1 and invcov is a matrix, it will be used as the
%      covariance matrix, if it is empty, it will be calaulated from input.
% LMETH: method for linkage function, 'average' by defalt, 'single',
%        'complete','centroid' and 'ward' are other choices
%        'additive' for additive tree building
% COMMENT: comment to be written into the output file

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

% Xiang Chen, 12/22/2004

if (length(name)~=length(feat)) 
     error('lengths of names and features arrays must match.')
end

if (~exist('meth', 'var') | isempty(meth))
    meth = 0;
end

if (~exist('lmeth', 'var') | isempty(lmeth))
     lmeth = 'average';
end
if (~isnumeric(meth) | (meth ~= 0))
    meth = 1;
end

if (~exist('comment', 'var'))
     comment = [];
end

for m = 1 : length(feat)
     feat{m} = double(feat{m});
end

fid = fopen(filename, 'w');

if (~exist('cycle', 'var') | isempty(cycle)) 
     cycle = 100;
end

features=[];
allf = [];
for i=1:length(feat)
    features(i,:)=mean(feat{i}, 1);
    allf = [allf
	     feat{i}];
end
zs=zscore(allf);
me = mean(allf);
sd = std(allf);
sd(find(sd==0)) = 1;
features = (features - repmat(me, [length(feat) 1])) ./ repmat(sd, [length(feat) 1]);
if (~exist('covs') & meth)
    for i = 1 : length(feat)
        feat{i} = (feat{i} - repmat(me, [size(feat{i}, 1) 1])) ./ repmat(sd, [size(feat{i}, 1) 1]);
    end
    covs = xc_cov(feat);
end



%Header of the output file
out = ['#NEXUS\n\n[\ntrees from\n\n' comment '\n]\n\n\n'];

%Taxa information
taxa = ['begin taxa;\ndimensions ntax = ' num2str(length(name)) ';\ntaxlabels'];
for m = 1 : length(name)
     %Remove '( , ) -' from the name
     %[s, t] = strtok(name{m}, '<');
     [s, t] = strtok(name{m}, '(');
     if (length(t))
	t(1) = []
     	%t = strtok(t, '>');
     	t = strtok(t, ')');
     	[t, t2] = strtok(t, ',');
     	while(~isempty(t2))
        	t2(1) = [];
		[t1, t2] = strtok(t2, ',');
        	t = [t '+' t1];
     	end	
     	%s = name{m}
     end;
     [s1, s2] = strtok(s, '-');
     if (isempty(s2))
         s = s1;
     else
         s2(1) = [];
         s = [s1 '_' s2];
     end
     name{m} = s;
     if (exist('t', 'var') & length(t))
        name{m} = [s '<' t '>'];
     end
     taxa = [taxa ' ' name{m}];
end
taxa = [taxa ';\nend;\n'];

%Original tree without bootstrap
if (meth == 0)
	dists = pdist(features,'euclid');
else
	dists = ml_pdist(features, 'mah', covs);
end

additive = strcmp(lmeth, 'additive');
if (additive)
     dists = squareform(dists);
     links = xc_neighborjoining(dists);
     out = [out taxa 'begin trees;\ntree fundamental_' treename ' = [&U][&W 1]'];
    
else
     links = linkage(dists, lmeth);
     out = [out taxa 'begin trees;\ntree fundamental_' treename ' = [&R][&W 1]'];
end

out = [out xc_link2tree(links, name) ';\nend;\n\n'];

%Bootstrapping
out = [out taxa 'begin trees;\n'];

%rand('state',sum(100*clock))
for m = 1 : cycle
    [set1, set2] = xc_randomhalf(feat);
    features=[];
    allf = [];
    for i=1:length(feat)
        features(i,:)=mean(set1{i}, 1);
        allf = [allf
	        set1{i}];
    end
    me = mean(allf);
    sd = std(allf);
    sd(find(sd==0)) = 1;
    features1 = (features - repmat(me, [length(feat) 1])) ./ repmat(sd, [length(feat) 1]);
    features=[];
    allf = [];
    for i=1:length(feat)
        features(i,:)=mean(set2{i}, 1);
        allf = [allf
	        set2{i}];
    end
    me = mean(allf);
    sd = std(allf);
    sd(find(sd==0)) = 1;
    features2 = (features - repmat(me, [length(feat) 1])) ./ repmat(sd, [length(feat) 1]);

    if (meth == 0)
    	dists1 = pdist(features1, 'euclid');
    	dists2 = pdist(features2, 'euclid');
    elseif (meth == 2)
        dists1 = ml_pdist(features1, 'city');
        dists2 = ml_pdist(features2, 'city');
    elseif (meth == 3)
        dists1 = ml_pdist(features1, 'ca');
        dists2 = ml_pdist(features2, 'ca');
    else
	if (~exist('covs', 'var')) 
	    dists1 = pdist(features1, 'mah');
	    dists2 = pdist(features2, 'mah');
	else
	    dists1 = ml_pdist(features1, 'mah', covs);
	    dists2 = ml_pdist(features2, 'mah', covs);
	end    
    end
    if (additive)
        links1 = xc_neighborjoining(squareform(dists1));
        links2 = xc_neighborjoining(squareform(dists2));
        out = [out 'tree B_1.' num2str(2 * m - 1) ' = [&U][&W 1]'];
        out = [out xc_link2tree(links1, name) ';\n'];
        out = [out 'tree B_1.' num2str(2 * m) ' = [&U][&W 1]'];
        out = [out xc_link2tree(links2, name) ';\n'];
    else
	links1 = linkage(dists1, lmeth);
        links2 = linkage(dists2, lmeth);
        out = [out 'tree B_1.' num2str(2 * m - 1) ' = [&R][&W 1]'];
        out = [out xc_link2tree(links1, name) ';\n'];
        out = [out 'tree B_1.' num2str(2 * m) ' = [&R][&W 1]'];
        out = [out xc_link2tree(links2, name) ';\n'];
    end
end

out = [out 'end;'];

fprintf(fid, out);
fclose(fid);

function cov_tmp = xc_cov(allf)
for m = 1 : length(allf)
    nocells(m) = size(allf{m}, 1);
end
s = min(nocells);
for m = 1 : 100
   tmpf = [];
   for n = 1 : length(allf) 
       x = randperm(size(allf{n}, 1));
       tmpf = [tmpf; allf{n}(x(1:s), :)];
   end
   cov_g(:,:,m) = cov(double(tmpf));
end
cov_tmp = mean(cov_g, 3);

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
     tmp{taxadimension + m} = ['(' tmp{links(m, 1)} ', ' tmp{links(m, 2)} ')'];
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

