function g = ml_consensus(feat, name, groupfile, distance, method, penc, niter)

% FUNCTION G = ML_CONSENSUS(FEAT, name, GROUPFILE, DISTANCE, METHOD, PENC, NITER)
% ml_consensus will perform k-means algorithm (if not done) followed by AIC or 
% BIC selection to choose the optimal k.  Then a set of dendrograms built on 
% randomly divided half set were created and used as input for majority
% consensus analysis. Finally a figure showing the consensus tree will be 
% displayed.
% feat: the input feature values. It is a cell array where each element is
%       the feature matrix for one clone
% name: name of each clone, a cell array with strings
% groupfile: filename to store the k-means result.  Use ml_readaicbic to read 
%            the information.
% distance: specify the distance function.  One of the following values: 
%           'euclidean' (by default), 'mahalanobis'.
% method: Either 'aic' or 'bic'.  Default 'aic'.
% penc: pencentage of the threshold.  If a cluster has the highest number of
%       observations from a clone and the pencentage is greater than the penc,
%       The clone will be considered to be in this cluster. 33 by default.
% niter: Number of iterations in creating half set dendrograms (100 by default,
%        which yields 200 dedrograms).
% g: the resturn membership array from k-means/AIC(BIC).  Each row presents a
%    clone and first k columns stands for the number of the cells in that
%    cluster.  The k+1 column is the highest percentage of cells for that clone
%    in any of the clusters.  The k+2 column is the ID of the cluster which has
%    that largest percentage, if it is not less than perc. 
% Xiang Chen, Jan 04, 2005

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

  if (~exist('feat', 'var'))
     error('ml_consensus must have the argument: feat.  Type "help ml_consensus" for more information');
  end
  if (~exist('name', 'var'))
     for m = 1 : length(feat);
         name{m} = num2str(m);
     end
  end
  if(~exist('groupfile', 'var'))
     groupfile = ['consensustmp' num2str(length(feat))];
  end

  currentdir = pwd;
  if (groupfile(1) ~= '/')
	groupfile = [currentdir '/' groupfile];
  end
  if (~exist('niter', 'var')|~isnumeric(niter))
     niter = 100;
  end
  if (~exist('distance', 'var'))
     distance = 'eu';
  end

  if (exist('method', 'var') & strcmp(method, 'bic'))
     bic = 1;
  else
     bic = 0;
     method = 'aic';
  end

  switch distance(1:2)
   case 'eu'
    eu = 1;
    ma = 0;
   case 'ma'
    ma = 1;
    eu = 0;
   otherwise
    error ('Method must be one of the following ''euclidean'', ''mahlanobis'' or ''both''.  Type ''help ml_consensus'' for details.');
end

  if (~exist('penc') | ~ isnumeric(penc))
     penc = 33;
  end

  if (~exist([groupfile '.mat'], 'file'))
      ml_kmeans_aicbic(feat, groupfile, distance);
  end
  [g, c] = ml_readaicbic(groupfile, distance, method, penc);
  features = [];
  for m = 1 : length(feat)
          if (c{m})
             features{length(features)+1} = double(feat{m}(c{m}, :));
	     gene{length(features)} = name{m};
	  end
  end

  %if (~exist([groupfile '.den'], 'file'))
      ml_bootstrap_half([groupfile '.den'], '', features, gene, niter, ma);
      %end
 % if (~exist([groupfile '.den_MRC'], 'file'))
      javadir = which('ml_consensus');
      javadir = [javadir(1:end-21) 'bin'];
      cd(javadir);
      unix(['java MajConsensus ' groupfile '.den']);
      cd(currentdir);
      %end
  [links, outarr] = ml_readconsensus([groupfile '.den_MRC'], features, ma);
  eval('ml_drawconsensus(links, gene, 9, ma);', ...
       'ml_drawconsensus(links(1:end-1), gene, 9, ma);');
  
