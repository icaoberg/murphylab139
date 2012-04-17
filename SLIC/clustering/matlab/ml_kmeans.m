function g = ml_kmeans(feat, groupfile, distance, method, penc, varargin)

% FUNCTION G = ML_KMEANS(FEAT, GROUPFILE, DISTANCE, METHOD, PENC)
% ml_kmeans will perform k-means algorithm followed by AIC or 
% BIC selection to choose the optimal k.
% feat: the input feature values. It is a cell array where each element is
%       the feature matrix for one clone
% groupfile: filename to store the k-means result.  Use ml_readaicbic to read 
%            the information.
% distance: specify the distance function.  One of the following values: 
%           'euclidean' (by default), 'mahalanobis'.
% method: Either 'aic' or 'bic'.  Default 'aic'.
% penc: pencentage of the threshold.  If a cluster has the highest number of
%       observations from a clone and the pencentage is greater than the penc,
%       The clone will be considered to be in this cluster. 33 by default.
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
     error('ml_kmeans must have the argument: feat.  Type "help ml_kmeans" for more information');
  end
  if(~exist('groupfile', 'var'))
     groupfile = ['kmeanstmp' num2str(length(feat))];
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
    error ('Method must be one of the following ''euclidean'', ''mahlanobis'' or ''both''.  Type ''help ml_kmeans'' for details.');
end

  if (~exist('penc') | isempty(penc) | ~ isnumeric(penc))
     penc = 33;
  end

%  if (~exist([groupfile '.mat'], 'file'))
      ml_kmeans_aicbic(feat, groupfile, distance, varargin{:});
%  end
   [g, c] = ml_readaicbic(groupfile, distance, method, penc);
