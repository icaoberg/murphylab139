function g = ml_kmeans_pslid(feat, groupfile, distance, method)

% FUNCTION G = ML_KMEANS_PSLID(FEAT, GROUPFILE, DISTANCE, METHOD)
% ml_kmeans_pslid will perform k-means algorithm (if not done) followed by 
% AIC or BIC selection to choose the optimal k.
% feat: the input feature values. It is a feature matrix, each column is a
%       feature and each row is an image
% groupfile: filename to store the cluster result for each k
% distance: specify the distance function.  One of the following values: 
%           'euclidean' (by default), 'mahalanobis'.
% method: Either 'aic' or 'bic'.  Default 'aic'.
% g: the resturn membership array from k-means/AIC(BIC).  Each element
%    represents the ID of the cluster which contains the corresponding image
% Xiang Chen, Aug 19, 2005

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
     error('ml_kmeans_pslid must have the argument: feat.  Type "help ml_kmeans_pslid" for more information');
  end
  if(~exist('groupfile', 'var'))
     htmlfile = ['kmeanstmp' num2str(size(feat, 1)) '.mat'];
  end
  currentdir = pwd;
  if (htmlfile(1) ~= '/')
	htmlfile = [currentdir '/' groupfile];
  end
  dots = findstr(groupfile, '.');
  if (isempty(dots) | ~strcmp(htmlfile(dots(end):end), '.mat'))
     groupfile = [groupfile '.html'];
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
    error ('Method must be one of the following ''euclidean'', ''mahlanobis'' or ''both''.  Type ''help ml_kmeans_pslid'' for details.');
  end

  feat = zscore(feat);
  cluster_m = [];
  min = inf;
  done = 0;
  cnt = 1;
  if (ma)
     cov_f = cov(feat);
     while (~done)
        cluster = ml_kmeans_mah(feat, [], cnt,cov_f);
	[a, b] = ml_aicbic_mah(feat, cluster, cov_f);
	if (bic)
            diff = b - min;
	else
   	    diff = a - min;
            b = a;
        end
	if (diff < 0)
            min = b;
	    cluster_m = cluster;
	elseif (diff >= (min * 0.2))
            done = 1;
	end
     end
  else
     while (~done)
        cluster = ml_kmeans_euc(feat, [], cnt);
	[a, b] = ml_aicbic_euc(feat, cluster, cov_f);
	if (bic)
            diff = b - min;
	else
   	    diff = a - min;
            b = a;
        end
	if (diff < 0)
            min = b;
	    cluster_m = cluster;
	elseif (diff >= (min * 0.2))
            done = 1;
	end
     end
  end

  g = zeros(1, size(feat, 1));
  for m = 1 : length(cluster_m)
     g(cluster_m{m}) = m;
  end
         
