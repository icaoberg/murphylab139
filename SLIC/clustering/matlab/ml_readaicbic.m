function [groupinfo, cluster, groupinfo2, cluster2] = ...
    ml_readaicbic(filename, distance, method, penc)
% FUNCTION [GROUPINFO, CLUSTER, GROUPINFO2, cluster2] = ...
%               ML_READAICBIC(FILENAME, DISTANCE, METHOD, PENC)
% find the optimal clustering result based on AIC or BIC.  The input to this
% function should be the result file from ml_kmeans_aicbic.m
% filename: the input file created by ml_kmeans_aicbic
% distance: one of the following: 'euclidean', 'mahlanobis', or 'both'.  The
%           specified distance must have already been calculated in
%           ml_kmeans_aicbic.m.  By default, 'euclidean'.
% method: Either 'aic' or 'bic'.  Default 'aic'.
% penc: pencentage of the threshold.  If a cluster has the highest number of
%       observations from a clone and the pencentage is greater than the penc,
%       The clone will be considered to be in this cluster. 33 by default.
% groupinfo: a m * (n + 2) matrix, for m clones and n clusters.  The first n
%            columns stands for the number of observations in each cluster.
%            The next column is the highest percentage of the observations in
%            each cluster.  The last column is the cluster ID with the highest
%            pencentage (0 if it is less than penc specified)
% cluster: a m-element cell array, each element contains the index of the
%          observations from the clone that are clustered into the cluster with
%          highest pencentage (empty if less than penc).
% groupinfo2, cluster2: same format as groupinfo and cluster, only used when
%                       the distance is 'both'.  In this case, the groupinfo 
%                       and cluster will be for Euclidean distance and the
%                       groupinfo2 and cluster2 represent Mahalanobis distance
% Xiang Chen, 12/22/2004

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

if (~exist('filename'))
    error('Not valid input!');
end

if (~exist('distance', 'var') | isempty(distance))
    distance = 'eu';
end

if (exist('method', 'var') & strcmp(method, 'bic'))
    bic = 1;
else
    bic = 0;
end

switch distance(1:2)
case 'eu'
    eu = 1;
    ma = 0;
case 'ma'
    ma = 1;
    eu = 0;
case 'bo'
    eu = 1;
    ma = 1;
otherwise
    error ('Method must be one of the following ''euclidean'', ''mahlanobis'' or ''both''.  Type ''help ml_kmeans_aic'' for details.');
end

if (~exist('penc') | ~ isnumeric(penc))
    penc = 33;
end

load(filename);
if (ma)
    if (bic)
        [Y, I] = sort(bic_mah);
    else
        [Y, I] = sort(aic_mah);
    end
    
    for m = 1 : I(1)+kmin-1
        tmp = floor(group_mah{I(1)}{m});
        u = unique(tmp);
        for o = 1 : length(u)
            groupinfo(u(o), m) = length(find(tmp == u(o)));
        end
    end
    [V, I2] = max(groupinfo, [], 2);
    groupinfo(:, m + 1) = uint8(round(V ./sum(groupinfo, 2) * 100));
    groupinfo(:, m + 2) = uint8(I2 .* (groupinfo(:, m + 1) >=penc));
    for o = 1 : size(groupinfo, 1)
        if groupinfo(o, m + 2)
            tmp = group_mah{I(1)}{groupinfo(o, m + 2)};
            idx = find(floor(tmp) == o);
            cluster{o} = uint8(round((tmp(idx) - o) * 1000));
        end
    end
    if (eu)
        groupinfo2 = groupinfo;
        cluster2 = cluster;
        clear groupinfo cluster
    end
end

if (eu)
    if (bic)
        [Y, I] = sort(bic_euclid);
    else
        [Y, I] = sort(aic_euclid);
    end
    
    for m = 1 : I(1)+kmin-1
        tmp = floor(group_euclid{I(1)}{m});
        u = unique(tmp);
        for o = 1 : length(u)
            groupinfo(u(o), m) = length(find(tmp == u(o)));
        end
    end
    [V, I2] = max(groupinfo, [], 2);
    groupinfo(:, m + 1) = uint8(round(V ./sum(groupinfo, 2) * 100));
    groupinfo(:, m + 2) = uint8(I2 .* (groupinfo(:, m + 1) >= penc));
    for o = 1 : size(groupinfo, 1)
        if groupinfo(o, m + 2)
            tmp = group_euclid{I(1)}{groupinfo(o, m + 2)};
            idx = find(floor(tmp) == o);
            cluster{o} = uint8(round((tmp(idx) - o) * 1000));
        end
    end
    
end
