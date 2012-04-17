function dist=ml_twomaha(X,Y,pooled,zscored)

%ML_TWOMAHA: Mahalanobisis distance between two populations
%   DIST=ML_TWOMAHA(X,Y,POOLED,ZSCORED) calculates 
%   the Mahalanobisis distance between two matrices X and Y.
%   If POOLED=0, unpooled distance will be calculated, otherwise
%   pooled distance will be calculated. If ZSCORED=0, there
%   is no preprocessing of the data. Otherwise, the pooled
%   data will be normalized to mean 0 and variance 1.

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

%   Written by T. Zhao
p=size(X,2);
feats=[X;Y];
m=size(X,1);
n=size(Y,1);

if(zscored~=0)
    feats=zscore(feats);
    X=feats(1:m,:);
    Y=feats((m+1):end,:);
end

if(pooled==0)
    S=cov(X)/m+cov(Y)/n; 
else
    S=(1/m+1/n)*((m-1)*cov(X)+(n-1)*cov(Y))/(m+n-2);
end

mux=mean(X);
muy=mean(Y);

dist=((mux-muy)/S)*(mux-muy)';

