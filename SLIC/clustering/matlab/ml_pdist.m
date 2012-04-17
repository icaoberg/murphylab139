function Y = ml_pdist(X,s,t)
%ML_PDIST Pairwise distance between observations.  Modified from matlab 
version
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

%of PDIST function.
%   Y = ML_PDIST(X,METRIC) returns a vector which contains all the
%   distances between each pair of observations in X computed using
%   the given METRIC.  X is a M by N matrix, treated as M observations
%   of N variables. Since there are M*(M-1)/2 pairs of observations in
%   X, the size of Y is M*(M-1)/2 by 1.  The default metric is
%   'EUCLID'.  The available metrics are:
%
%      'euclid'    --- Euclidean metric
%      'seuclid'   --- Standardized Euclid metric
%      'cityblock' --- City Block metric
%      'mahal'     --- Mahalanobis metric
%      'minkowski' --- Minkowski metric
%
%   Y = ML_PDIST(X, 'minkowski', p) specifies the exponents in the
%   'Minkowski' computation. When p is not given, p is set to 2.
%
%   Y = ML_PDIST(X, 'mah', p) specifies the covariance matrix in the
%   'Minkowski' computation. When p is not given, it will be estimated 
%   from X.
%
%   The output Y is arranged in the order of ((1,2),(1,3),..., (1,M),
%   (2,3),...(2,M),.....(M-1,M)).  i.e. the upper right triangle of
%   the M by M square matrix. To get the distance between observation
%   i and observation j, either use the formula Y((i-1)*(M-i/2)+j-i)
%   or use the helper function Z = SQUAREFORM(Y), which will return a
%   M by M symmetric square matrix, with Z(i,j) equaling the distance
%   between observation i and observation j.
%
%   See also SQUAREFORM, LINKAGE

%   ZP You, 3-10-98
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 1.2 $
%   Modified Mahalanobis distance calculation for the option of getting a
%   covariance matrix passed in.
%   Xiang Chen, 12/7/2004

if nargin >= 2
   if length(s) < 2
      error('Unrecognized metric');
   else 
      s = lower(s(1:2));
   end
else
   s = 'eu';
end

if s == 'mi' % Minkowski distance need a third argument
   if nargin < 3
      t = 2; 
   elseif t <= 0
      error('The third argument has to be positive.');
   end
end

[m, n] = size(X);

if m < 2
   error('The first argument has to be a numerical matrix with at least two rows');
end

p = (m-1):-1:2;
I = zeros(m*(m-1)/2,1);
I(cumsum([1 p])) = 1;
I = cumsum(I);
J = ones(m*(m-1)/2,1);
J(cumsum(p)+1) = 2-p;
J(1)=2;
J = cumsum(J);

if (strcmp(s, 'ca'))
     Y1 = X(I,:)';
     Y2 = X(J,:)';
     clear I J p
else
     Y = (X(I,:)-X(J,:))';
     I = []; J = []; p = [];  % no need for I J p any more.
end

switch s
case 'eu' % Euclidean
   Y = sum(Y.^2);
   Y = sqrt(Y);
case 'se' % Standadized Euclidean
   D = diag(var(X));
   Y = sum(D*(Y.^2));
   Y = sqrt(Y);
case 'ci' % City Block
   Y = sum(abs(Y));
case 'ma' % Mahalanobis
   if (exist('t', 'var'))
         v = inv(t);
   else
         v = inv(cov(X));
   end
   Y = sqrt(sum((v*Y).*Y));
case 'mi' % Minkowski
   Y = sum(abs(Y).^t).^(1/t);
case 'ch' % Chebychov
   Y = max(abs(Y));
case 'ca' %CosAngle
   Y = acos(dot(Y1, Y2) ./ sqrt(dot(Y1, Y1) .* dot(Y2, Y2)));
otherwise
   error('no such method.');
end
