function [pvalue,ts]=ml_lvntest2(X,Y)

%ML_LVNTEST2: Levene's test for variances
%   [PVALUE,TS]=ML_KSTEST2 performs a Levene's test 
%   to determine if independent random multvariate samples,
%   X and Y, have the same variance.
%   The funciton will return p-value PVALUE and test statistic TS. 

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

X=zscore(X);
m=size(X,1);
Y=zscore(Y);
n=size(Y,1);

mux=mean(X);
muy=mean(Y);

X=abs(X-repmat(mux,m,1));
Y=abs(Y-repmat(muy,n,1));

[pvalue,ts]=ml_ht2test2(X,Y);