function [pvalue,ts]=ml_compfeats_m(feats1,feats2,method,t)

%ML_COMPFEATS_M Multivariate hypothesis test: compare whether two matrices
%   are statistically the same.
%   [PVALUE,TS]=ML_COMPFEATS_M(X1,X2,METHOD,T) performs a multivariate
%   test to determine whether the two samples X1 and X2 could be from the
%   same population by METHOD with additonal parameters T.
%   [PVALUE,TS]=ML_COMPFEATS_M(X1,X2,METHOD) uses the default parameters.
%   X1 and X2 are feature matrices that should have the same number of columns.
%   METHOD available now are:
%   'ml_ht2test2': Hotelling T2 test, with parameters T={0} (unpooled) or 
%       {1} (pooled)
%   'ml_wwtest2': Multivaraite Wald-Woldfowitz test (or FR test) with paramters
%       for distance function. For example T = {'eu'} is for using Euclidean
%       funtion.
%   'ml_knntest2': Nearest neighbor test, with parameters T{1} neighbors 
%       and distance function T{2}. For example, T = {5,'eu'} is for 5
%       neighbors and Euclidean fucntion.
%   PVALUE is the p-value. Emprically, PVALUE<0.05 indicates strong evidence
%   of that the two matrices are from different population. TS is the test statistic.

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

if nargin < 3
    error('Two matrices and one test method are required for a comparison');
end

if ( size(feats1, 2) ~= size(feats2, 2) )
    error('Matrices must have the same number of features');
end


allmethods={'ml_ht2test2' 'ml_wwtest2' 'ml_knntest2'};

if all(strcmp(allmethods,method)==0)
    error('Unrecognized comparison method');
end

if nargin<4
    t={};
end

%Check and remove rows with NaN values
nansample1=find(isnan(sum(feats1,2)));
nansample2=find(isnan(sum(feats2,2)));

if ~isempty(nansample1)
    warning(['Row ' num2str(nansample1') ' in the 1st group has NaN features. Removed.']);
    feats1(nansample1,:)=[];
end


if ~isempty(nansample2)
    warning(['Row ' num2str(nansample2') ' in the 2nd group has NaN features. Removed.']);
    feats2(nansample2,:)=[];
end


param={feats1,feats2,t{:}};
[pvalue,ts]=eval([method '(param{:})']);
