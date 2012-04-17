function both = ml_compfeats_u(feature_matrix1,feature_matrix2,method,issorted)

%ML_COMPFEATS_U Uivariate hypothesis test: compare two samples feature by feature.
%   BOTH=ML_COMPFEATS_U(FEATS1,FEATS2,METHOD,ISSORTED) performs 
%   a uivariate test for individual features of FEATS1 and FEATS2 by METHOD.
%   FEATS1 and FEATS2 are feature matrices should have the same number of columns.
%   METHOD available now are:
%   'ml_ttest2': t-test
%   'ml_kstest2': ks test
%   'ml_waldtest2': Wald test
%   'ml_pmttest2': permutation test
%   BOTH have two rows. The first row contains the indices of features, 
%   and the second row contains corresponding p-values. If ISSORTED is 1, which is 
%   its default value, p-values in BOTH will be sorted in an asending order.
%   BOTH=ML_COMPFEATS_U(FEATS1,FEATS2,METHOD) uses the default value of ISSORTED

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

if ( size(feature_matrix1, 2) ~= size(feature_matrix2, 2) )
    error('Matrices must have the same number of features');
end

allmethods={'ml_ttest2' 'ml_kstest2' 'ml_waldtest2'};

if all(strcmp(allmethods,method)==0)
    error('Unrecognized comparison method');
end

if ~exist('issorted','var')
    issorted=1;
end

both=1:size(feature_matrix1,2);
pvalues=[];

%Calculate p-values
for featindex=1:size(feature_matrix1,2)
    s1=feature_matrix1(:,featindex);
    s2=feature_matrix2(:,featindex);
    pvalue=eval(strcat(method,'(s1,s2)'));
    pvalues=[pvalues pvalue];
end

both=[both;pvalues];

%Sort p-values
if issorted==1
    conf = both(2,:);
    [sorted_conf index] = sort(conf');
    both = [index sorted_conf]';
end

