function [pvalues,ts]=ml_compfeats_m_pw(features,method,t)

%ML_COMPFEATS_M Pairwise multivariate hypothesis test
%   [PVALUES,TS]=ML_COMPFEATS_M_PW(FEATURES,METHOD,T) performs a pariwise multivariate
%   test to compare all pairs in FEATURES, which is a cell array of feature matrices.
%   See ML_COMPFEATS_M for details of the parameters METHOD and T.
%   The jth row and ith column (j>i) in PALUES and TS is the pvalue and test statistic 
%   from the jth and ith feature matrices.

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

if nargin<3
    t={};
end

prot_no=length(features);

for i=1:prot_no
    for j=(i+1):prot_no
        [pvalues(j,i),ts(j,i)]=ml_compfeats_m(features{i},features{j},method,t);  
    end
end
