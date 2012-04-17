function [pvalue,ts] = ml_ht2test2(X, Y, pooled)
%ML_HT2TEST2 Hotelling T2 test
%   PVALUE = ML_HT2TEST2(X,Y) performs the pooled Hotelling T2 test to
%   compare two feature matrices X and Y and returns p-value.
%   
%   PVALUE = ML_HT2TEST2(X,Y,POOLED) does unpooled Hotelling T2 test if
%   POOLED is 0 and pooled Hotelling T2 test if POOLED is not 0.
%   Pooled Hotelling T2 test: 
%       (mean(X)-mean(Y))*inv(Sp*sqrt(1/m+1/n))*(mean(X)-mean(Y))' ~ 
%           (m+n-2)*p/(m+n-p-1)*Fp,n+m-p-1
%        (Fp,n+m-p-1 means an F districtuion with the degrees of freedom p
%        and n+m-p-1)
%   Unpooled Hotelling T2 test: 
%       (mean(X)-mean(Y))*inv(Sx/m+Sy/n)*(mean(X)-mean(Y))' ~ X2p 
%       (X2p means chi square with the freedom p)
%   
%   [PVALUE,TS] = ML_HT2TEST2(...) also returns test statistic.
%   
%   See also

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

%   ??-???-200? Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU


if size(X,1)+size(Y,1)<=size(X,2)+1
    error(['Hotelling T2 test: The total number of samples minus ' ...
        ' the number of features must be greater than 1']);
end

%Initialize p-value
pvalue=1;

if ~exist('pooled','var')
    pooled=1;
end

%Find features that are constant
constidx=ml_constidx([X;Y]);;

if ~isempty(constidx)
    %If the constant features in X and Y have different values,
    %p-value is set to 0.
    warning('constant included'); 
    if any(X(1,constidx)~=Y(1,constidx))
        pvalue=0;
        ts=Inf;
    else
        %If the constant features in X and Y have the same value,
        %these features are removed.
        warning(['feature ' num2str(constidx) ' removed'])
        X(:,constidx)=[];
        Y(:,constidx)=[];
    end
end

%Calclulate p-value
if pvalue~=0
    if(pooled==0)
        %unpooled
        ts=ml_twomaha(X,Y,0,1);
        pvalue=1-chi2cdf(ts,size(X,2));
    else
        %pooled
        [ts,df1,df2]=ml_fvalue(X,Y);
        if isnan(ts)
            error('More samples may be needed');
        end
        pvalue=1-fcdf(ts,df1,df2);
    end
end