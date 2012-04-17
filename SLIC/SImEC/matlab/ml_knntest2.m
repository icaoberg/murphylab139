function [pvalue,ts]=ml_knntest2(X1,X2,k,s)

%ML_KNNTEST2: Nearest neighbor two smaple test 
%   [PVALUE,TS]=ML_KNNTEST2(X1,X2,K,S) performs hypothesis test on X1 and X2
%   to test whethe they could have the same distribution. K is the number
%   of neariest neighbors and S is the distance function. For more options
%   of S, please refer to pdist.
%   [PVALUE,TS]=ML_KNNTEST2(X1,X2,K) uses the default distance function 'eu'
%   (Euclidean distance).
%   [PVALUE,TS]=ML_KNNTEST2(X1,X2) uses the defult K.
%   The funciton will return p-value PVALUE and test statistic TS.  
%   reference:
%   A multivariate two-sample test based on the number of nearest neighbor type coincidences
%                    Norbert Henze, 1988

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

if ~exist('s','var')
    s='eu';
end

if nargin<3
    k=-1;
end

if k==-1
    k=min([size(X1,1)-1,size(X2,1)-1,5]);
end

X=zscore([X1;X2]);

n1=size(X1,1);
n2=size(X2,1);
n=n1+n2;

pwdist=squareform(pdist(X,s,n1));
mindist=min(pwdist(pwdist~=0));
[index1,index2]=find(pwdist==0);
for i=1:length(index1)
    if index1(i)~=index2(i)
        pwdist(index1(i),index2(i))=unifrnd(0,mindist);
    end
end

knnmatrix=ml_rank(pwdist);
knnmatrix=knnmatrix-1;
B=[[ones(n1,n1),zeros(n1,n2)];[zeros(n2,n1),ones(n2,n2)]];
t=zeros(size(B));
for r=1:k
    a(:,:,r)=(knnmatrix==r);
    t=t+a(:,:,r).*B;
end

ts=sum(t(:));

ap=sum(a,3);

d=sum(ap,1);
c=sum((d-k).^2)/n/k;
v=sum(sum(ap.*ap'))/n/k;
m=(n1*(n1-1)+n2*(n2-1))/(n-1);
q=4*(n1-1)*(n2-1)/((n-2)*(n-3));
EL=k*m;
VarL=k*n1*n2*(q*(1+v-2*k/(n-1))+(1-q)*c)/(n-1);

pvalue=1-normcdf((ts-EL)/sqrt(VarL));

ts=ts/n/k;