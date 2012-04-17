function [pvalue,ts]=ml_wwtest2(X1,X2,s)

%ML_WWTEST2: Multivariate Wald-Wolfowitz test (Alos known as FR test)
%   [PVALUE,TS]=ML_WWTEST2(X1,X2,S) performs the FR test 
%   on two samples X1 and X2. Currently only Euclidean function is
%   supported, therefore S must be 'eu'.
%   The funciton will return p-value PVALUE and test statistic TS.
%   Reference:
%   Jerome H. Friedman and Lawrence C. Rafsky, 1979

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
    
X=zscore([X1;X2]);

m=size(X1,1);
n=size(X2,1);
N=m+n;
v=ml_mstree(X,s,m);

R=sum((v(1,:)<=m) & (v(2,:)>m))+sum((v(1,:)>m) & (v(2,:)<=m))+1;
Er=2*m*n/N+1;
C=0;

for i=1:N
    deg=sum(v(1,:)==i)+sum(v(2,:)==i);
    C=C+deg*(deg-1)/2;
end

Vr=2*m*n*((2*m*n-N)/N+(C-N+2)*(N*(N-1)-4*m*n+2)/(N-2)/(N-3))/N/(N-1);
W=(R-Er)/sqrt(Vr);
ts=R;
pvalue=normcdf(W,0,1);
