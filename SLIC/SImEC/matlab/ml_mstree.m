function v=ml_mstree(X,s,t)

%ML_MSTREE: build a minimal spanning tree
%   V=ML_MSTREE(X,S,T) build the minimal panning tree for multivariate data X.
%   S the distance funcion and T is the additional parameter for distance
%   function. See ml_pdist for details.
%   The minimal spanning tree is stored in V as [node1;node2;dist]

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

if size(X,1)<=1
    v=[];
    return;
end

%Calculate pairwise distance
G=squareform(pdist(X,s,t));

%Keep the lower triangular part of G
trilG=triu(NaN*ones(size(G)))+G;
symG=tril(G)+tril(G)';

pnum=size(G,1);
m=1:(pnum-1);
n=m;
c=m;

%Initialize
pset=1:pnum;
pnewset=[];

%The first node
[c(1),m(1),n(1)]=ml_min(trilG);
pnewset=[pnewset,m(1),n(1)];
pset([m(1),n(1)])=[];
pnewset=sort(pnewset);
pset=sort(pset);

%Keep adding nodes
for i=2:(pnum-1)
    m1=symG(pnewset,:);
    m1(:,pnewset)=[];
    [c(i),mm,nn]=ml_min(m1);
    m(i)=pnewset(mm);
    n(i)=pset(nn);
    pnewset=[pnewset,n(i)];
    pset(nn)=[];   
    pnewset=sort(pnewset);
    pset=sort(pset);
end

v=[m;n;c];