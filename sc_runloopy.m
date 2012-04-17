function belief = sc_runloopy(dnapot, edgepot, param)
%SC_RUNLOOPY performs the prior updating with voting potential
%   BELIEF = SC_RUNLOOPY(DNAPOT, EDGEPOT, PARAM)
%   
% Input:
%   DNAPOT is the dna potential, which is a 3D matrix with demensions 3 X row X col
%   EDGEPOT is the edge potential, which is a 4D matrix with dimensions row X col X W X W
%   PARAM is a struct which has the parameter settings for the prior updating with voting potential
%
% Output:
%   BELIEF is a 3D matrix with demensions 3 X row X col
%
% usage:
%   [belief, message] = sc_runloopy(message, dnapot, edgepot, param);
%
%   See also: sc_subregion
%
%   09-Aug-2006 Written by S. C. Chen

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

ncell = param.ncell;
nclass = ncell + 1;
alpha = param.alpha;  
beta = param.beta;
wndsz = param.wndsz;
row = param.row;
col = param.col;
stopthreshold = param.stopthreshold;
count = 0;

% message is the passed information between the neighboring pixels, which
% is a 5D matrix with demensions 3 X row X col X W X W.  It is intialized
% to be uniform across the classes and will the updated and considered as
% an updated prior
message = ones(nclass, row, col, wndsz, wndsz)/nclass;

while 1
    tic;
    count=count+1;
    
    % messages from the factor nodes to center node x_i 
    messagenew = message;   
    sigmaext = reshape(edgepot,[1 size(edgepot)]);
    sigmaext = repmat(sigmaext,nclass,1);
    mext = message.*(1-sigmaext)*alpha + (1-message) .* sigmaext*beta;       
    message_sum = 1+squeeze(sum(sum(mext,4),5));
            
    % messages from the center node node x_i to other factor nodes
    phy = message_sum./repmat(sum(message_sum,1),size(message_sum,1),1);
    phi = dnapot.*phy;    
    phi = phi./repmat(sum(phi,1),size(phi,1),1);

    for nrrel = -ncell:ncell
        for ncrel = -ncell:ncell
            if nrrel<=0
                mr = 1;
                mer = -nrrel;
                pr = 1-nrrel;
                per = 0;
            else
                mr = nrrel+1;
                mer = 0;
                pr = 1;
                per = nrrel;
            end
            
            if ncrel<=0
                mc = 1;
                mec = -ncrel;
                pc = 1-ncrel;
                pec = 0;
            else
                mc = ncrel+1;
                mec = 0;
                pc = 1;
                pec = ncrel;
            end
            message(:,mr:end-mer,mc:end-mec,-nrrel+nclass,-ncrel+nclass)...
                = phi(:,pr:end-per,pc:end-pec);
        end
    end    
    m = abs(message - messagenew);
    elapsed_time = toc;
    fprintf(1,'%2d iteration\telapsed_time:%.3f\tmessage difference:%.3f\n',count,elapsed_time, sum(m(:))); 
    if sum(m(:)) < stopthreshold
        break;
    else
        messagenew = message;
    end            
end

% calculate belief
belief= dnapot.* phy;
normFactor = sum(belief,1);
for k=1:size(belief,1)
    belief(k,:,:) = belief(k,:,:)./normFactor;
end
