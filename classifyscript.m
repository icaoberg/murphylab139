function classifyscript(flag, classnum, randseed, allflag, maximgperfold, idxfold, sdaflag, weightflag, cidx, gidx, nbrs, classifierid)     
%CLASSIFYSCRIPT runs cell level classification on unique location classes with purality voting scheme at idxfold-th fold
%   CLASSIFYSCRIPT(FLAG, CLASSNUM, RANDSEED, ALLFLAG, MAXIMGPERFOLD, IDXFOLD, SDAFLAG, WEIGHTING, CIDX, GIDX, NBRS, CLASSIFIERID)
%
% Input:
%   FLAG specifiy which feature matrix to use
%   CLASSNUM specify the number of class to classify
%   RANDSEED specify the rand seed to split the training/testing data
%   ALLFLAG specify if we want to use all the images in the fold for training (if 1, maximgperfold will not be useful)
%   MAXIMGPERFOLD specifiy the maximum number of images to for trainin at each fold
%   IDXFOLD specify the n-th fold to use
%   SDAFLAG specifiy if sda idx is used
%   WEIGHTFLAG specifiy if the class-dependent weighting is used
%   CIDX specifiy the parameter index used for the slack penalty in LIBSVM
%   GIDX specifiy the parameter index used for the width of the gaussian kernel
%   NBRS specifiy the number of neighbors in KNN
%   CLASSIFIERID specifiy which type of classifier is used (1 for LIBSVM and 2 for KNN)
%
% Output:
%   savefilestr = sprintf('celllevel/result/RES_F%dC%dN%dI%dA%dM%dR%d_S%dW%d_C%dG%dN%d%s.mat', ...^M
%             flag, classnum, nfold, idxfold, allflag, maximgperfold, randseed, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);
%
%   See also sc_classifyyeast, createnfoldscriptR
%
%   06-June-2007 Written by S. C. Chen


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


addpath(genpath([pwd filesep 'SLIC']));

param.fold = 6;
param.flag = flag;
param.classnum = classnum;
param.allflag = allflag;
param.maximgperfold = maximgperfold;
param.randseed = randseed;

param.sdaflag = sdaflag;
param.weightflag = weightflag;
param.cidx = cidx;
param.gidx = gidx;
param.nbrs = nbrs;
param.classifierid = classifierid;

param.idxfold = idxfold;
sc_classifyyeast(param);
