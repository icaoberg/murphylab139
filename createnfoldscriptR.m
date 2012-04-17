function createnfoldscriptR(flag, classnum, randseed, allflag, maximgperfold)
%CREATENFOLDSCRIPTR creates the index for training images, training features, and sda idx for all folds
%   CREATENFOLDSCRIPTR(FLAG, CLASSNUM, RANDSEED, ALLFLAG, MAXIMGPERFOLD)
%   
% Input: 
%   FLAG specifiy which feature matrix to use (either with unique location, ambiguous or punctate_composite class)
%   CLASSNUM specify the number of class to classify
%   RANDSEED specify the rand seed to split the training/testing data
%   ALLFLAG specify if we want to use all the images in the fold for training (if 1, maximgperfold will not be useful)
%   MAXIMGPERFOLD specifiy the maximum number of images to for trainin at each fold
%
%   See also createnfoldrandperm, sc_makegfm
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

for idxfold=1:param.fold
    param.idxfold = idxfold;
    createnfoldrandperm(param);
end
