function comebine_foldseed(param)
%COMBINE_FOLDSEED combine results from each fold and each seed, only used for ambiguous and punctate_composite class
%   COMBINE_FOLDSEED(PARAM)
%
% Input:
%   PARAM.FLAG specifiy which feature matrix to use
%   PARAM.CLASSNUM specify the number of class to classify
%   PARAM.NFOLD specify the number of fold to use for cross-validation
%   PARAM.IDXFOLD specify the n-th fold to use
%   PARAM.RANDSEED specify the rand seed to split the training/testing data
%   PARAM.ALLFLAG specify if we want to use all the images in the fold for training (if 1, maximgperfold will not be useful)
%   PARAM.MAXIMGPERFOLD specifiy the maximum number of images to for trainin at each fold
%   PARAM.SDAFLAG specifiy if sda idx is used
%   PARAM.WEIGHTFLAG specifiy if the class-dependent weighting is used
%   PARAM.CIDX specifiy the parameter index used for the slack penalty in LIBSVM
%   PARAM.GIDX specifiy the parameter index used for the width of the gaussian kernel
%   PARAM.NBRS specifiy the number of neighbors in KNN
%   PARAM.CLASSIFIERID specifiy which type of classifier is used (1 for LIBSVM and 2 for KNN)
%						     
% Output:
%   savefilestr = sprintf('celllevel/result/SEED_F%dC%dN%dA%dM%dR%d_S%dW%d_C%dG%dN%d%s.mat', ...
%	      flag, classnum, nfold, allflag, maximgperfold, seedidx, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);
%
%   See also sc_classifyyeast, combineseedscript, combinescript, combine_fold, combine_fold_all
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


tic;

flag = param.flag;
classnum = param.classnum;
nfold = param.fold;

allflag = param.allflag;
maximgperfold = param.maximgperfold;
randseed = param.randseed;
sdaflag = param.sdaflag;
weightflag = param.weightflag;

cidx = param.cidx;
gidx = param.gidx;
nbrs = param.nbrs;

if param.classifierid == 1
    classifiername = 'libsvm';
else
    classifiername = 'nn';
end

tmp =load([pwd filesep 'data' filesep 'yeastfeat_d' num2str(flag) '.mat']);
yeastdata = tmp.yeastdata;

for seedidx = 1:randseed % up to how many seed
    combinefoldfilestr = sprintf('celllevel/result/COMBINE_F%dC%dN%dA%dM%dR%d_S%dW%d_C%dG%dN%d%s.mat', ...
		      flag, classnum, nfold, allflag, maximgperfold, seedidx, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);

    oneclassnum = 1;
    % show the accuracy at i th seedidx
    imgconfvoting{seedidx} = zeros(oneclassnum,classnum);

    % load giant feature matrix file
    fprintf(1,' Loading classifier [%s] at %d th random seed...\n',combinefoldfilestr, seedidx);
    tmp1=load(combinefoldfilestr);
    
    for t=1:oneclassnum               
        for e=1:yeastdata.ucsfresult.classSize
            yeastdata.ucsfresult.imganot{t}{e}.voteclass = ...
                    [yeastdata.ucsfresult.imganot{t}{e}.voteclass ...
                     tmp1.yeastdata.ucsfresult.imganot{t}{e}.voteclass];             
            thisrow = zeros(1,classnum);
            for t2=1:classnum
                thisrow(t2) = sum(tmp1.yeastdata.ucsfresult.imganot{t}{e}.voteclass==t2);
            end
            [votecount,voteclass] = max(thisrow);             
            imgconfvoting{seedidx}(t,voteclass) = imgconfvoting{seedidx}(t,voteclass) + 1;
        end
    end
    imgconfvotingnorm{seedidx} = imgconfvoting{seedidx}./(sum(imgconfvoting{seedidx},2)*ones(1,classnum));
end      

elapsetime = toc


YFL034Widx = find(tmp1.yeastdata.ucsfresult.classidx{1} == strmatch('YFL034W', yeastdata.geneNames));

clear tmp tmp1;
savefilestr = sprintf('celllevel/result/SEED_F%dC%dN%dA%dM%dR%d_S%dW%d_C%dG%dN%d%s.mat', ...
		      flag, classnum, nfold, allflag, maximgperfold, randseed, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);

savefilestr
save(savefilestr);

imgconfvotingseed = zeros(classnum,classnum);

for seedidx = 1:randseed % up to how many seed
    imgconfvotingseed = zeros(oneclassnum,classnum);   
    for t=1:oneclassnum            
        for e=1:yeastdata.ucsfresult.classSize
            thisrow = zeros(1,classnum);                         
            for t2=1:classnum
                thisrow(t2) = sum(yeastdata.ucsfresult.imganot{t}{e}.voteclass(1:seedidx*nfold)==t2);
            end
            [votecount,voteclass] = max(thisrow);

            % vote for each image
            imgconfvotingseed(t,voteclass) = imgconfvotingseed(t,voteclass) + 1;
            
        end
    end
    imgconfvotingnormseed = imgconfvotingseed./(sum(imgconfvotingseed,2)*ones(1,classnum));
end


elapsetime = toc
yeastdata.ucsfresult.cellfeat = [];
save(savefilestr);

imgconfvotingnormseed

%keyboard;
