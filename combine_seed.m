function comebine_seed(param)
%COMBINE_SEED combines the result of different seeds (the result has to be combined by fold before) 
%   COMBINE_SEED(PARAM)
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
%             flag, classnum, nfold, allflag, maximgperfold, seedidx, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);
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


tic
addpath(genpath([pwd filesep 'SLIC']));

flag = param.flag;
classnum = param.classnum;
nfold = param.fold;

%idxfold = param.idxfold;
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

    % show the accuracy at i th seedidx
    cellconf{seedidx} = zeros(classnum,classnum);
    cellconfvoting{seedidx} = zeros(classnum,classnum);
    imgconfvoting{seedidx} = zeros(classnum,classnum);

    % load giant feature matrix file
    fprintf(1,' Loading classifier [%s] at %d th random seed...\n',combinefoldfilestr, seedidx);
    tmp1=load(combinefoldfilestr);
   
    [imagenums,order]=sort(yeastdata.ucsfresult.classSize, 'descend');
    testimgperclass = imagenums(1:classnum);
    testorder = order(1:classnum);    
 
    for idxfold = 1:nfold               
        fprintf(1,'    Loading rand perm index file and sda idx at %dth fold...\n',idxfold);
        nfoldrandpermfilestr = sprintf('celllevel/result/PERM_F%dC%dN%dI%dA%dM%dR%d.mat', ...
       	      flag, classnum, nfold, idxfold, allflag, maximgperfold, randseed);
        tmp2 = load(nfoldrandpermfilestr);
        dataidx = tmp2.dataidx;
        dataidxarray{seedidx} = dataidx;
        clear tmp2;        
        for t=1:classnum     
            testimgnum = length(dataidx{t}{idxfold});
            for e=1:testimgnum
                yeastdata.ucsfresult.imganot{testorder(t)}{dataidx{t}{idxfold}(e)}.correctcell(seedidx) = ...
                    tmp1.yeastdata.ucsfresult.imganot{testorder(t)}{dataidx{t}{idxfold}(e)}.correctcell;
                yeastdata.ucsfresult.imganot{testorder(t)}{dataidx{t}{idxfold}(e)}.voteclass(seedidx) = ...
                    tmp1.yeastdata.ucsfresult.imganot{testorder(t)}{dataidx{t}{idxfold}(e)}.voteclass;
                yeastdata.ucsfresult.imganot{testorder(t)}{dataidx{t}{idxfold}(e)}.votecorrect(seedidx) = ...
		    tmp1.yeastdata.ucsfresult.imganot{testorder(t)}{dataidx{t}{idxfold}(e)}.votecorrect;

                if ~iscell(yeastdata.ucsfresult.cellanot{testorder(t)}{dataidx{t}{idxfold}(e)}.classify(1))
                    yeastdata.ucsfresult.cellanot{testorder(t)}{dataidx{t}{idxfold}(e)}.classify = [];
                end
                yeastdata.ucsfresult.cellanot{testorder(t)}{dataidx{t}{idxfold}(e)}.classify{seedidx} = ...
		    tmp1.yeastdata.ucsfresult.cellanot{testorder(t)}{dataidx{t}{idxfold}(e)}.classify;

                thisrow = zeros(1,classnum);

                for t2=1:classnum
                    thisrow(t2) = sum(tmp1.yeastdata.ucsfresult.cellanot{testorder(t)}{dataidx{t}{idxfold}(e)}.classify==t2);
                end
	       
                cellconf{seedidx}(t,:) = cellconf{seedidx}(t,:) + thisrow;
                [votecount,voteclass] = max(thisrow);
                cellconfvoting{seedidx}(t,voteclass) = cellconfvoting{seedidx}(t,voteclass) + sum(thisrow);
                imgconfvoting{seedidx}(t,voteclass) = imgconfvoting{seedidx}(t,voteclass) + 1;
            end
        end
    end

    imgconfvotingnorm{seedidx} = imgconfvoting{seedidx}./(sum(imgconfvoting{seedidx},2)*ones(1,classnum));
    fprintf(1,'image accuracy (u)(%.2f%%)(n)(%.2f%%)\n',100*sum(diag(imgconfvoting{seedidx}))/sum(imgconfvoting{seedidx}(:)), ...
                                                        100*sum(diag(imgconfvotingnorm{seedidx}))/sum(imgconfvotingnorm{seedidx}(:)));
    acc(seedidx) = 100*sum(diag(imgconfvoting{seedidx}))/sum(imgconfvoting{seedidx}(:));
    normacc(seedidx) = 100*sum(diag(imgconfvotingnorm{seedidx}))/sum(imgconfvotingnorm{seedidx}(:)); 
end       

elapsetime = toc

clear tmp tmp1;
savefilestr = sprintf('celllevel/result/SEED_F%dC%dN%dA%dM%dR%d_S%dW%d_C%dG%dN%d%s.mat', ...
		      flag, classnum, nfold, allflag, maximgperfold, randseed, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);


imgconfvotingseed = zeros(classnum,classnum);
imgconfvotingseedcell = zeros(classnum,classnum);

for seedidx = 1:randseed % up to how many seed
    imgconfvotingseed = zeros(classnum,classnum);
    for idxfold = 1:nfold
        for t=1:classnum
            testimgnum = length(dataidxarray{1}{t}{idxfold});
            for e=1:testimgnum
                thisrow = zeros(1,classnum);
                for t2=1:classnum
                    thisrow(t2) = sum(yeastdata.ucsfresult.imganot{testorder(t)}{dataidxarray{1}{t}{idxfold}(e)}.voteclass(1:seedidx)==t2);
                end
                [votecount,voteclass] = max(thisrow);

                % vote for each image
                imgconfvotingseed(t,voteclass) = imgconfvotingseed(t,voteclass) + 1;

            end
        end
    end

    imgconfvotingnormseed = imgconfvotingseed./(sum(imgconfvotingseed,2)*ones(1,classnum));

fprintf(1,'[%d seed]: image accuracy (u)(%.2f%%)(n)(%.2f%%)\n',seedidx, 100*sum(diag(imgconfvotingseed))/sum(imgconfvotingseed(:)), ...
        100*sum(diag(imgconfvotingnormseed))/sum(imgconfvotingnormseed(:)))

     accvote(seedidx) = 100*sum(diag(imgconfvotingseed))/sum(imgconfvotingseed(:));  % accumulate
end

yeastdata.ucsfresult.cellfeat = [];


elapsetime = toc
save(savefilestr);

[acc; normacc]
mean(acc)
mean(normacc)
accvote
