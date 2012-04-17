function comebine_fold(param)
%COMBINE_FOLD classify the images in the test fold and combine results for each fold
%   COMBINE_FOLD(PARAM)
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
%   savefilestr = sprintf('celllevel/result/COMBINE_F%dC%dN%dA%dM%dR%d_S%dW%d_C%dG%dN%d%s.mat', ...
%             flag, classnum, nfold, allflag, maximgperfold, randseed, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);
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

if param.allflag == 1
    combine_fold_all(param);
    return;
end

tic
addpath(genpath([pwd filesep 'SLIC']));

flag = param.flag;
classnum = param.classnum;
nfold = param.fold;


allflag = param.allflag;
if param.allflag == 3
    % set all flag to be 3 when calculating the probabilities:
    allflag = 1;    % for loading, will be changed back
end

maximgperfold = param.maximgperfold;
randseed = param.randseed;
sdaflag = param.sdaflag;
weightflag = param.weightflag;

cidx = param.cidx;
gidx = param.gidx;
nbrs = param.nbrs;

if param.classifierid == 1
    classifiername = 'libsvm';
end

if param.classifierid == 2
    classifiername = 'nn';
end

if param.classifierid == 3
   classifiername = 'libsvm';  % but use DAG for prediction
end

tmp =load([pwd filesep 'data' filesep 'yeastfeat_d' num2str(flag) '.mat']);
yeastdata = tmp.yeastdata;

cellconf = zeros(classnum,classnum);
cellconfvoting = cellconf;
imgconfvoting = cellconf;

for idxfold = 1:nfold
    modelfilestr = sprintf('celllevel/result/RES_F%dC%dN%dI%dA%dM%dR%d_S%dW%d_C%dG%dN%d%s.mat', ...
       flag, classnum, nfold, idxfold, allflag, maximgperfold, randseed, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);

    % load giant feature matrix file
    fprintf(1,' Loading classifier [%s] at %d th fold...\n',modelfilestr, idxfold);
    tmp=load(modelfilestr);
    model = tmp.model;

    [imagenums,order]=sort(yeastdata.ucsfresult.classSize, 'descend');
    testimgperclass = imagenums(1:classnum);
    testorder = order(1:classnum);
    imgperfold = min(floor(imagenums/nfold), maximgperfold);

    fprintf(1,' Loading rand perm index file and sda idx at %th fold...\n',idxfold);

    nfoldrandpermfilestr = sprintf('celllevel/result/PERM_F%dC%dN%dI%dA%dM%dR%d.mat', ...
       	      flag, classnum, nfold, idxfold, allflag, maximgperfold, randseed);

    tmp = load(nfoldrandpermfilestr);
    dataidx = tmp.dataidx;
    idx_sda = tmp.idx_sda;
    
    trainfeat = tmp.trainfeat;
    trainlabel = tmp.trainlabel;
    cellregnum = tmp.cellregnum;
    clear tmp;

    fprintf(1,'Testing classifier...\n\n')
        
    for t=1:classnum     
        testimgnum = length(dataidx{t}{idxfold});
        for e=1:testimgnum
            testfeat = yeastdata.ucsfresult.cellfeat{testorder(t)}{dataidx{t}{idxfold}(e)};
            testlabel = repmat(t, [size(yeastdata.ucsfresult.cellfeat{testorder(t)}{dataidx{t}{idxfold}(e)},1) 1]);
                
            testfeat = testfeat(:,idx_sda);
 
            [trainfeat_n, testfeat_n] = sc_featurenorm(trainfeat(:,idx_sda), testfeat);
            if param.classifierid == 3
                [predict_label, accuracy, prob_estimates] = svmpredictDAG(testlabel, testfeat_n, model); % test the training data        
            else
                [predict_label, accuracy, prob_estimates] = svmpredict(testlabel, testfeat_n, model); % test the training data     
            end   
            thisrow = zeros(1,classnum);

            for t2=1:classnum                       
                thisrow(t2) = sum(predict_label==t2);
            end  

            cellconf(t,:) = cellconf(t,:) + thisrow;   
            [votecount,voteclass] = max(thisrow);
            cellconfvoting(t,voteclass) = cellconfvoting(t,voteclass) + sum(thisrow);
            imgconfvoting(t,voteclass) = imgconfvoting(t,voteclass) + 1; 

            yeastdata.ucsfresult.imganot{testorder(t)}{dataidx{t}{idxfold}(e)}.correctcell = ...
                    [yeastdata.ucsfresult.imganot{testorder(t)}{dataidx{t}{idxfold}(e)}.correctcell sum(predict_label==t)];
            yeastdata.ucsfresult.imganot{testorder(t)}{dataidx{t}{idxfold}(e)}.voteclass = ...
                    [yeastdata.ucsfresult.imganot{testorder(t)}{dataidx{t}{idxfold}(e)}.voteclass voteclass];
            yeastdata.ucsfresult.imganot{testorder(t)}{dataidx{t}{idxfold}(e)}.votecorrect = ...
                    [yeastdata.ucsfresult.imganot{testorder(t)}{dataidx{t}{idxfold}(e)}.votecorrect (voteclass == t)];             
            yeastdata.ucsfresult.cellanot{testorder(t)}{dataidx{t}{idxfold}(e)}.classify = predict_label';  
            yeastdata.ucsfresult.cellanot{testorder(t)}{dataidx{t}{idxfold}(e)}.prob = prob_estimates;
    
            imgconfvotingnorm = imgconfvoting./(sum(imgconfvoting,2)*ones(1,classnum));
            fprintf(1,'image accuracy (unnormalized) = %.2f%%\n',100*sum(diag(imgconfvoting))/sum(imgconfvoting(:)));
        end
    end
end       

elapsetime = toc
yeastdata.ucsfresult.cellfeat = [];
yeastdata.ucsf = [];
yeastdata.procfiles = [];
yeastdata.geneNames=[];
yeastdata.imageIndices = [];


if param.allflag == 3
    % set all flag to be 3 when calculating the probabilities:
    allflag = 3;    % for loading, will be changed back
end

if param.classifierid == 3
    savefilestr = sprintf('celllevel/result/COMBINE_F%dC%dN%dA%dM%dR%d_S%dW%d_C%dG%dN%d%sDAG.mat', ...
       flag, classnum, nfold, allflag, maximgperfold, randseed, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);
    save(savefilestr,'yeastdata','cellconf','cellconfvoting','imgconfvoting','imgconfvotingnorm','elapsetime');
else
    savefilestr = sprintf('celllevel/result/COMBINE_F%dC%dN%dA%dM%dR%d_S%dW%d_C%dG%dN%d%s.mat', ...
       flag, classnum, nfold, allflag, maximgperfold, randseed, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);
    save(savefilestr,'yeastdata','cellconf','cellconfvoting','imgconfvoting','imgconfvotingnorm','elapsetime');
end

cellconf
cellconfvoting
imgconfvoting

fprintf(1,'\nNumber of classes = %d\n',classnum);
for i=1:classnum
    fprintf(1,'Class %d(%d): %s\n',i,yeastdata.ucsfresult.classSize(testorder(i)),yeastdata.ucsfresult.classNames{testorder(i)});  
end

fprintf(1,'image accuracy(normalized) = %.2f%%  ',100*sum(diag(imgconfvotingnorm))/sum(imgconfvotingnorm(:)));
fprintf(1,'(unnormalized) = %.2f%%\n\n',100*sum(diag(imgconfvoting))/sum(imgconfvoting(:)));

fprintf(1,'File %s created using %.4f minutes',savefilestr, elapsetime/60);
