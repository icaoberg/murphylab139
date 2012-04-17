function classify_other(param)
%CLASSIFY_OTHER uses classifier (with flag) to train image with otherflag, only used for ambiguous or punctate_composite class
%   CLASSIFY_OTHER(PARAM)
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
%   See also 
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

% New field... what 'other' data source, either 41 (ambiguous) or 42 (punctate_composite)
otherflag = param.otherflag


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
end

if param.classifierid == 2
    classifiername = 'nn';
end

if param.classifierid == 3
   classifiername = 'libsvm';  % but use DAG for prediction
end


tmp =load([pwd filesep 'data' filesep 'yeastfeat_d' num2str(otherflag) '.mat']);
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
        
    testimgnum = yeastdata.ucsfresult.classSize;
    for e=1:testimgnum
        testfeat = yeastdata.ucsfresult.cellfeat{1}{e};
        testlabel = repmat(1, [size(yeastdata.ucsfresult.cellfeat{1}{e},1) 1]);
                
        testfeat = testfeat(:,idx_sda);
 
        [trainfeat_n, testfeat_n] = sc_featurenorm(trainfeat(:,idx_sda), testfeat);
        if param.randseed <= 8
            [predict_label, accuracy, prob_estimates] = svmpredict(testlabel, testfeat_n, model, '-b 1'); % test the training data        
        else
            [predict_label, accuracy, prob_estimates] = svmpredict(testlabel, testfeat_n, model); % test the training data     
        end   
        thisrow = zeros(1,classnum);

        for t2=1:classnum                       
            thisrow(t2) = sum(predict_label==t2);
        end  

        cellconf(1,:) = cellconf(1,:) + thisrow;   
        [votecount,voteclass] = max(thisrow);
        cellconfvoting(1,voteclass) = cellconfvoting(1,voteclass) + sum(thisrow);
        imgconfvoting(1,voteclass) = imgconfvoting(1,voteclass) + 1; 

        yeastdata.ucsfresult.imganot{1}{e}.correctcell = ...
	            [yeastdata.ucsfresult.imganot{1}{e}.correctcell sum(predict_label==1)];
        yeastdata.ucsfresult.imganot{1}{e}.voteclass = ...
                    [yeastdata.ucsfresult.imganot{1}{e}.voteclass voteclass];
        yeastdata.ucsfresult.imganot{1}{e}.votecorrect = ...
                    [yeastdata.ucsfresult.imganot{1}{e}.votecorrect (voteclass == 1)];
        if ~iscell(yeastdata.ucsfresult.cellanot{1}{e}.classify(1))
            yeastdata.ucsfresult.cellanot{1}{e}.classify = [];
        end	
        yeastdata.ucsfresult.cellanot{1}{e}.classify{idxfold} = predict_label';  
%        yeastdata.ucsfresult.cellanot{1}{e}.prob{idxfold} = prob_estimates;
        imgconfvotingnorm = imgconfvoting./(sum(imgconfvoting,2)*ones(1,classnum));
        fprintf(1,'image accuracy (unnormalized) = %.2f%%\n',100*sum(diag(imgconfvoting))/sum(imgconfvoting(:)));
    end
end       

elapsetime = toc

yeastdata.ucsfresult.cellfeat = [];
yeastdata.ucsf = [];
yeastdata.procfiles = [];
yeastdata.geneNames=[];
yeastdata.imageIndices = [];


if param.classifierid == 3
    savefilestr = sprintf('celllevel/result/COMBINE_F%dC%dN%dA%dM%dR%d_S%dW%d_C%dG%dN%d%sDAG.mat', ...
       otherflag, classnum, nfold, allflag, maximgperfold, randseed, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);
    save(savefilestr,'yeastdata','cellconf','cellconfvoting','imgconfvoting','imgconfvotingnorm','elapsetime');
else
    savefilestr = sprintf('celllevel/result/COMBINE_F%dC%dN%dA%dM%dR%d_S%dW%d_C%dG%dN%d%s.mat', ...
       otherflag, classnum, nfold, allflag, maximgperfold, randseed, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);
    save(savefilestr,'yeastdata','cellconf','cellconfvoting','imgconfvoting','imgconfvotingnorm','elapsetime');
end

cellconf
cellconfvoting
imgconfvoting


fprintf(1,'File %s created using %.4f minutes',savefilestr, elapsetime/60);
