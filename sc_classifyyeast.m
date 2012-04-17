function sc_classifyyeast(param)
%function sc_classifyyeast
%SC_CLASSIFYYEAST runs cell level classification on unique location classes with purality voting scheme
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
%   savefilestr = sprintf('celllevel/result/RES_F%dC%dN%dI%dA%dM%dR%d_S%dW%d_C%dG%dN%d%s.mat', ...^M
%             flag, classnum, nfold, idxfold, allflag, maximgperfold, randseed, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);^M
%
%   See also createnfoldscriptR, createnfoldrandperm
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

%   11-Apr-2007 Initial write S.C. Chen
%   Copyright (c) Center for Bioimage Informatics, CMU

tic
addpath(genpath([pwd filesep 'SLIC']));

flag = param.flag;
classnum = param.classnum;
nfold = param.fold;
idxfold = param.idxfold;
allflag = param.allflag;
maximgperfold = param.maximgperfold;
randseed = param.randseed;
sdaflag = param.sdaflag;
weightflag = param.weightflag;

cidx = param.cidx;
gidx = param.gidx;
nbrs = param.nbrs;

%if param.classifier == 'nn'
%param.cidx=0; param.gidx=0;

%if param.classifier == 'libsvm'
%param.nbrs = 0;
%param.suffix

if param.classifierid == 1
    classifiername = 'libsvm';
else
    classifiername = 'nn';
end

savefilestr = sprintf('celllevel/result/RES_F%dC%dN%dI%dA%dM%dR%d_S%dW%d_C%dG%dN%d%s.mat', ...
     flag, classnum, nfold, idxfold, allflag, maximgperfold, randseed, sdaflag, weightflag, cidx, gidx, nbrs, classifiername);

if exist(savefilestr,'file')
     fprintf(1,'%s exists, continue......\n',savefilestr);
     return;    
end

savefilestr 

% load giant feature matrix file
fprintf(1,' Loading feature matrix...\n');
tmp =load([pwd filesep 'data' filesep 'yeastfeat_d' num2str(flag) '.mat']);
yeastdata = tmp.yeastdata;
clear tmp;

[imagenums,order]=sort(yeastdata.ucsfresult.classSize, 'descend');
testimgperclass = imagenums(1:classnum);
testorder = order(1:classnum);
imgperfold = min(floor(imagenums/nfold), maximgperfold);

fprintf(1,' Loading rand perm index file and sda idx for nfold cross-validation...\n');

nfoldrandpermfilestr = sprintf('celllevel/result/PERM_F%dC%dN%dI%dA%dM%dR%d.mat', ...
       	      flag, classnum, nfold, idxfold, allflag, maximgperfold, randseed);
if exist(nfoldrandpermfilestr,'file')
     tmp = load(nfoldrandpermfilestr);
     dataidx = tmp.dataidx;
     idx_sda = tmp.idx_sda;
     if param.sdaflag == 0
         idx_sda = [1:245];
     end
     trainfeat = tmp.trainfeat;
     trainlabel = tmp.trainlabel;
     cellregnum = tmp.cellregnum;

     clear tmp;
else
     fprintf(1,'====The randperm file and SDA index %s does not exist====\n', nfoldrandpermfilestr);
     return;
end


fprintf(1,'Training classifier...\n');

trainfeat = trainfeat(:,idx_sda);    tmp = trainfeat(1,:);
[trainfeat_n, tmp] = sc_featurenorm(trainfeat, tmp);

% down weighted
scalec = cellregnum / cellregnum(end);
scalec = 1./ scalec;

param.c=2.^[-5:2:15];  %11
param.g=2.^[-15:2:3];  %10

option = ['-s 0 -c ' num2str(param.c(cidx)) ' -g ' num2str(param.g(gidx)) ' -b 1 '];
if weightflag == 1
    for i=1:length(scalec)
         op = sprintf(' -w%d %f',i,scalec(i));
         option = [option op];
    end
end
option

model = svmtrain(trainlabel, trainfeat_n, option);

cellconf = zeros(classnum,classnum);
cellconfvoting = cellconf;
imgconfvoting = cellconf;

fprintf(1,'Testing classifier...\n\n')
        
for t=1:classnum     
    if allflag == 1
        testimgnum = length(dataidx{t}{idxfold});
    else
        testimgnum = imgperfold(t);
    end
    for e=1:testimgnum

        testfeat = yeastdata.ucsfresult.cellfeat{testorder(t)}{dataidx{t}{idxfold}(e)};
        testlabel = repmat(t, [size(yeastdata.ucsfresult.cellfeat{testorder(t)}{dataidx{t}{idxfold}(e)},1) 1]);
                
        testfeat = testfeat(:,idx_sda);
  
        [trainfeat_n, testfeat_n] = sc_featurenorm(trainfeat, testfeat);
        [predict_label, accuracy, dec_values] = svmpredict(testlabel, testfeat_n, model); % test the training data        
                
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
    
        imgconfvotingnorm = imgconfvoting./(sum(imgconfvoting,2)*ones(1,classnum));

        fprintf(1,'image accuracy = %.2f%%\n',100*sum(diag(imgconfvotingnorm))/sum(imgconfvotingnorm(:)));
        fprintf(1,'image accuracy (adjusted) = %.2f%%\n',100*sum(diag(imgconfvoting))/sum(imgconfvoting(:)));
    end
end       

%clear allcellidx allfeat allfeatcell allfeatsdatrain alllabel alllabelsdatrain ans
elapsetime = toc
%clear trainfeat testfeat trainlabel tmp

yeastdata.ucsfresult.cellfeat = [];
save(savefilestr,'yeastdata','model','elapsetime');

cellconf
cellconfvoting
imgconfvoting

fprintf(1,'\nNumber of classes = %d\n',classnum);
for i=1:classnum
    fprintf(1,'Class %d(%d): %s\n',i,yeastdata.ucsfresult.classSize(testorder(i)),yeastdata.ucsfresult.classNames{testorder(i)});  
end

fprintf(1,'image accuracy(normalized) = %.2f%%  ',100*sum(diag(imgconfvotingnorm))/sum(imgconfvotingnorm(:)));
fprintf(1,'(unnormalized) = %.2f%%\n\n',100*sum(diag(imgconfvoting))/sum(imgconfvoting(:)));

whos trainlabel

fprintf(1,'File %s created using %.4f minutes',savefilestr, elapsetime);
