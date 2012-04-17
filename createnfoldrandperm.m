function createnfoldrandperm(param)
%CREATENFOLDRANDPERM creates the index for training images, training features, and sda idx for all folds
%   CREATENFOLDRANDPERM(PARAM)
%
% Input:
%   PARAM.FLAG specifiy which feature matrix to use (either with unique location, ambiguous or punctate_composite class)
%   PARAM.CLASSNUM specify the number of class to classify
%   PARAM.NFOLD specify the number of fold to use for cross-validation
%   PARAM.IDXFOLD specify the n-th fold to use
%   PARAM.RANDSEED specify the rand seed to split the training/testing data
%   PARAM.ALLFLAG specify if we want to use all the images in the fold for training (if 1, maximgperfold will not be useful)
%   PARAM.MAXIMGPERFOLD specifiy the maximum number of images to for trainin at each fold
%
% Output:
%   nfoldrandpermfilestr = sprintf('celllevel/result/PERM_F%dC%dN%dI%dA%dM%dR%d.mat', ...^M
%	     flag, classnum, nfold, idxfold, allflag, maximgperfold, randseed);^M
%
%   See also createnfoldscriptR, sc_makegfm
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
tic
flag = param.flag;
classnum = param.classnum;
nfold = param.fold;
idxfold = param.idxfold;
allflag = param.allflag;
maximgperfold = param.maximgperfold;
randseed = param.randseed;

rand('seed',randseed);

if ~exist('celllevel/result','dir')
    mkdir('celllevel/result');
end

nfoldrandpermfilestr = sprintf('celllevel/result/PERM_F%dC%dN%dI%dA%dM%dR%d.mat', ...
    flag, classnum, nfold, idxfold, allflag, maximgperfold, randseed);

if exist(nfoldrandpermfilestr, 'file');
    fprintf(1, '%s exist...\n\n', nfoldrandpermfilestr);
    return;
end


fprintf(1,'Loading feature matrix...\n');
tmp =load([pwd filesep 'data' filesep 'yeastfeat_d' num2str(flag) '.mat']);
yeastdata = tmp.yeastdata;

[imagenums,order]=sort(yeastdata.ucsfresult.classSize, 'descend');


testimgperclass = imagenums(1:classnum);
testorder = order(1:classnum);
imgperfold = min(floor(imagenums/nfold), maximgperfold);

[imagenums,order]=sort(yeastdata.ucsfresult.classSize, 'descend');

for j=1:classnum
    idx{j} = randperm(imagenums(j));
    binsize = floor(imagenums(j)/nfold);

    start_idx = 1;
    for w=1:(nfold-mod(imagenums(j), nfold))
        dataidx{j}{w} = idx{j}(start_idx:(start_idx+binsize-1));
        start_idx  = start_idx+binsize;
    end

    for w=(nfold+1-mod(imagenums(j), nfold)):nfold
        dataidx{j}{w} = idx{j}(start_idx:(start_idx+binsize));
        start_idx  = start_idx+binsize+1;
    end
end

allfeatcell = cell(1,classnum);
for t=1:classnum
    fprintf(1,'Loading %d th class...', t);
    imnum = 0;
    if allflag == 0
        for e=1:imgperfold(t)
            for w=1:nfold
                if w ~= idxfold
                    imnum = imnum + 1;
                    allfeatcell{t} = [allfeatcell{t}; yeastdata.ucsfresult.cellfeat{testorder(t)}{dataidx{t}{w}(e)}];
                end
            end
        end
    else
        for w=1:nfold
            if w ~= idxfold
                for e=1:length(dataidx{t}{w})
                    imnum = imnum + 1;
                    allfeatcell{t} = [allfeatcell{t}; yeastdata.ucsfresult.cellfeat{testorder(t)}{dataidx{t}{w}(e)}];
                end
            end
        end
    end
    fprintf(1,'(%d cell regions in %d images were loaded)\n', size(allfeatcell{t},1), imnum);
end

trainfeat = [];   trainlabel = [];  cellregnum = [];
for t=1:classnum
    trainfeat = [trainfeat; allfeatcell{t}];
    trainlabel = [trainlabel; repmat(t,[size(allfeatcell{t},1),1])];
    cellregnum = [cellregnum size(allfeatcell{t},1)];
end

fprintf(1,' Doing SDA...\n');
while( license('checkout','statistics_toolbox') == 0 )
    fprintf(1,'Waiting for license of statistics_toolbox...\n');
    pause(120);
end
[idx_sda, ignoreidx] = ml_stepdisc(allfeatcell,[nfoldrandpermfilestr '.sdalog']);

elapsetime = toc;
save(nfoldrandpermfilestr, 'dataidx', 'idx_sda','trainfeat','trainlabel','cellregnum', 'elapsetime','param');

fprintf(1,'File %s created using %.4f minutes',nfoldrandpermfilestr, elapsetime);
