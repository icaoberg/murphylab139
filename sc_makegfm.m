function sc_makegfm(classflag)
%SC_MAKEGFM loads the calculated cell level features and store them in yeastdata.ucsfresult 
%   When classflag == 1, prepare features for 21 unique location patterns
%   When classflag == 2, prepare features for mixture location patterns 
%   When classflag == 3, prepare features for 20 unique location patterns (no punctate_composite)
%   When classflag == 41, prepare features for ambiguous class
%   When classflag == 42, prepare features for punctate_composite class
%
%   See also makedata4

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

%   06-June-2007 Initial written by S.C. Chen
%   Copyright (c) Center for Bioimage Informatics, CMU

har_intbins = 128;
downsamplerate = 1;

resultdir = [pwd filesep 'data'];
featrootroot = [pwd filesep 'celllevel' filesep 'feat'];
featroot = ['feat_b' num2str(har_intbins) '_d' num2str(downsamplerate) '_pos'];
featdir = [featrootroot filesep featroot];

load([pwd filesep 'data' filesep 'yeastdata.mat']);
% We use unique gene name for cell level classification

% classflag == 1, preparing features stored in yeastdata.ucsfresult for cell level classification on unique location patterns 
if classflag == 1
    allclassLabels = yeastdata.ucsf.classLabels;
    locnums = yeastdata.ucsf.locNums;
    classLabels = allclassLabels(find(locnums==1),:); % unique labels
    uniquelabels = unique(classLabels); % this is our new label 1:17 or 1:22
    % Do not use ambiguous
    uniquelabels(end) = [];
end

% classflag == 2, preparing features stored in yeastdata.ucsfresult for cell level classification on mixture class
if classflag == 2
    allclassLabels = yeastdata.ucsf.classLabels;
    locnums = yeastdata.ucsf.locNums;
    classLabels = allclassLabels(find(locnums~=1),:); % non-unique labels excluding unvisualized data
    uniquelabels = unique(classLabels);
    uniquelabels = [uniquelabels];
    % exclude unvisualized category
    uniquelabels(1) = [];
end


% classflag == 3, preparing features stored in yeastdata.ucsfresult for cell level classification on unique location patterns
if classflag == 3  % 20 classes, no 'punctate_composite'  
    % use unique genes
    yeastdata = yeastdataunigene;

    allclassLabels = yeastdata.ucsf.classLabels;
    locnums = yeastdata.ucsf.locNums;

    classLabels = allclassLabels(find(locnums==1),:); % unique labels
    uniquelabels = unique(classLabels); % this is our new label 1:17 or 1:22

e    setempty = [];
    removesetname = {'ambiguous','punctate_composite'};
    for i=1:length(removesetname)
        setempty = [setempty; strmatch(removesetname{i},yeastdata.ucsf.strClassNames(uniquelabels))];
    end
    uniquelabels(setempty) = [];
end

if classflag == 41
    % use unique genes
    yeastdata = yeastdataunigene;

    allclassLabels = yeastdata.ucsf.classLabels;
    locnums = yeastdata.ucsf.locNums;

    % The ambiguous class
    uniquelabels = 97;
end

if classflag == 42
    % use unique genes
    yeastdata = yeastdataunigene;  

    allclassLabels = yeastdata.ucsf.classLabels;
    locnums = yeastdata.ucsf.locNums;

    % The punctate_composite class
    uniquelabels = 46;
end

strClassNames = yeastdata.ucsf.strClassNames;


savefilename = [pwd filesep 'data' filesep 'yeastfeat_d' num2str(classflag) '.mat'];


classNum = size(uniquelabels, 1);
classSize = zeros(1, classNum);
classNames = cell(1, classNum);
classidx = cell(1, classNum);

for i=1:classNum
    classidx{i} = find(allclassLabels == uniquelabels(i));
    classSize(i) = length(classidx{i});
    classNames{i} = strClassNames{uniquelabels(i)};
end

% add by Sam on 11/03/06, for complete analysis
yeastdata.ucsfresult.classNames = classNames;
yeastdata.ucsfresult.classSize = classSize;
yeastdata.ucsfresult.classidx = classidx;
yeastdata.ucsfresult.classNum = classNum;

yeastdata.ucsfresult.imganot = cell(1, classNum);
yeastdata.ucsfresult.cellanot = cell(1, classNum);
yeastdata.ucsfresult.cellfeat = cell(1, classNum);

%  yeastdata.ucsfresult.imganot{i}{j}, class, img
failimg = [];    % class, imgidx

for i=1:classNum
    yeastdata.ucsfresult.imganot{i} = cell(1, yeastdata.ucsfresult.classSize(i));
    yeastdata.ucsfresult.cellanot{i} = cell(1, yeastdata.ucsfresult.classSize(i));
    yeastdata.ucsfresult.cellfeat{i} = cell(1, yeastdata.ucsfresult.classSize(i));

    for j=1:yeastdata.ucsfresult.classSize(i)
        % read in the segmentation & the original image

        idx = yeastdata.ucsfresult.classidx{i}(j);
        featfile = tz_yeastimg2featname(yeastdata.procfiles{2}{idx});
        fprintf(1,'(c%d-i%d)',i,j);

        featfile = featfile(1: size(featfile,2)-4);
        featfile = [featfile '.feat.mat'];
        featpath = [featdir filesep featfile];

        if exist(featpath)
            imgfeat = load(featpath);
	    if isempty(sc_removenan(imgfeat.feats))
                % this feature file is empty
                yeastdata.ucsfresult.imganot{i}{j} = NaN;
                yeastdata.ucsfresult.cellanot{i}{j} = NaN;
                yeastdata.ucsfresult.cellfeat{i}{j} = NaN;
                failimg = [failimg; i j];
                continue;
            end
	else
            % this feature file does not exist
            yeastdata.ucsfresult.imganot{i}{j} = NaN;
            yeastdata.ucsfresult.cellanot{i}{j} = NaN; 
            yeastdata.ucsfresult.cellfeat{i}{j} = NaN;
            failimg = [failimg; i j];
            continue;
        end

	feats = imgfeat.feats;
	newfeatpath = featpath;   newfeatpath(strfind(newfeatpath, 'd1')+1)='2';    newfeat = load(newfeatpath);
        feats = [feats newfeat.feats(:,73:85)];
        newfeatpath = featpath;   newfeatpath(strfind(newfeatpath, 'd1')+1)='3';    newfeat = load(newfeatpath);
        feats = [feats newfeat.feats(:,73:85)];
        newfeatpath = featpath;   newfeatpath(strfind(newfeatpath, 'd1')+1)='4';    newfeat = load(newfeatpath);
        feats = [feats newfeat.feats(:,73:85)];
        newfeatpath = featpath;   newfeatpath(strfind(newfeatpath, 'd1')+1)='5';    newfeat = load(newfeatpath);
        feats = [feats newfeat.feats(:,73:85)];
        newfeatpath = featpath;   newfeatpath(strfind(newfeatpath, 'd1')+1)='6';    newfeat = load(newfeatpath);
        feats = [feats newfeat.feats(:,73:85)];

        [feats validrow, anarow] = sc_removenan(feats);      
        yeastdata.ucsfresult.cellfeat{i}{j} = feats;

        imginfo.numcell = size(feats,1);
        if imginfo.numcell == 0
            failimg = [failimg; i j];
            continue;
        end
        imginfo.correctcell = [];        
        imginfo.voteclass = [];
        imginfo.votecorrect = [];
        yeastdata.ucsfresult.imganot{i}{j} = imginfo;

        cellinfo.classify = NaN * zeros(1, imginfo.numcell);
        cellinfo.validrow = validrow;
        cellinfo.anarow = anarow;
      
	yeastdata.ucsfresult.cellanot{i}{j} = cellinfo;
    end
    fprintf(1,'\n\n');
end

% do post-processing using failimg to remove imganot{i}{j}, cellanot{i}{j}, cellfeat{i}{j}
for i=size(failimg,1):-1:1
    yeastdata.ucsfresult.classSize(failimg(i,1)) = yeastdata.ucsfresult.classSize(failimg(i,1))-1;
    yeastdata.ucsfresult.classidx{failimg(i,1)}(failimg(i,2)) = [];
    yeastdata.ucsfresult.imganot{failimg(i,1)}(failimg(i,2)) = [];
    yeastdata.ucsfresult.cellanot{failimg(i,1)}(failimg(i,2)) = [];
    yeastdata.ucsfresult.cellfeat{failimg(i,1)}(failimg(i,2)) = [];
end

save(savefilename ,'yeastdata','uniquelabels','allclassLabels','strClassNames'); 
fprintf(1,'File %s is saved\n', savefilename);
