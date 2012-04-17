% This script generates feature matrix for field-level classification

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

clear;

har_intbins = 128;
downsamplerate = 1;
whichdatabase = 2;

resultdir = [pwd filesep 'data'];
featrootroot = [pwd filesep 'fieldlevel' filesep 'feat'];
featroot = ['feat_b' num2str(har_intbins) '_d' num2str(downsamplerate)];
featdir = [featrootroot filesep featroot];

load([pwd filesep 'data' filesep 'yeastdata.mat']);

%yeastdata = yeastdataunigene;

if whichdatabase == 2
    featfilename = [featroot '_ucsf.mat'];
    allclassLabels = yeastdata.ucsf.classLabels;
    locnums = yeastdata.ucsf.locNums;
    strClassNames = yeastdata.ucsf.strClassNames;
else
    fprintf(1,'Please use UCSF database as the ground truth\n');
    return;
end












savefilename = [pwd filesep 'data' filesep 'yeastfieldfeat.mat'];


classLabels = allclassLabels(find(locnums==1),:);
uniquelabels = unique(classLabels); % this is our new label 1:17 or 1:22
uniquelabels(end) = [];


classNum = size(uniquelabels, 1);
classSize = zeros(1, classNum);
classNames = cell(1, classNum);
classidx = cell(1, classNum);

for i=1:classNum
    classidx{i} = find(allclassLabels == uniquelabels(i));
    classSize(i) = length(classidx{i});
    classNames{i} = strClassNames{uniquelabels(i)};
end

featbyclass = cell(classNum, max(classSize));
countbyclass = ones(1, classNum);

fprintf(1,'Executing %s\n', featdir);
procfiles = yeastdata.procfiles;

prefix = 'feats';
feat = [];

for k=1:length(procfiles{1})
    thisclass = find(uniquelabels == allclassLabels(k));
    if length(thisclass) == 0 
        continue;
    end        
    
    if ~exist('procfiles','var')
        featfile = [prefix num2str(k)];
    else
        featfile = tz_yeastimg2featname(procfiles{2}{k});
    end
    fprintf(1,'%d.',k);
    featfile = featfile(1: size(featfile,2)-4);
    featfile = [featfile '.feat.mat'];    
    featpath = [featdir filesep featfile];

    if exist(featpath)
        load(featpath);    
    else
        continue;
    end 
    if ~isempty(feats)        
        featbyclass{thisclass}{countbyclass(thisclass)} = sc_removenan(feats);
        countbyclass(thisclass) = countbyclass(thisclass) + 1;
    end
end
classSize = countbyclass - 1;
save(savefilename);
fprintf(1,'\nLoading is done and %s is created\n\n', savefilename);
