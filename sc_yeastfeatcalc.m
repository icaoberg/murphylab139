function sc_yeastfeatcalc(method)
%SC_YEASTFEATCALC  A function for cell/field level feature calculation on yeast images
%   when method.celllevel == 1, do cell level feature calculation 
%   when method.celllevel == 0, do field level feature calculation 
%   22-Nov-2005 Initial write T. Zhao
%   11-Apr-2007 Modified S.C. Chen

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

imglist = load(method.procfilesmatname);
nprocfiles = length(imglist.procfiles{1});

if method.celllevel == 1
    resultdir180feat = [method.resultdir filesep 'feat_b256_d1'];
    resultdir = [method.resultdir filesep 'feat_b' num2str(method.har_intbins) '_d' num2str(method.downsamplerate)];
    %mkdir(resultdir);    % Do not store the features from masks without postprocessing
    resultposdir = [resultdir '_pos'];   
    segdirpath = [pwd filesep 'celllevel' filesep 'mask' filesep ...
		 'loopy' num2str(method.alpha) '_' num2str(method.beta)];

    % load the parameters of logistic regression for post-processing
    load([pwd filesep 'data' filesep 'logisparam.mat']);  
else
    resultposdir = [pwd filesep 'fieldlevel' filesep 'feat' filesep 'feat_b' num2str(method.har_intbins) '_d' num2str(method.downsamplerate)];
end
mkdir(resultposdir);

for n=1:nprocfiles
    segImagePath = imglist.procfiles{1}{n};
    [segdir,segfile,ext] = fileparts(segImagePath);
    [chr,fileprefix,channel] = tz_parseyeastname(segfile);
    featFileName = [segdir(findstr('C',segdir):end) '_' chr '_' fileprefix];
    
    controlDirectory = [resultposdir filesep featFileName];
    [s,msg] = mkdir(resultposdir,featFileName);

    % Check if the job has been taken over
    if strfind(msg,'exist')
        continue;
    end
    
    fprintf(1, 'Processing %dth image(%s)...\n',n,segImagePath);
    resultposfile = [resultposdir filesep featFileName '.feat.mat'];

    if ~exist(resultposfile,'file') & method.celllevel == 1
        if ~exist([segdirpath '/' featFileName '_DIC.seg.mat'],'file')            
            fprintf(1,'The segmentation file does not exist\n');	    
            continue;
        else
            load([segdirpath '/' featFileName '_DIC.seg.mat']);
            fprintf(1,'Reading segmentation file...');      
            normreg=tz_normreg(segimg);
            if sum(segimg(:)) == 0
                fprintf(1, 'No feature can be calculated from empty masks\n');
                continue;
            end
        end

        gfpimgname = imglist.procfiles{3}{n};
        if ~exist(gfpimgname,'file')
            fprintf(1,'GFP image %s does not exist. Skip it...\n',gfpimgname);
            continue;
	else
            fprintf(1,'Reading GFP image...');
	    gfpimg=imread(gfpimgname);
        end

        dnaimgname = imglist.procfiles{1}{n};
	if ~exist(dnaimgname,'file')
            fprintf(1,'DNA image %s does not exist. Skip it...\n',dnaimgname);
            continue;  
	else
            fprintf(1,'Reading DNA image...');
            dnaimg=imread(dnaimgname);
        end	

        feats = [];
        if method.har_intbins == 256 
            feats=sc_regionfeatmat(gfpimg,dnaimg,normreg,'all180',method.har_intbins,method.downsamplerate);
        else	  
   	    resultfile180feat = [resultdir180feat filesep featFileName '.feat.mat'];
            if ~exist(resultfile180feat,'file')
                feats=sc_regionfeatmat(gfpimg,dnaimg,normreg,'all180',method.har_intbins,method.downsamplerate);
	    else
    	        harfeats=sc_regionfeatmat(gfpimg,dnaimg,normreg,'har',method.har_intbins,method.downsamplerate);
                load(resultfile180feat);
                feats(:,73:85) = harfeats;
            end
        end

        fprintf(1,'Before post-processing: %d masks',size(feats,1));

        [tmpfeat feat] = ml_featurenorm(logisparam.trainfeat, feat);
    	augcellfeat = [ones(size(feat,1),1) feat];

        augcellscore = augcellfeat * logisparam.wts;
        posflag = find(augcellscore >= logisparam.stopscore);
        negflag = find(augcellscore < logisparam.stopscore);
        feats(negflag,:) = [];      
        
        augcellscore(negflag) = [];
        [tmp, idx] = sort(augcellscore);
        idx = flipud(idx);
        feats = feats(idx,:);
                      
        fprintf(1,'\tAfter post-processing: %d masks\n',size(feats,1));       
        save(resultposfile,'feats','posflag','negflag');        
    end

    if ~exist(resultposfile,'file') & method.celllevel == 0
        gfpimgname = imglist.procfiles{3}{n};
        if ~exist(gfpimgname,'file')
            fprintf(1,'GFP image %s does not exist. Skip it...\n',gfpimgname);
            continue;
        else
            fprintf(1,'Reading GFP image...');
            gfpimg=imread(gfpimgname);
        end

        dnaimgname = imglist.procfiles{1}{n};
        if ~exist(dnaimgname,'file')
            fprintf(1,'DNA image %s does not exist. Skip it...\n',dnaimgname);
            continue;
        else
            fprintf(1,'Reading DNA image...');
            dnaimg=imread(dnaimgname);
        end
	feats=sc_regionfeatmat(gfpimg,dnaimg, -1, 'mcell_all', method.har_intbins,method.downsamplerate);
        save(resultposfile,'feats');        
    end
    fprintf(1,'%s saved\n\n', resultposfile);
    unix(['rm -rf ' controlDirectory]);
end
