function sc_segyeast(datafile, method, resultdir)
%SC_SEGYEAST segment yeast cell images into regions
%   SC_SEGYEAST(DATAFILE,METHOD,RESULTDIR)
%   
% Input:
%   DATAFILE is the filename of the list of images to be processed
%   METHOD is a struct which has the information the type and parameters of
%   three segmentation methods.
%   To run spectral-clustering, you have to set method.id=1
%   To run loopy belief propagation, you have to set method.id=2, and
%   specify the two strength parameters of the voting potential, method.alpha
%   and method.beta.  
%   To run seeded watershed, you have to set method.id=3 and specify the
%   scale to the RC threshold.  
%   RESULTDIR specifies the path to store the segmented masks
%
% Output:
%   This function stores the segmented masks into the RESULTDIR
%
% usage: 
%   For spectral-clustering:
%      method.id=1;
%      sc_segyeast(procfilesmatname, method);
%   For loopy belief propagation with voting potential:
%      method.id=2; method.alpha=2; method.beta=5;
%      sc_segyeast(procfilesmatname, method);
%   For seeded watershed segmentation:
%      method.id=3; method.scale=1.5;
%      sc_segyeast(procfilesmatname, method);
%
%   See also: segscript, segscriptallyeast
%
%   22-Nov-2005 Initial written by T. Zhao
%   09-Aug-2006 Modified by S. C. Chen

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

if nargin < 1
    error('At least one argument is required')
end

if license('checkout','statistics_toolbox')==0
  % the segmentation algorithm needs to get the license handle of the statistics toolbox
  error('The statistics toolbox is not available.  Please try again later.')
  return;
end

imglist = load(datafile);
nprocfiles = length(imglist.procfiles{1});

fprintf(1, 'Total %d image will be segmented...\n',nprocfiles);

if method.id == 1
    % spectral-clustering
    subdir = ['spec' num2str(method.seed)];
elseif method.id == 2
    % loopy
    subdir = ['loopy' num2str(method.alpha) '_' num2str(method.beta)];
elseif method.id == 3
    % seeded watershed
    subdir = ['ws' num2str(10*method.scale)];    
end

if nargin < 3
    resultdir = imglist.rootdir;
end

% These images in the yeast database cannot be segmented
skiplist = [86, 538, 880, 1707, 2909, 3000, 2643, 3407, 2056];

for n=1:nprocfiles

    dnaImagePath = imglist.procfiles{1}{n};
    dicImagePath = imglist.procfiles{2}{n};
    if find(n == skiplist)
        fprintf(1, '%dth image(%s) cannot be segmented, so skip it...\n',n,dicImagePath);   
        continue;
    end
    
    % specify the segmentation file name (for saving)
    segFileName = [tz_yeastimg2featname(dicImagePath) '.seg'];
    % specify the path for the controlDirectory
    controlDirectory = [resultdir filesep subdir filesep segFileName];
    
    % make a directory similar to the segFileName; if this directory
    % exists, the segmentation job has been taken over by another node    
    % skip this image and work on the next image
    [s,msg] = mkdir([resultdir filesep subdir], segFileName); 
    if strfind(msg,'exist')
        continue;
    end
    
    fprintf(1, 'Processing %dth image(%s)...\n',n,dicImagePath);

    % specify the absolute path to the result file name (for saving)
    resultfile = [resultdir filesep subdir filesep segFileName '.mat'];
    
    if ~exist(resultfile,'file')
        %calculate the processing time
        t0 = clock;
        % use controlDirectory to save the intermediate results

        method.imgidx = n;
        [segimg, feat] = sc_cutcells(dnaImagePath,dicImagePath,method,controlDirectory);
        
        t = etime(clock, t0);
        
        if findstr(version,'6.')
            save(resultfile,'segimg','feat','t');
        else
            save(resultfile,'segimg','feat','t','-V6');
        end              
    end
    unix(['rm -r ' controlDirectory]);   
end
