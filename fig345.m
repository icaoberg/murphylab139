% This script creates Figure 3, Figure 4, and Figure 5

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

clear all;

fprintf(1,'Generating Figure 3, Figure 4, and Figure 5...\n');

load('data/yeastdata.mat');
addpath(genpath([pwd filesep 'SLIC']));

orfs = {'YAL009W','YGR130C','YFL034W'};
regions = {[390 510;290 410];[75 195;70 190];[30 150;100 220]};
dnathrs = [20000 20000 30000];
protthrs = [20000 Inf Inf];


savedir = [pwd filesep 'data'];
segdir = [pwd filesep 'celllevel/mask/loopy2_5/'];

for i=1:length(orfs)
    orf = orfs{i};
    idx = strmatch(orf,yeastdata.geneNames);
    dnaimg = imread(yeastdata.procfiles{1}{idx});
    dnaimg = dnaimg(regions{i}(1,1):regions{i}(1,2), ...
                    regions{i}(2,1):regions{i}(2,2));
    dnaimg = imresize(dnaimg,2,'bilinear');            
    pixelidx = find(dnaimg>dnathrs(i));
    if ~isempty(pixelidx)
        dnaimg2 = dnaimg;
        dnaimg2(pixelidx) = 0;
        dnaimg(pixelidx) = max(dnaimg2(:));
    end      

    protimg = imread(yeastdata.procfiles{3}{idx});
    protimg = protimg(regions{i}(1,1):regions{i}(1,2), ...
                      regions{i}(2,1):regions{i}(2,2));
    protimg = imresize(protimg,2,'bilinear');
    pixelidx = find(protimg>protthrs(i));
    if ~isempty(pixelidx)
        protimg2 = protimg;
        protimg2(pixelidx) = 0;
        protimg(pixelidx) = max(protimg2(:));
    end
    filename = tz_yeastimg2featname(yeastdata.procfiles{2}{idx});
    segfilepath = [segdir filesep filename '.seg.mat'];
    tmp = load(segfilepath);
    cellimg = edge(imfill(tmp.segimg,'hole'),'sobel',0);
    cellimg = cellimg(regions{i}(1,1):regions{i}(1,2), ...
                      regions{i}(2,1):regions{i}(2,2));
    cellimg = imresize(cellimg,2,'bilinear');
    rgbimg = ml_synrgbimg(dnaimg,protimg,cellimg,'isc');
    rgbimg = mat2gray(rgbimg);
    
    rgbimg = imresize(rgbimg,[1039/242]);
    imwrite(rgbimg,[savedir filesep 'Fig' num2str(i+2) '.tif'],'tif','Resolution',300);
end

fprintf(1,'Generation of data/Fig3.tif data/Fig4.tif data/Fig5.tif is done!\n\n');
