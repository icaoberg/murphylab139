function [segimg, feat] = sc_cutcells(dnaImagePath,dicImagePath,method,controlDirectory)
%SC_CUTCELLS A function to identify a subimage of interest in the DIC image
%   and run loopy belief segmentation on the subimage. After segmentation,
%   masks caused by artifects are removed and features of each mask are calculated.
%   [SEGIMG, FEAT] = SC_CUTCELLS(DNAIMAGEPATH,DICIMAGEPATH,METHOD,CONTROLDIRECTORY)
%   
% Input:
%   DNAIMAGEPATH is the path for DNA image
%   DICIMAGEPATH is the path for DIC image
%   METHOD is a struct which has the information the type and parameters of
%   three segmentation methods.  Please see sc_segyeast for details
%   CONTROLDIRECTORY is the path for the control directory
%   The control directory is used to indicate if the image has been
%   processing and also store the intermediate segmentation results
%
% Output:
%   SEGIMG is the index image for the masks
%   FEAT is the 10 features 
%
% usage: 
%   [segimg, feat] = sc_cutcells(dnaImagePath,dicImagePath,method,controlDirectory)
%
%   See also: sc_segyeast, segscriptallyeast, sc_beliefprop
%
%   22-Nov-2005 Initial written by A. Shanbhag
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

if method.id == 1
    % spectral clustering
    fprintf(1,'Spectral clustering is not available in this release\n');
    feat = [];  segimg = [];    
    return;
end

if method.id == 3
    % run seeded watershed
    if isempty( which('mmcwatershed') )
        fprintf(1,'Seeded watershed is not available in this release\n');
        feat = [];  segimg = [];    
        return;
    else        
	dnaim = double(imread(dnaImagePath));
        seedim = ml_preprocess_scale(dnaim, [], 'ml', 'quantile','rc',method.scale) > 0;
        dnaim = uint16(int_invert(dnaim));
        segimg = mmcwatershed(dnaim, seedim, mmsebox, 'regions');

        feat = [];
        for i=1:max(segimg(:))
            maski = (segimg==i);
            if sum(sum(dicbgsub.*(maski)))/sum(maski(:)) < 0.01
                segimg(find(segimg==i)) = 0;
            else
                feat = [feat; sc_maskfeatcalc(maski, dnaim)];
            end
        end 
	segimg = uint16(segimg);        
        return;
    end
end

dnaim = double(imread(dnaImagePath));    
dnaim(dnaim>65000) = min(dnaim(:));	% remove pixels with very large intensity
dnaim = max(0, dnaim - median(dnaim(:))); % assumes <= half image      

im = double(imread(dicImagePath));   % im: dicim
% will be used to calculate dic content per pixel after background substraction
% the artifect has low dic content per pixel
dicbgsub = ml_preprocess_scale(im, [], 'ml', 'quantile','rc',1.4);     

im =  (im - min(min(im))) / max(max(im)) - min(min(im)); 
segimg = zeros(size(im));;

% background subtraction
bgsub = max(0,abs(im-median(im(:)))-.05);
bgval = median(im(:))+.05;

% filtered and thresholded image
rad = 12;
[xxx,yyy]=meshgrid(0:2*rad);
circ = ((xxx-rad).^2+(yyy-rad).^2 <= rad^2);
circ = circ / sum(circ(:));
bwim = conv2(bgsub, circ, 'same')>0.03;

bwim = bwmorph(bwim, 'dilate');
[L , ncomponents]= bwlabel(bwim);

stats = regionprops(L, 'area');
idx = find([stats.Area] <= 25);
L(ismember(L, idx)) =0;

reducedn = max(max(L));
idxim = [];

subregcount = 1;
for (j = 1:reducedn)    
    if (any(L(:)==j) )
        mask = L; %save intermed result as temp var!
        cond = ( (L~=j) );
        mask(cond) = 0;
        
        %we want gray original mask , not the bwmask
        newmask = im;
        newmask(find(mask ~= j)) = 0;
        
        [xx , yy] = find(newmask ~=0);
        
        [w h] = size(im);
        
        if min(xx) > 6
            ulx = min(xx) -5;
        else 
            ulx = min(xx);
        end        
        
        if min(yy) > 6
            uly = min(yy) -5;
        else 
            uly = min(yy);
        end        
        
        if max(xx) + 6 < w
            lrx = max(xx) +5;
        else 
            lrx = max(xx);
        end
        
        if max(yy) + 6 < h
            lry = max(yy) +5;
        else 
            lry = max(yy);
        end
        
        rx = uly:lry;
        ry = ulx:lrx;       
        
        if length(rx)*length(ry)>250*400
            rx = rx(1:min(250,length(rx)));
            ry = ry(1:min(400,length(ry)));            
        end
        
        dicim = double(imread(dicImagePath));
        dicsub = dicim(ry, rx);
        edgesub = double(ml_mmthin(tz_imstdedge(dicsub)>0));            
        
        fprintf(1, '%dth subregion\n',subregcount);
        
        if method.id == 1
            % spectral clustering                                
            dnasub = dnaim (ry,rx);
            idxim = segment(edgesub, dnasub, method.seed);            
        else
            % loopy
            dnaimmask = dnaim .* (newmask~=0);
            dnasub = dnaimmask (ry,rx);      
            
            if sum(dnasub)<1 
                continue;
            end
            idxim = sc_subregion(edgesub, dnasub, method, subregcount, controlDirectory);
        end
        
        if idxim == 0
            continue;
        end        
        
        saveim = zeros(size(im));   saveim(ry,rx) = idxim;         
        figure(1);  imshow(saveim, []);   saveas(1,[controlDirectory filesep 'lbl' num2str(subregcount)], 'png');
        
        if (subregcount>1)
            previdx = max(max(segimg));
            idxim(find(idxim ~=0)) = idxim(find(idxim ~=0)) + previdx;  % 1 + previdx
        end        
        maskregion = zeros(size(im));        
        maskregion(ry,rx) = idxim;
        
        segimg = segimg + maskregion;
        subregcount = subregcount+1;
    end
end

feat = [];
for i=1:max(segimg(:))
     maski = (segimg==i);
     if sum(sum(dicbgsub.*(maski)))/sum(maski(:)) < 0.01
         segimg(find(segimg==i)) = 0;
     else
         feat = [feat; sc_maskfeatcalc(maski, dnaim)];
     end
end

% find nonempty classes (excluding background)
nonempty = unique(segimg);
if (nonempty(1) == 0)
     nonempty(1) = [];
end

% renumber nonempty classes: newnum is an array whose i-th element is
% the new number for the i-th class (or 0 if the i-th class was empty)
newnum = zeros(max(nonempty),1);
newnum(nonempty) = 1:length(nonempty);

% remove empty classes
segimg(segimg>0) = newnum(segimg(segimg>0));

segimg = uint16(segimg);

max(segimg(:))
