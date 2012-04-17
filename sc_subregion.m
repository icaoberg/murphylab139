function idxim = sc_subregion(edgeimg, dnaimg, method, subregcount, controlDirectory)
%SC_SUBREGION sets up the edge and dna potential for loopy belief
%propagation with voting potential. It saves intermediate results in the
%control directory and returns the index image of the masks
%   IDXIM = SC_SUBREGION(EDGEIMG, DNAIMG, METHOD, SUBREGCOUNT, CONTROLDIRECTORY)
%   
% Input:
%   EDGEIMG is the smoothed edge image in a subregion
%   DNAIMG is the dna image in a subregion
%   METHOD is a struct which has the information the type and parameters of
%   three segmentation methods.  Please see sc_segyeast for details
%   SUBREGCOUNT is the index of the subregion in this image.  It is used
%   for saving the intermediate results in the control directory.
%   CONTROLDIRECTORY is the path for the control directory
%   The control directory is used to indicate if the image has been
%   processing and also store the intermediate segmentation results
%
% Output:
%   SEGIMG is the index image for the masks
%
% usage:
%   idxim = sc_subregion(edgesub, dnasub, method, subregcount, controlDirectory)
%
%   See also: sc_cutcells, sc_runloopy, tz_imnbdens, sc_imdnapot
%
%   09-Aug-2006 Written by S. C. Chen

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

% number of foreground classes per pixel
ncell = 2;  
% number of classes per pixel (ncell + 1 background class)
nclass = ncell + 1; 
% W, the neighborhood window size;  each pixel has W X W - 1 neighbors
wndsz = 7;  
% the size of the image
[row, col] = size(dnaimg);  
% the returned index image; intialized to 0
idxim = 0; 

% use a 7 X 7 gaussian window to smooth the edge image
gauss = -3:1:3; gauss = exp(-gauss.^2./2); gauss = gauss/sum(gauss);
edgeimg = min(.8, conv2(edgeimg, gauss'*gauss, 'same'));
% calculate the edge potential with tz_imnbdens.m
% edgepot is a 4-D matrix with dimension row X col X W X W. 
edgepot = tz_imnbdens(edgeimg,(wndsz-1)/2);

% normalize the edgepot using a sigmoid function
% the threshold and slope are empirically determined
edgethres = 1;    edgeslope = 5;
edgepot = tz_sigmoid(edgepot,edgethres,edgeslope);

% use the dnaimg, smoothed edgeimg, ncell and theta to determine the dna potential
% dnapot is a 3-D matrix with dimension nclass X row X col
theta = 0.1;
dnapot = sc_imdnapot(dnaimg,edgeimg,ncell,theta);
if dnapot == 0
    % sometimes the subregion does not have any DNA intensity; ignore this region
    fprintf(1, 'no dna content in this subregion\n\n');
    return;
end

if method.imgidx == 1 & subregcount == 11
    % this is to show the Fig. 4 in the CBICB paper:  An example of a
    % boundary contour superimiposed on the DNA potential of the foreground label. 
    figure(1); imshow(1.3*mat2gray(double(edgeimg)) + mat2gray(squeeze(dnapot(2,:,:))), []);     
    print('-depsc2', 'figures/Fig4.eps');
end
% %%%%%%%%%%%% generate a colored image to check if the dna and edge potential are generated correctly
% figure(1); imshow(1.3*mat2gray(double(edgeimg)) + mat2gray(squeeze(dnapot(2,:,:))), []);
% while 1    
%     [y, x] = ginput(1);
%     if y<0
%         break;
%     end
%     fprintf(1, 'click x = %d, y = %d\n', round(y), round(x));    
%     fprintf(1, 'dnapot:');
%     dnapot(:, round(x), round(y))
%     fprintf(1, 'edgepot:');
%     squeeze(edgepot(round(x),round(y),:,:)) 
% end

% set the parameters for sc_runloopy
param.ncell = ncell;
param.alpha = method.alpha;  
param.beta = method.beta;
param.wndsz = wndsz;
param.row = row;
param.col = col;
param.stopthreshold = 0.1;

% run loopy belief propagation
belief = sc_runloopy(dnapot, edgepot, param);

NumCell = 0;
% find which of the foreground class has the largest belief
[tmpy, labels] = max(belief(2:nclass,:,:));
labels = squeeze(labels);
labels = labels + 1;    % class index starting from 2
% if the background class has the largest belief (saved in bgidx), set the
% labels of the background pixels to be 1
bgidx = find(squeeze(belief(1,:,:)) >= 1/nclass);     labels(bgidx) = 1;    

% firstiteration is a flag to identify if the following loopy is at its first iteration
firstiteration = 1;
while 1
    % First we want to find what is the pixel that is most likely to be at
    % the foreground; i.e. its background potential is lowest
    % find the backgroundbelief for each pixel
    backgroundbelief = squeeze(belief(1,:,:));
      
    if firstiteration == 0
        % if not at the first iteration, set the backgroundbelief of the
        % pixels already classified as the background to be 1 (highest) 
        backgroundbelief(find(labels==1)) = 1;
    end
    
    % sort the backgrondbelief in ascending order
    [bmin,index] = sort(backgroundbelief(:));   
    % find [x,y], the cooridinates of the pixel with the lowest backgroundbelief 
    x = mod(index(1),size(dnaimg,1));
    y = (index(1) - x)/size(dnaimg,1);
    x = x + 1;    y = y + 1;    % let the minimal possible cooridinate is (1,1)
    
    % check if [x,y] is out of the boundary (or very close to)
    if x > size(belief,2)-2 | y > size(belief,3)-2 | x < 3 | y < 3
        fprintf(1, 'Auto-clicked coordinate is out of bound...\n');    
        fprintf(1, 'Click (%d,%d) in image of size [%dx%d]\n', x, y, size(belief,2), size(belief,3)');                
        return;
    end
    
    if firstiteration
        % it is at the first iteration, set firstiteration to 0
        firstiteration = 0;
        figure(1);  imshow(labels, []);
        % the bgidx are the pixels with labels = 1
        bgidx = find(labels == 1);    
        dnapotrs = reshape(dnapot, nclass, row*col);
        bgvec = [1; zeros(ncell,1)];
        
        % set the dnapot of the bgidx pixels to be [1;0;0]
        dnapotrs(:, bgidx) = repmat(bgvec, [1 length(bgidx)]);
        dnapot = reshape(dnapotrs, nclass, row, col);
        
        fgidx = find(labels >1);
        if isempty(fgidx)
            NumCell = 0;
            fprintf('There are no cells in the image\n');
            break;
        else
            NumCell = 1;
        end
    else        
        NumCell = NumCell + 1;
    end
    
    % set the dnapot of the auto-clicked pixel and its 8 neighbors to be [0;0;1]    
    dnapot(:, round(x),round(y-1))   = [zeros(ncell,1);1];
    dnapot(:, round(x),round(y))     = [zeros(ncell,1);1];
    dnapot(:, round(x),round(y+1))   = [zeros(ncell,1);1];
    dnapot(:, round(x-1),round(y-1)) = [zeros(ncell,1);1];
    dnapot(:, round(x-1),round(y))   = [zeros(ncell,1);1];
    dnapot(:, round(x-1),round(y+1)) = [zeros(ncell,1);1];
    dnapot(:, round(x+1),round(y-1)) = [zeros(ncell,1);1];
    dnapot(:, round(x+1),round(y))   = [zeros(ncell,1);1];
    dnapot(:, round(x+1),round(y+1)) = [zeros(ncell,1);1];
        
    % run belief propagation again with new dna potential
    belief = sc_runloopy(dnapot, edgepot, param);     
    
    % show the labeled images
    [tmpy, labels] = max(belief);
    labels = squeeze(labels);
    figure(2); imagesc(labels);
    
    % find those pixels with the label == nclass
    L = bwlabel(labels==nclass);   
    lblidx = find(L == L(round(x),round(y)));   
    % To avoid the problem of the effect of the numerical precision, ignore the pixels which does not have
    % significant 2nd foreground belief over the 1nd foreground belief
    lblidx(find(belief(nclass,lblidx) - belief(nclass-1,lblidx)<1e-10)) = [];
    
    if length(lblidx) > 150
        fprintf(1, '========%d mask is generated========\n', NumCell);        
    else
        fprintf(1, 'The segmented mask is too small; \n %d masks are generated and the segmentation stops...', NumCell-1);
        if NumCell == 1
            return;
        else
            break;
        end
    end   
    
    % diliate this cell a little bit
    segcell = zeros(size(labels));    
    segcell(lblidx) = 1;
    se = strel('ball',4,4);    
    segcell = double(imdilate(segcell,se)>4.2);
    
    % store the cells at this subregion
    segcells{NumCell} = segcell;    
    
    % save intermediate results in the control directory
    imwrite(segcell, [controlDirectory filesep 'sub' num2str(subregcount) '_cell' num2str(NumCell) '.png']);                   
    
    dnapotrs = reshape(dnapot, nclass, row*col);
    dnapotrs(:, lblidx) = repmat(bgvec, [1 length(lblidx)]);  % set all pixels of this label to be background
    dnapot = reshape(dnapotrs, nclass, row, col);
    
    labels(lblidx) = 1;   
    rgbimg(:,:,1)=mat2gray(labels);
    rgbimg(:,:,2)=mat2gray(double(edgeimg));
    rgbimg(:,:,3)=mat2gray(squeeze(dnapot(2,:,:)));
    % save intermediate results in the control directory    
    imwrite(rgbimg, [controlDirectory filesep 'sub' num2str(subregcount) '_totcell' num2str(NumCell) '.png']);                     
end

if NumCell == 0
    return;
end

idxim = zeros(size(segcells{1}));

for j = 1:size(segcells,2)
    idxim(segcells{j}>.9) = j;
end
idxim = reshape(idxim, size(edgeimg));

class = idxim;

% remove masks that don't pass heuristic tests
for j = 1:max(class(:))
    classok = sc_testmask(class==j, dnaimg);
    if (classok<0)
        class(class==j) = 0;
    end
end

% find nonempty classes (excluding background)
nonempty = unique(class(:));
if (nonempty(1) == 0)
    nonempty = nonempty(2:end);
end

% renumber nonempty classes: newnum is an array whose i-th element is
% the new number for the i-th class (or 0 if the i-th class was empty)
newnum = zeros(max(nonempty),1);
newnum(nonempty) = 1:length(nonempty);

% remove empty classes
class(class>0) = newnum(class(class>0));
idxim = class;
