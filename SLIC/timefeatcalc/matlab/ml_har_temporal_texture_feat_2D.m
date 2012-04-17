function [feat,tnames] = ml_har_temporal_texture_feat_2D (dir_name,mask_name,pattern,level,graylevel,har_pixsize,scale)

%ML_HAR_TEMPORAL_TEXTURE_FEAT_2D calculate 2D haralick temporal texture features
%   FEAT = ML_HAR_TEMPORAL_TEXTURE_FEAT_2D(DIR_NAME,MASK_NAME,PATTERN,LEVEL,GRAYLEVEL,HAR_PIXSIZE,SCALE)
%DIR_NAME: directory name where one movie is stored
%MASK_NAME: file name for the one mask of all the images in the movie
%PATTERN: regular expression for the 2D image file name
%LEVEL: how many different time intervals to consider
%GRAYLEVEL: default is 255
%HAR_PIXSIZE: image will be rescaled to har_pixsize per pixel,
%             default is 1.1 micron/pixel
%SCALE: the real resolution of the image, 0.11 is default for spinning disk we use

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

%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%  
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%  
%  Yanhua Hu, Dec 2005. The version now is to use the same background for all the images and no other preprocessing like thresholding. 



DEFAULT_LEVEL = 5;
DEFAULT_INTENSITY_BINS = 255;
DEFAULT_HAR_PIXEL_SIZE = 1.1;
DEFAULT_SCALE = 0.11;   % micron/pixel

if ~exist('level', 'var')
     level=DEFAULT_LEVEL;
end
if ~exist('graylevel','var')
     graylevel=DEFAULT_INTENSITY_BINS;
end
if ~exist('har_pixsize', 'var')
     har_pixsize=DEFAULT_HAR_PIXEL_SIZE;
end
if ~exist('scale', 'var')
     scale=DEFAULT_SCALE;
end

har_rescale_factor = scale / har_pixsize;
tnames = {'angular_second_moment' 'contrast' 'correlation' 'sum_of_squares' 'inverse_diff_moment' 'sum_avg' 'sum_var' 'sum_entropy' 'entropy' 'diff_var' 'diff_entropy' 'info_measure_corr_1' 'info_measure_corr_2'};

d=ml_dir(fullfile(dir_name,pattern));
img_first=double(ml_readimage(fullfile(dir_name, d{1})));
if isempty(mask_name)
   mask_img = ones(size(img_first));
else
   mask_img = double(imread(mask_name));
end

% find out how many time points
num_time = length(d);

feat=[]; 

common = ml_imgcommonpixel(img_first) ;
img_first = mask_img.*(img_first-common); 

% find the largest pixel value in the all the time point
Max_pixel= max(img_first(:));
Rescale_factor=graylevel/Max_pixel;

for k=1:level % different time interval
 feat_vals=[];
 for i=1:num_time
   if i+k<=num_time
      img1=double(ml_readimage(fullfile(dir_name, d{i})));
      img1 = max(0,mask_img.*( img1-common));
      img2=double(ml_readimage(fullfile(dir_name, d{i+k})));
      img2 = max(0,mask_img.*( img2-common));
      if (~isempty(img1) & ~isempty(img2))
         %%%% in mb_imgscle = 
         %scaledimage = image - min(image(:)) ;
         %scaledimage = scaledimage/max(scaledimage(:)) ;
         img1 = imresize(img1, har_rescale_factor, 'bilinear');
         img1= uint8(min(graylevel,img1*Rescale_factor));
         img2 = imresize(img2, har_rescale_factor, 'bilinear'); 
         img2= uint8(min(graylevel,img2*Rescale_factor));
   
         cocmat = ml_updtCOCMAT(img1,img2,[],graylevel);
     
         % make it sparse
         h=1;
         while h<size(cocmat,1)
           if sum(cocmat(h,:))==0
            cocmat(h,:)=[];
            cocmat(:,h)=[];
            h=h-1;
           end
           h=h+1;
         end

         cocmat=single(cocmat/sum(cocmat(:)));
         [row,col]=size(cocmat);
         if not(isempty(cocmat)|row==1|col==1)
           vals=ml_Har_Temporal_Texture(cocmat);
           vals=(vals(1:13))';
           feat_vals=[feat_vals;vals]; % for the same time interval
         else 
           warn=d{i}
         end
     end 
   end
 end

 if ~isempty(feat_vals)
  feat_valsize = size(feat_vals)
  feat_vals=double(feat_vals);
  feat_mean= mean(feat_vals,1);
  feat_var= var (feat_vals);
  feat = [feat,feat_mean,feat_var];
 end

end 
