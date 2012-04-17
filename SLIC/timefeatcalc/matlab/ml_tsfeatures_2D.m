function feat_vals = ml_tsfeatures_2D(dir_name,img_pattern,featsets,imagemask,scale)

% ML_TSFEATURES_2D Calculate time series features for 2D movies
% FEAT_VALS = ML_TSFEATURES_2D(DIR_NAME,IMG_PATTERN,FEATSETS,IMAGEMASK,SCALE)
% DIR_NAME is the directory where all images are stored
% IMG_PATTERN is the file pattern of images
% FEATSETS identifies the sets of features to be calculated. There
% is a three letter code for each feature set and they are calculated 
% in the order specified. 
%     har - haralick temporal texture features (the only set for now)
% IMAGEMASK is the file name for the mask
% SCALE is distance in mircometers that each pixel represent. 
%        By default it is 0.11um/pixel. 
% Omitting the SCALE and IMAGEMASK arguments will have the same effect as
% entering them as 0 or [].

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

% By Yanhua Hu, July 2005. Code format follows the ml_features.m as the previous work in the lab.

% Check arguments & use defaults where necessary


DEFAULT_SCALE = 0.11;   % micron/pixel

if( nargin < 5| isempty(scale)|scale == 0 )
        scale = DEFAULT_SCALE;
end
if( nargin < 4 | imagemask==0)
	imagemask = [];
end
if( nargin < 3)
	error('This function takes at least 3 arguments');
end

if ~exist('dir_name','var')
	error('You must supply a directory path where the movie is stored');
end

if ~exist('img_pattern','var')
	error('You must supply the regular expression of the image file names');
end

if ~exist('featsets', 'var')
        error('You must choose some feature sets to calculate');
end


  
  % init output variables
  %feat_slf = [];
  %feat_names = [];
   feat_vals = [];

for i = 1 : length( featsets)
  switch featsets{i}
    case 'har' % haralick Temporal Texture features
       DEFAULT_LEVEL = 5;
       DEFAULT_INTENSITY_BINS = 255;
       DEFAULT_HAR_PIXEL_SIZE = 1.1;
       [tvalues,tnames] = ml_har_temporal_texture_feat_2D (dir_name,imagemask,img_pattern,DEFAULT_LEVEL,DEFAULT_INTENSITY_BINS,DEFAULT_HAR_PIXEL_SIZE,scale);
       feat_vals = tvalues;
  end
end
