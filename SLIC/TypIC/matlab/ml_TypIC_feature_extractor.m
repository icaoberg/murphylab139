function [numfeat, numfiles, imagenames, AllFeatures] = ...
    ml_TypIC_feature_extractor(dir_name, dir_crop, featsetname);

% Extracts features from a given set of images for use of TypIC
% Taken from beginning of TypIX.m by Robert Murphy
% Modified Superbowl Sunday Februaray 2, 2002 by Jennifer Lin
% Modified February 25, 2002 by Jennifer Lin to rename from TypIC_feature_extractor.m to jl_TypIC_feature_extractor.m

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

% addpath multout
%addpath /home/lin4/matlab

imagenames = ml_dir(dir_name);
% Get rid of . and ..
imagenames(1:2)=[];

if (isempty(dir_crop))
  crop_imagenames = [];
else
  crop_imagenames = ml_dir(dir_crop); 
  % Get rid of . and ..
  crop_imagenames(1:2) = [];
end

n_images = length(imagenames);

AllFeatures = [];

for N = 1 : n_images

  file = [dir_name '/' imagenames{N}]; 
  % Read the images
  im = ml_readimage(file);
 
  if ((isempty(dir_crop)) == 0)
     crop_file = [dir_crop '/' crop_imagenames{N}];
     % Read crop images
     im_crop = ml_readimage(crop_file);
  else
     im_crop = [];
  end 

  % Calculate the features
  %which ml_featset;
  [featnames, features] = ml_featset(double(im), double(im_crop), [], featsetname);

  AllFeatures = [AllFeatures; features];

end

[numfiles, numfeat] = size(AllFeatures);
