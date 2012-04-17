function ml_make_sasstepdiscinp( feature_matrix, fname)
% function ML_MAKE_SASSTEPDISCINP( FEATURE_MATRIX, FNAME)
% This function is an OBSOLETE SDA code.  It loads the features 
% from the FEATURE_MATRIX name, standardizes the features, 
% and save the matrix in comma separated text file with FNAME

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

% Copyright (C) 2006  Murphy Lab
% Carnegie Mellon University  


% load features if necessary
if( ischar(feature_matrix))
    f=load(feature_matrix);
    feature_matrix = f.all_features;
end

% Standardize features
if( iscell(feature_matrix))
    feats = [];
    classes = [];
    for cls = 1:length(feature_matrix)
	cls_feats = double(feature_matrix{cls});
	feats = [feats; cls_feats];
	classes = [classes; repmat(cls,[size(cls_feats,1) 1])];
    end
else
    feats = double(feature_matrix(:,1:end-1));
    classes = double(feature_matrix(:,end));
end
MEAN = mean(feats);
STD = std(feats);
MEAN_MATRIX = repmat(MEAN,[size(feats,1) 1]);
feats = feats - MEAN_MATRIX;
STD_MATRIX = repmat(STD,[size(feats,1) 1]);
feats = feats ./ STD_MATRIX;
% Set z-scored features to zero where std was zero
feats(find(STD_MATRIX==0))=0;

std_feature_matrix = [classes feats];
% Save them in comma separated text file
ml_dlmwrite(fname,std_feature_matrix,',');

