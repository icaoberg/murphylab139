function [feats,names,slfnames] = sc_regionfeatmat(img,dnaimg,masks,featset,har_intbins,downsamplerate)
%SC_REGIONFEATMAT Calculate feature matrix of regions in an image.
%   FEATS = SC_REGIONFEATMAT(IMG,DNAIMG,MASKS) returns a feature matrix of 
%   the regions in IMG. DNAIMG is an image from dna channel and used in 
%   feature calculation. See TZ_MASKIMGS for more information about the two
%   parameters IMG and MASKS.
%   
%   [FEATS,NAMES,SLFNAMES] = TZ_REGIONFEATMAT(...) also returns feature
%   names and SLF names.
%
%   See also

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

%   11-Nov-2005 Initial write T. Zhao
%   14-Mar-2007 Modified by S.-C. Chen
%   Copyright (c) Center for Bioimage Informatics, CMU

if masks == -1
    % Calculate field level features
    regions{1} = img;   dnas{1} = dnaimg;
else
    regions = tz_maskimgs(img,masks);

    if ~isempty(dnaimg)
        dnas = tz_maskimgs(dnaimg,masks);
    end
end

feats = [];
names = [];
slfnames = [];

for k=1:length(regions)
    if ~isempty(dnas)
        dnaRegion = dnas{k};
    else
        dnaRegion = [];
    end
    
    % scale = 0.23;
    % scale_factor = scale / default_scale;   scale_factor = 1;
    % radius_pixel = radius / scale; => input radius = radius_pixel * scale = 15 * 0.23
    % 1 / downsamplerate = scale / har_pixsize; =>  har_pixsize = scale(0.23) * downsamplerate 

    scale = 0.23;
    har_pixsize = scale * downsamplerate;
    radius =  15 * scale;  % the radius of yeast is 15 pixels
    [names, feats(k,:), slfnames] = ml_featset( double(regions{k}), [], double(dnaRegion),featset,scale,radius,'yesbgsub','rc', har_pixsize, har_intbins);  
end
