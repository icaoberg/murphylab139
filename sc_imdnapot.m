function pots = sc_imdnapot(dnaimg,edgeimg,ncell,theta)
%SC_IMDNAPOT calculates the dna potential in the yeast image
%   POTS = SC_IMDNAPOT(DNAIMAGE,EDGEIMG,NCELL,THETA)
%   Tune the background potentail (the first entity of the potential for
%   each pixel) so that the at the meidan intensity of the edge is 0.263
%   and 95th percentile of the intensity is 0.11
%   POTS = SC_IMDNAPOT(DNAIMG,EDGEIMG,NCELL,THETA)
% 
% Input:
%   DNAIMG is the dna image
%   EDGEIMG is the smoothed edge image
%   NCELL is the number of foreground classes
%   THETA is the parameter 
%
% Output:
%   POTS is a 3-D matrix with dimension nclass X row X col
%
% usage: 
%   theta = 0.1;
%   dnapot = sc_imdnapot(dnaimg,edgeimg,ncell,theta);
%
%   See also sc_paramlogistic, sc_beliefprop, tz_imnbdens, tz_sigmoid
%
%   9-Aug-2006 Initial written by S. C. Chen

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

if nargin < 4
    error('Exactly 4 arguments are required')
end

nclass = ncell + 1;

a = theta;   % the low background potential 
b = 0.5;     % the median background potential

%targetpot = isdna() * a + (1-isdna()) * b;
%targetpot = isdna() * (a - b) + b;
%isdna() = sigmoid( dnaslope * ( dnaimg - dnathres ) );
%targetpot = sigmoid( dnaslope * ( dnaimg - dnathres ) ) * (a - b) + b;
%sigout = sigmoid( dnaslope * ( dnaimg - dnathres ) );

targetpot = [0.11;1/3-0.07]; % background pot, 0.11 if 95th, 0.263 if at the meidan of the edge
sigout = (targetpot - b)./(a - b);

% use a gaussian filter to smooth the dna image
dnafilter = fspecial('gaussian',5,2);
dnaimg = imfilter(dnaimg,dnafilter,'symmetric');
dnaimg = uint8(round(mat2gray(double(dnaimg))*255));
% calculate the slope and threshold of the logistic regression
pth = 0.95;
sortdna = sort(dnaimg(:));          
xi = double(sortdna(floor(pth*length(sortdna))));    
xi = [xi; median(double(dnaimg(find(edgeimg>0))))];
bi = log((1-sigout)./sigout);
if xi(1) == xi(2)
     pots = 0;
     return;
end

dnaslope = (bi(2) - bi(1))/(xi(1) - xi(2));
dnathres = bi(1)/dnaslope + xi(1);

if dnathres == 0
    pots = 0;
    return;
end

isdna = tz_sigmoid(double(dnaimg),dnathres,dnaslope);
pots = zeros(nclass,size(dnaimg,1),size(dnaimg,2));
% calculate the background potential value
pots(1,:,:) = isdna * a + (1 - isdna) * b;

% foreground is calculated as (1-background)/ncell;
for k=2:nclass
    pots(k,:,:) = (1-pots(1,:,:))/ncell;
end
