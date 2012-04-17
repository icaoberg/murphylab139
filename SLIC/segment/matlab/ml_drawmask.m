function ml_drawmask(img)
%ML_DRAWMASK Draw masks on an image
%   ML_DRAWMASK(IMG) lets a user draw polygon on the image IMG by clicking 
%   left mouse button. Click the right button or double click the left 
%   button to finish one polygon. Read the title of the figure to get hints.
%
%   IMPORTANT: If you finds any bad English in this program, please fix it. 
%   Thanks.

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

%   09-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

while(1)  
    img = double(img);
    imshow(img,[]);
    title('Draw a polygon ...');
    bw = roipoly;%(mat2gray(img));
    bw = bwperim(bw);
    img(find(bw==1)) = max(img(:));
    imshow(img,[]);  
    title('click LEFT BUTTON to continue or RIGHT BUTTON to stop');
    [x,y,button] = ginput(1);  
    if(button == 3)
        title('Drawing done.');
        break; 
    end
end
