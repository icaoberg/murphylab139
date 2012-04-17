function img2=ml_draw_line2(r1,c1,r2,c2,img,lineVal)
% function img2=ml_draw_line2(r1,c1,r2,c2,img,lineVal)
% This function draws a line on img starting from point
% [r1,c1] to point [r2,c2] with intensity lineVal

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

if(isempty(img))
     img2=[];
     return;
end

pts=ml_getlinepts(round([r1 c1]),round([r2 c2]));            

pts=[pts,lineVal+zeros(size(pts,1),1)];
ng=sum(pts<=0,2)+(pts(:,1)>size(img,1))+(pts(:,2)>size(img,2));
pts(ng>0,:)=[];

if(~isempty(pts))
    img2 = ml_setimgptspixel(img,pts);
else
    img2=img;
end
