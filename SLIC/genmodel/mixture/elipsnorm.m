function y = elipsnorm(m,cov,level,dashed);
%
% draws one contour plot of a bivariate
% gaussian density with mean "m" and covariance
% matrix "cov". 
% The level is controled by "level".
% If "dashed==1", the line is dashed.
% 
if nargin<4 
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

   dashed=0;
end
[uu,ei,vv]=svd(cov);
a = sqrt(ei(1,1)*level*level);
b = sqrt(ei(2,2)*level*level);
theta = [0:0.01:2*pi];
xx = a*cos(theta);
yy = b*sin(theta);
cord = [xx' yy']';
cord = uu*cord;
if dashed==1
   plot(cord(1,:)+m(1),cord(2,:)+m(2),'--k','LineWidth',2)
else
   plot(cord(1,:)+m(1),cord(2,:)+m(2),'k','LineWidth',2)
end
