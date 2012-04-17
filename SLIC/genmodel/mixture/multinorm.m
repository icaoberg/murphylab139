function y = multinorm(x,m,covar)
% Evaluates a multidimensional Gaussian
% of mean m and covariance matrix covar
% at the array of points x
%
[dim npoints] = size(x);
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

dd = det(covar+ realmin*eye(dim));
in = inv(covar+ realmin*eye(dim));
ff = ((2*pi)^(-dim/2))*((dd)^(-0.5));
quadform = zeros(1,npoints);
centered = (x-m*ones(1,npoints));
if dim ~= 1
   y = ff * exp(-0.5*sum(centered.*(in*centered)));
else
   y = ff * exp(-0.5*in*centered.^2 );
end



